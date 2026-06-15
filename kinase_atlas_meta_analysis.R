


#!/usr/bin/env Rscript

## ===============================================================
## ESR1 heterogeneity + kinase correlation + kinase driver pipeline
##
## Core steps:
##   01) ESR1 QC / heterogeneity per atlas sample
##   02) Per-sample metacell ESR1~kinase-signature correlations
##   03) Weighted meta-analysis of kinase correlations
##   04) Per-sample gene-level driver analysis within selected kinases
##   05) Weighted meta-driver summary across samples
##
## Seurat v5 compatible: uses GetAssayData(..., layer=)
## ===============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(Matrix)
  library(UCell)
  library(BiocParallel)
  library(FNN)
  library(stringr)
  library(qs)
  library(scales)
  library(tools)
})

## ---------------------------------------------------------------
## CONFIG
## ---------------------------------------------------------------
CFG <- list(
  sample_list_csv = "/path/to/sample_list_cancer_cells_102_samples.csv",
  kinase_signature_file = "/path/to/ARCHS4_Kinases_Coexp.txt",
  out_dir = "/path/to/out_dir",
  
  assay_preference = c("SCT", "RNA"),
  reduction_for_metacells = "umap",
  cells_per_metacell = 25L,
  min_genes_per_signature = 5L,
  n_workers_ucell = 4L,
  kinase_suffix_pattern = "kinase ARCHS4 coexpression$",
  
  qc_min_cells = 250L,
  qc_max_esr1_frac_zero = 0.95,
  qc_min_esr1_sd = 0.05,
  qc_esr1_mean_clip = 2,
  
  meta_min_samples_per_signature = 3L,
  meta_hardness = 1.0,
  
  K_SIGN = 100L,
  driver_alpha = 1.0
)

## ---------------------------------------------------------------
## DIRS / FILES
## ---------------------------------------------------------------
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CFG$out_dir, "per_sample_corr"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CFG$out_dir, "per_sample_meta"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CFG$out_dir, "per_sample_drivers"), recursive = TRUE, showWarnings = FALSE)

FILES <- list(
  log = file.path(CFG$out_dir, "esr1_kinase_pipeline.log"),
  qc_summary = file.path(CFG$out_dir, "cancer_samples_ESR1_QC_summary.csv"),
  qc_weights = file.path(CFG$out_dir, "cancer_samples_ESR1_QC_with_weights.csv"),
  corr_all = file.path(CFG$out_dir, "ESR1_kinase_corr_all_samples.csv"),
  sig_cov = file.path(CFG$out_dir, "kinase_signature_coverage_per_sample.csv"),
  corr_meta = file.path(CFG$out_dir, "ESR1_kinase_corr_meta_weighted.csv"),
  meta_drivers = file.path(CFG$out_dir, "ESR1_kinase_meta_drivers_gene_scores.csv")
)

log_con <- file(FILES$log, open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE)

log_msg <- function(..., .sep = " ") {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = .sep))
  writeLines(msg, con = log_con)
  message(msg)
}

log_msg("==== Starting ESR1 heterogeneity / kinase pipeline ====")

## ---------------------------------------------------------------
## HELPERS
## ---------------------------------------------------------------

safe_layer <- function(object, assay = NULL, layer_name) {
  tryCatch(GetAssayData(object, assay = assay, layer = layer_name), error = function(e) NULL)
}

choose_expression_assay <- function(so, preferred = c("SCT", "RNA")) {
  available <- names(so@assays)
  for (a in preferred) {
    if (a %in% available) return(a)
  }
  if (length(available) == 0) stop("No assays found in Seurat object")
  available[[1]]
}

ensure_umap_reduction <- function(so, reduction_name = "umap", dims = 1:30) {
  if (reduction_name %in% names(so@reductions)) return(so)
  log_msg("  UMAP not found; attempting to compute it")
  if (!("pca" %in% names(so@reductions))) {
    so <- RunPCA(so, verbose = FALSE)
  }
  RunUMAP(
    so,
    reduction = "pca",
    dims = dims,
    reduction.name = reduction_name,
    verbose = FALSE
  )
}

normalize_sig_name <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "_UCell$", "")
  x <- str_replace(x, " human kinase ARCHS4 coexpression$", "")
  x <- str_replace(x, " human kinase ARCHS4 coexpression_UCell$", "")
  x <- str_replace(x, "^Kinase_", "")
  x <- str_replace(x, "^ARCHS4_Kinase_", "")
  x <- gsub("[^A-Za-z0-9]+", "", x)
  toupper(x)
}

safe_cor <- function(x, y, method = "pearson") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || length(y) < 3L) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

safe_cortest_p <- function(x, y, method = "pearson") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || length(y) < 3L) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  out <- suppressWarnings(try(cor.test(x, y, method = method), silent = TRUE))
  if (inherits(out, "try-error")) return(NA_real_)
  out$p.value
}

load_seurat_any <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  ext <- tolower(file_ext(path))
  obj <- if (ext == "qs") qs::qread(path) else readRDS(path)
  if (!inherits(obj, "Seurat")) stop("Loaded object is not Seurat: ", paste(class(obj), collapse = ";"))
  obj
}

make_umap_metacells_scRNA <- function(seurat_object,
                                      reduction = "umap",
                                      cells_per_metacell = 25L) {
  if (!(reduction %in% names(seurat_object@reductions))) {
    stop("Reduction '", reduction, "' not found")
  }
  emb <- Embeddings(seurat_object, reduction)
  all_cells <- rownames(emb)
  unassigned <- all_cells
  metacell_id <- setNames(rep(NA_character_, length(all_cells)), all_cells)
  mc_index <- 1L
  
  while (length(unassigned) > 0L) {
    seed <- sample(unassigned, 1)
    k <- min(cells_per_metacell, length(unassigned))
    nn <- FNN::get.knnx(
      data = emb[unassigned, , drop = FALSE],
      query = emb[seed, , drop = FALSE],
      k = k
    )
    grp <- unassigned[nn$nn.index[1, ]]
    metacell_id[grp] <- paste0("MC", mc_index)
    unassigned <- setdiff(unassigned, grp)
    mc_index <- mc_index + 1L
  }
  
  seurat_object$metacell_id <- metacell_id[Cells(seurat_object)]
  seurat_object
}

meta_one_signature <- function(df_sig, min_samples_per_signature) {
  k <- nrow(df_sig)
  if (k < min_samples_per_signature) {
    return(tibble(
      n_samples = k,
      n_samples_eff = NA_integer_,
      r_meta = NA_real_,
      z_meta = NA_real_,
      p_meta = NA_real_,
      sum_w_eff = sum(df_sig$w_eff),
      mean_r = mean(df_sig$r),
      sd_r = ifelse(k > 1, sd(df_sig$r), NA_real_),
      n_pos = sum(df_sig$r > 0),
      n_neg = sum(df_sig$r < 0)
    ))
  }
  z_i <- atanh(df_sig$r)
  w_i <- df_sig$w_eff
  sum_w <- sum(w_i)
  z_bar <- sum(w_i * z_i) / sum_w
  r_meta <- tanh(z_bar)
  z_meta <- z_bar * sqrt(sum_w)
  p_meta <- 2 * pnorm(-abs(z_meta))
  tibble(
    n_samples = k,
    n_samples_eff = k,
    r_meta = r_meta,
    z_meta = z_meta,
    p_meta = p_meta,
    sum_w_eff = sum_w,
    mean_r = mean(df_sig$r),
    sd_r = ifelse(k > 1, sd(df_sig$r), NA_real_),
    n_pos = sum(df_sig$r > 0),
    n_neg = sum(df_sig$r < 0)
  )
}

## ---------------------------------------------------------------
## LOAD SAMPLE LIST
## ---------------------------------------------------------------
log_msg("Reading sample list:", CFG$sample_list_csv)
samples <- read_csv(CFG$sample_list_csv, show_col_types = FALSE) %>%
  mutate(
    Path_raw = Path,
    rds_path = gsub('(^"|"$)', "", Path_raw),
    file_ext = tools::file_ext(rds_path),
    loader = if_else(tolower(file_ext) == "qs", "qs", "rds")
  )
log_msg("Loaded", nrow(samples), "samples")

## ---------------------------------------------------------------
## LOAD KINASE SIGNATURES
## ---------------------------------------------------------------
log_msg("Reading kinase signatures:", CFG$kinase_signature_file)
kinase_signatures_raw <- read.table(
  CFG$kinase_signature_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)
kinase_names <- as.character(kinase_signatures_raw[, 1])
kinase_signatures_all <- lapply(seq_len(nrow(kinase_signatures_raw)), function(i) {
  row <- as.character(kinase_signatures_raw[i, -1])
  row <- row[!is.na(row) & row != "" & row != "NA"]
  unique(row)
})
names(kinase_signatures_all) <- kinase_names
log_msg("Loaded", length(kinase_signatures_all), "kinase signatures")

## ---------------------------------------------------------------
## STEP 01: ESR1 QC / HETEROGENEITY
## ---------------------------------------------------------------
compute_esr1_qc <- function(rds_path, sample_id, treatment = NA_character_, loader = "rds", file_ext = NA_character_) {
  fail_row <- function(msg, loaded_ok = FALSE,
                       n_cells = NA_integer_, n_features = NA_integer_, assay_used = NA_character_,
                       median_nCount_RNA = NA_real_, median_nFeature_RNA = NA_real_, median_pct_mt = NA_real_) {
    tibble(
      Sample = sample_id,
      Treatment = treatment,
      rds_path = rds_path,
      file_ext = file_ext,
      loader = loader,
      loaded_ok = loaded_ok,
      n_cells = n_cells,
      n_features = n_features,
      assay_used = assay_used,
      esr1_mean = NA_real_,
      esr1_median = NA_real_,
      esr1_sd = NA_real_,
      esr1_iqr = NA_real_,
      esr1_min = NA_real_,
      esr1_max = NA_real_,
      esr1_frac_zero = NA_real_,
      esr1_frac_na = NA_real_,
      median_nCount_RNA = median_nCount_RNA,
      median_nFeature_RNA = median_nFeature_RNA,
      median_pct_mt = median_pct_mt,
      error_msg = msg
    )
  }
  
  tryCatch({
    so <- load_seurat_any(rds_path)
    meta <- so@meta.data
    
    first_existing_col <- function(meta_df, candidates) {
      hit <- candidates[candidates %in% colnames(meta_df)]
      if (length(hit) == 0) return(NA_character_)
      hit[[1]]
    }
    safe_median_candidates <- function(meta_df, candidates) {
      col <- first_existing_col(meta_df, candidates)
      if (is.na(col)) return(NA_real_)
      median(meta_df[[col]], na.rm = TRUE)
    }
    
    median_nCount_RNA <- safe_median_candidates(meta, c("nCount_RNA", "nCount_SCT", "nCount"))
    median_nFeature_RNA <- safe_median_candidates(meta, c("nFeature_RNA", "nFeature_SCT", "nFeature"))
    median_pct_mt <- safe_median_candidates(meta, c("percent.mt", "pct_mt", "percent_mt", "pct.mt"))
    
    assay_used <- choose_expression_assay(so, CFG$assay_preference)
    DefaultAssay(so) <- assay_used
    
    data_mat <- safe_layer(so, assay = assay_used, layer_name = "data")
    counts_mat <- safe_layer(so, assay = assay_used, layer_name = "counts")
    
    esr1_vec <- NULL
    if (!is.null(data_mat) && "ESR1" %in% rownames(data_mat)) {
      esr1_vec <- data_mat["ESR1", ]
    } else if (!is.null(counts_mat) && "ESR1" %in% rownames(counts_mat)) {
      esr1_vec <- counts_mat["ESR1", ]
    }
    
    if (is.null(esr1_vec)) {
      return(fail_row(
        "ESR1_not_found",
        loaded_ok = TRUE,
        n_cells = ncol(so),
        n_features = nrow(so),
        assay_used = assay_used,
        median_nCount_RNA = median_nCount_RNA,
        median_nFeature_RNA = median_nFeature_RNA,
        median_pct_mt = median_pct_mt
      ))
    }
    
    esr1_vals <- as.numeric(esr1_vec)
    frac_zero <- mean(esr1_vals == 0, na.rm = TRUE)
    frac_na <- mean(is.na(esr1_vals))
    
    tibble(
      Sample = sample_id,
      Treatment = treatment,
      rds_path = rds_path,
      file_ext = file_ext,
      loader = loader,
      loaded_ok = TRUE,
      n_cells = ncol(so),
      n_features = nrow(so),
      assay_used = assay_used,
      esr1_mean = mean(esr1_vals, na.rm = TRUE),
      esr1_median = median(esr1_vals, na.rm = TRUE),
      esr1_sd = sd(esr1_vals, na.rm = TRUE),
      esr1_iqr = IQR(esr1_vals, na.rm = TRUE),
      esr1_min = suppressWarnings(min(esr1_vals, na.rm = TRUE)),
      esr1_max = suppressWarnings(max(esr1_vals, na.rm = TRUE)),
      esr1_frac_zero = frac_zero,
      esr1_frac_na = frac_na,
      median_nCount_RNA = median_nCount_RNA,
      median_nFeature_RNA = median_nFeature_RNA,
      median_pct_mt = median_pct_mt,
      error_msg = NA_character_
    )
  }, error = function(e) {
    fail_row(paste0("ERROR: ", conditionMessage(e)))
  })
}

log_msg("STEP 01: ESR1 QC")
qc_tbl <- map_dfr(seq_len(nrow(samples)), function(i) {
  row <- samples[i, ]
  compute_esr1_qc(
    rds_path = row$rds_path,
    sample_id = row$Sample,
    treatment = row$Treatment,
    loader = row$loader,
    file_ext = row$file_ext
  )
})
write_csv(qc_tbl, FILES$qc_summary)
log_msg("Wrote QC summary:", FILES$qc_summary)

qc_good <- qc_tbl %>%
  filter(
    loaded_ok,
    esr1_frac_na == 0,
    n_cells >= CFG$qc_min_cells,
    esr1_frac_zero < CFG$qc_max_esr1_frac_zero,
    esr1_sd > CFG$qc_min_esr1_sd
  )

qc_weights <- qc_good %>%
  mutate(
    log_n_cells = log10(pmax(n_cells, 10)),
    size_score = rescale(log_n_cells, to = c(0, 1)),
    nonzero_score = 1 - esr1_frac_zero,
    hetero_score = rescale(esr1_iqr, to = c(0, 1)),
    esr1_mean_clipped = pmin(esr1_mean, CFG$qc_esr1_mean_clip),
    mean_score = rescale(esr1_mean_clipped, to = c(0, 1)),
    quality_score = (size_score + nonzero_score + hetero_score + mean_score) / 4,
    sample_weight = quality_score / sum(quality_score)
  )
write_csv(qc_weights, FILES$qc_weights)
log_msg("Wrote QC weights:", FILES$qc_weights)

## ---------------------------------------------------------------
## STEP 02: PER-SAMPLE METACELL ESR1 ~ KINASE CORRELATIONS
## ---------------------------------------------------------------
process_one_sample_corr <- function(row, kinase_signatures_all, cfg) {
  sample_id <- row$Sample
  treatment <- row$Treatment
  path <- row$rds_path
  weight <- row$sample_weight
  
  log_msg("[CORR] Sample:", sample_id, "| Treatment:", treatment)
  
  if (!file.exists(path)) {
    return(list(
      corr = tibble(Sample = sample_id, Treatment = treatment, signature = character(), signature_full = character(),
                    r = numeric(), p_value = numeric(), n_metacells = integer(), n_cells = integer(),
                    sample_weight = weight, error_msg = "file_not_found"),
      sig_coverage = tibble()
    ))
  }
  
  so <- tryCatch(load_seurat_any(path), error = function(e) NULL)
  if (is.null(so)) {
    return(list(
      corr = tibble(Sample = sample_id, Treatment = treatment, signature = character(), signature_full = character(),
                    r = numeric(), p_value = numeric(), n_metacells = integer(), n_cells = integer(),
                    sample_weight = weight, error_msg = "load_failed"),
      sig_coverage = tibble()
    ))
  }
  
  assay_used <- choose_expression_assay(so, cfg$assay_preference)
  DefaultAssay(so) <- assay_used
  so <- ensure_umap_reduction(so, reduction_name = cfg$reduction)
  
  counts_mat <- safe_layer(so, assay = assay_used, layer_name = "counts")
  if (is.null(counts_mat)) {
    return(list(
      corr = tibble(Sample = sample_id, Treatment = treatment, signature = character(), signature_full = character(),
                    r = numeric(), p_value = numeric(), n_metacells = integer(), n_cells = ncol(so),
                    sample_weight = weight, error_msg = "counts_missing"),
      sig_coverage = tibble()
    ))
  }
  genes_universe <- rownames(counts_mat)
  
  sig_filtered <- lapply(kinase_signatures_all, function(gs) intersect(gs, genes_universe))
  sig_lengths <- lengths(sig_filtered)
  keep_idx <- which(sig_lengths >= cfg$min_genes_per_signature)
  sig_filtered <- sig_filtered[keep_idx]
  sig_lengths <- sig_lengths[keep_idx]
  
  sig_cov_tbl <- tibble(
    Sample = sample_id,
    Treatment = treatment,
    signature_full = names(sig_filtered),
    n_genes_in_sig = as.integer(sig_lengths)
  )
  
  if (length(sig_filtered) == 0) {
    return(list(
      corr = tibble(Sample = sample_id, Treatment = treatment, signature = character(), signature_full = character(),
                    r = numeric(), p_value = numeric(), n_metacells = integer(), n_cells = ncol(so),
                    sample_weight = weight, error_msg = "no_signatures_with_genes"),
      sig_coverage = sig_cov_tbl
    ))
  }
  
  meta_cols <- colnames(so@meta.data)
  existing_kin_cols <- grep(cfg$kinase_suffix_pattern, meta_cols, value = TRUE)
  
  if (length(existing_kin_cols) == 0) {
    bp <- BiocParallel::MulticoreParam(workers = cfg$n_workers_ucell)
    so <- AddModuleScore_UCell(
      obj = so,
      features = sig_filtered,
      assay = assay_used,
      slot = "counts",
      BPPARAM = bp,
      name = NULL
    )
    BiocParallel::bpstop(bp)
    meta_cols <- colnames(so@meta.data)
    existing_kin_cols <- grep(cfg$kinase_suffix_pattern, meta_cols, value = TRUE)
  }
  
  if (length(existing_kin_cols) == 0) {
    return(list(
      corr = tibble(Sample = sample_id, Treatment = treatment, signature = character(), signature_full = character(),
                    r = numeric(), p_value = numeric(), n_metacells = integer(), n_cells = ncol(so),
                    sample_weight = weight, error_msg = "no_kinase_columns"),
      sig_coverage = sig_cov_tbl
    ))
  }
  
  sig_names_full <- names(sig_filtered)
  kin_map <- setNames(existing_kin_cols, normalize_sig_name(existing_kin_cols))
  matched_idx <- match(normalize_sig_name(sig_names_full), names(kin_map))
  sig_cols <- unname(kin_map[matched_idx])
  names(sig_cols) <- sig_names_full
  sig_cols <- sig_cols[!is.na(sig_cols)]
  
  if (length(sig_cols) == 0) {
    return(list(
      corr = tibble(Sample = sample_id, Treatment = treatment, signature = character(), signature_full = character(),
                    r = numeric(), p_value = numeric(), n_metacells = integer(), n_cells = ncol(so),
                    sample_weight = weight, error_msg = "no_matching_signature_columns"),
      sig_coverage = sig_cov_tbl
    ))
  }
  
  so_mc <- make_umap_metacells_scRNA(
    seurat_object = so,
    reduction = cfg$reduction,
    cells_per_metacell = cfg$cells_per_metacell
  )
  
  Idents(so_mc) <- "metacell_id"
  agg_expr <- AggregateExpression(
    object = so_mc,
    assays = assay_used,
    group.by = "metacell_id",
    layer = "data",
    verbose = FALSE
  )
  metacell_expression <- agg_expr[[assay_used]]
  
  if (!("ESR1" %in% rownames(metacell_expression))) {
    return(list(
      corr = tibble(Sample = sample_id, Treatment = treatment, signature = character(), signature_full = character(),
                    r = numeric(), p_value = numeric(), n_metacells = dplyr::n_distinct(so_mc$metacell_id), n_cells = ncol(so),
                    sample_weight = weight, error_msg = "ESR1_missing"),
      sig_coverage = sig_cov_tbl
    ))
  }
  
  esr1_expr <- as.numeric(metacell_expression["ESR1", ])
  names(esr1_expr) <- colnames(metacell_expression)
  
  meta_mc <- so_mc@meta.data
  sig_cols <- intersect(sig_cols, colnames(meta_mc))
  sig_split <- split(meta_mc[, sig_cols, drop = FALSE], meta_mc$metacell_id)
  metacell_signatures <- t(vapply(
    sig_split,
    function(df) colMeans(df, na.rm = TRUE),
    numeric(length(sig_cols))
  ))
  rownames(metacell_signatures) <- names(sig_split)
  sig_df <- as.data.frame(metacell_signatures)
  sig_df <- sig_df[names(esr1_expr), , drop = FALSE]
  
  r_vec <- numeric(ncol(sig_df))
  p_vec <- numeric(ncol(sig_df))
  for (i in seq_len(ncol(sig_df))) {
    vec <- sig_df[[i]]
    ct <- suppressWarnings(try(cor.test(esr1_expr, vec, method = "pearson"), silent = TRUE))
    if (inherits(ct, "try-error")) {
      r_vec[i] <- NA_real_
      p_vec[i] <- NA_real_
    } else {
      r_vec[i] <- unname(ct$estimate)
      p_vec[i] <- ct$p.value
    }
  }
  
  sig_full_names <- colnames(sig_df)
  sig_short <- sig_full_names %>%
    str_replace("_UCell$", "") %>%
    str_replace(" human kinase ARCHS4 coexpression$", "") %>%
    str_replace(" human kinase ARCHS4 coexpression_UCell$", "")
  
  corr_tbl <- tibble(
    Sample = sample_id,
    Treatment = treatment,
    signature = sig_short,
    signature_full = sig_full_names,
    r = r_vec,
    p_value = p_vec,
    n_metacells = dplyr::n_distinct(so_mc$metacell_id),
    n_cells = ncol(so),
    sample_weight = weight,
    error_msg = NA_character_
  )
  
  write_csv(corr_tbl, file.path(CFG$out_dir, "per_sample_corr", paste0(sample_id, "_ESR1_kinase_corr.csv")))
  write_csv(
    tibble(
      Sample = sample_id,
      Treatment = treatment,
      assay_used = assay_used,
      n_cells = ncol(so),
      n_features = nrow(so),
      n_metacells = dplyr::n_distinct(so_mc$metacell_id),
      sample_weight = weight,
      ESR1_mean_metacell = mean(esr1_expr),
      ESR1_sd_metacell = sd(esr1_expr)
    ),
    file.path(CFG$out_dir, "per_sample_meta", paste0(sample_id, "_metacell_meta_summary.csv"))
  )
  
  list(corr = corr_tbl, sig_coverage = sig_cov_tbl)
}

log_msg("STEP 02: per-sample metacell ESR1~kinase correlations")
all_corr_list <- list()
all_sigcov_list <- list()
for (i in seq_len(nrow(qc_weights))) {
  row_i <- qc_weights[i, ]
  res_i <- tryCatch(
    process_one_sample_corr(
      row = row_i,
      kinase_signatures_all = kinase_signatures_all,
      cfg = list(
        assay_preference = CFG$assay_preference,
        reduction = CFG$reduction_for_metacells,
        cells_per_metacell = CFG$cells_per_metacell,
        min_genes_per_signature = CFG$min_genes_per_signature,
        n_workers_ucell = CFG$n_workers_ucell,
        kinase_suffix_pattern = CFG$kinase_suffix_pattern
      )
    ),
    error = function(e) {
      log_msg("  Unexpected error in correlation step for", row_i$Sample, ":", conditionMessage(e))
      list(corr = tibble(), sig_coverage = tibble())
    }
  )
  all_corr_list[[i]] <- res_i$corr
  all_sigcov_list[[i]] <- res_i$sig_coverage
}
all_corr <- bind_rows(all_corr_list)
all_sigcov <- bind_rows(all_sigcov_list)
write_csv(all_corr, FILES$corr_all)
write_csv(all_sigcov, FILES$sig_cov)
log_msg("Wrote all-sample correlation table:", FILES$corr_all)
log_msg("Wrote signature coverage table:", FILES$sig_cov)

## ---------------------------------------------------------------
## STEP 03: WEIGHTED META-ANALYSIS OF KINASE CORRELATIONS
## ---------------------------------------------------------------
log_msg("STEP 03: weighted kinase meta-analysis")
qc_sel <- qc_weights %>% select(Sample, sample_weight)
corr_all2 <- all_corr %>%
  left_join(qc_sel, by = "Sample", suffix = c("", ".qc")) %>%
  mutate(sample_weight = coalesce(sample_weight, sample_weight.qc)) %>%
  select(-sample_weight.qc) %>%
  filter(!is.na(r), is.finite(r), r > -0.999999, r < 0.999999) %>%
  mutate(
    sample_weight = if_else(is.na(sample_weight), 1.0, sample_weight),
    n_cells = if_else(is.na(n_cells), 0L, n_cells),
    w_eff = n_cells * (sample_weight ^ CFG$meta_hardness)
  ) %>%
  filter(w_eff > 0)

meta_tbl <- corr_all2 %>%
  group_by(signature, signature_full) %>%
  group_modify(~ meta_one_signature(.x, CFG$meta_min_samples_per_signature)) %>%
  ungroup() %>%
  arrange(desc(r_meta)) %>%
  mutate(rank_desc = row_number())

write_csv(meta_tbl, FILES$corr_meta)
log_msg("Wrote weighted kinase meta-analysis:", FILES$corr_meta)

## ---------------------------------------------------------------
## STEP 04: PER-SAMPLE DRIVER ANALYSIS ACROSS TOP ESR1-NEGATIVE KINASES
## ---------------------------------------------------------------
get_signature_genes <- function(signature_name, signatures_all) {
  i <- which(tolower(names(signatures_all)) == tolower(signature_name))[1]
  if (is.na(i)) return(character(0))
  signatures_all[[i]]
}

process_one_sample_drivers <- function(row, top_sig_names, signatures_all, cfg) {
  sample_id <- row$Sample
  treatment <- row$Treatment
  path <- row$rds_path
  out_file <- file.path(CFG$out_dir, "per_sample_drivers", paste0(sample_id, "_ESR1_kinase_drivers_raw.csv"))
  
  log_msg("[DRIVERS] Sample:", sample_id)
  
  if (!file.exists(path)) return(invisible(NULL))
  so <- tryCatch(load_seurat_any(path), error = function(e) NULL)
  if (is.null(so)) return(invisible(NULL))
  
  assay_used <- choose_expression_assay(so, CFG$assay_preference)
  DefaultAssay(so) <- assay_used
  so <- ensure_umap_reduction(so, reduction_name = cfg$reduction)
  
  so_mc <- make_umap_metacells_scRNA(
    seurat_object = so,
    reduction = cfg$reduction,
    cells_per_metacell = cfg$cells_per_metacell
  )
  
  Idents(so_mc) <- "metacell_id"
  agg_expr <- AggregateExpression(
    object = so_mc,
    assays = assay_used,
    group.by = "metacell_id",
    layer = "data",
    verbose = FALSE
  )
  expression_data <- agg_expr[[assay_used]]
  if (!("ESR1" %in% rownames(expression_data))) return(invisible(NULL))
  
  esr1_vec <- as.numeric(expression_data["ESR1", ])
  names(esr1_vec) <- colnames(expression_data)
  genes_universe <- rownames(expression_data)
  
  all_sig_res <- list()
  for (sig_name in top_sig_names) {
    genes_all <- get_signature_genes(sig_name, signatures_all)
    genes <- intersect(genes_all, genes_universe)
    if (length(genes) < cfg$min_genes_per_sig) next
    
    X <- t(expression_data[genes, , drop = FALSE])
    common_mc <- intersect(names(esr1_vec), rownames(X))
    if (length(common_mc) < 3L) next
    y <- esr1_vec[common_mc]
    X <- X[common_mc, , drop = FALSE]
    
    sig_mean <- rowMeans(X, na.rm = TRUE)
    r_sig <- safe_cor(y, sig_mean, method = "pearson")
    
    per_gene <- map_dfr(colnames(X), function(g) {
      tibble(
        gene = g,
        r_gene = safe_cor(y, X[, g], method = "pearson"),
        p_gene = safe_cortest_p(y, X[, g], method = "pearson")
      )
    })
    
    r_loo <- map_dbl(colnames(X), function(g) {
      if (ncol(X) == 1L) return(NA_real_)
      sig_minus_g <- rowMeans(X[, setdiff(colnames(X), g), drop = FALSE], na.rm = TRUE)
      safe_cor(y, sig_minus_g, method = "pearson")
    })
    
    sig_tbl <- per_gene %>%
      mutate(
        Sample = sample_id,
        Treatment = treatment,
        signature_full = sig_name,
        r_sig = r_sig,
        delta_r = r_sig - r_loo,
        n_genes_sig = length(genes),
        n_metacells = length(common_mc)
      )
    
    all_sig_res[[length(all_sig_res) + 1L]] <- sig_tbl
  }
  
  if (length(all_sig_res) == 0L) return(invisible(NULL))
  drivers_sample <- bind_rows(all_sig_res)
  write_csv(drivers_sample, out_file)
  invisible(NULL)
}

log_msg("STEP 04: per-sample kinase driver analysis")
top_sig_names <- meta_tbl %>% arrange(r_meta) %>% slice_head(n = CFG$K_SIGN) %>% pull(signature_full)
for (i in seq_len(nrow(qc_weights))) {
  row_i <- qc_weights[i, ]
  tryCatch(
    process_one_sample_drivers(
      row = row_i,
      top_sig_names = top_sig_names,
      signatures_all = kinase_signatures_all,
      cfg = list(
        reduction = CFG$reduction_for_metacells,
        cells_per_metacell = CFG$cells_per_metacell,
        min_genes_per_sig = CFG$min_genes_per_signature
      )
    ),
    error = function(e) log_msg("  Unexpected driver error for", row_i$Sample, ":", conditionMessage(e))
  )
}

## ---------------------------------------------------------------
## STEP 05: META-DRIVER SUMMARY ACROSS SAMPLES
## ---------------------------------------------------------------
log_msg("STEP 05: meta-driver summary")
driver_files <- list.files(file.path(CFG$out_dir, "per_sample_drivers"), pattern = "\\.csv$", full.names = TRUE)
if (length(driver_files) > 0L) {
  drivers_all <- map_dfr(driver_files, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    if (!("Sample" %in% colnames(df))) {
      df$Sample <- sub("_ESR1.*$", "", basename(f))
    }
    if (!("signature_full" %in% colnames(df)) && ("signature" %in% colnames(df))) {
      df <- df %>% rename(signature_full = signature)
    }
    df
  })
  
  corr_join <- corr_all2 %>% select(Sample, signature_full, r_sample = r, w_eff)
  
  drivers_join <- drivers_all %>%
    inner_join(corr_join, by = c("Sample", "signature_full")) %>%
    filter(is.finite(r_gene), is.finite(delta_r), is.finite(r_sample), is.finite(w_eff)) %>%
    mutate(
      r_gene_neg = pmax(0, -r_gene),
      delta_neg = pmax(0, -delta_r),
      sig_neg = pmax(0, -r_sample),
      contrib = w_eff * sig_neg * (r_gene_neg + CFG$driver_alpha * delta_neg)
    )
  
  meta_drivers <- drivers_join %>%
    group_by(signature_full, gene) %>%
    summarise(
      score = sum(contrib, na.rm = TRUE),
      n_samples_hit = n_distinct(Sample),
      mean_r_gene = mean(r_gene, na.rm = TRUE),
      mean_delta_r = mean(delta_r, na.rm = TRUE),
      mean_r_sample = mean(r_sample, na.rm = TRUE),
      sum_w_eff = sum(w_eff, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(meta_drivers, FILES$meta_drivers)
  log_msg("Wrote meta-driver table:", FILES$meta_drivers)
} else {
  log_msg("No per-sample driver files found; skipping meta-driver summary")
}

log_msg("==== Finished ESR1 heterogeneity / kinase pipeline ====")