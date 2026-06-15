


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# ============================================================
# PURPOSE
# ============================================================
# Train a cancer-cell-specific low_ESR1 signature using a seed
# signature derived from PDX1 / HCl005 cancer cells.
#
# The seed signature is projected onto a panel of matched
# cancer/TME sample pairs listed in a CSV file. Samples with the
# highest abundance of low_ESR1-like cancer cells are prioritized.
#
# For each selected sample:
#   1) score the cancer-only object with the HCl005 seed signature
#   2) define target low_ESR1 cancer cells using the per-sample frac_top05
#   3) map those target cancer cells into the matched full TME Seurat
#   4) perform DE: target cancer cells vs all remaining TME cells
#   5) collect UP/DOWN genes across samples
#
# Finally, derive a consensus low_ESR1 signature across selected samples.
# ============================================================


# ============================================================
# INPUTS
# ============================================================
csv_file <- "route/to/seurats_for_signature_low_ESR1.csv"
out_root <- "low_ESR1_DE_signature_panel"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# HCl005 / PDX1-derived seed signature
# These should already exist in memory or be loaded here.
# Example:
# low_ESR1_up   <- scan("route/to/HCl005_low_ESR1_up.txt", what = character())
# low_ESR1_down <- scan("route/to/HCl005_low_ESR1_down.txt", what = character())

stopifnot(exists("low_ESR1_up"), exists("low_ESR1_down"))


# ============================================================
# PARAMETERS
# ============================================================
score_col <- "low_ESR1_difference"

# sample selection
use_top_n_samples <- 10
min_frac_top05 <- 0.15

# DE parameters
de_logfc_threshold <- 0.10
de_min_pct         <- 0.05

# signature cutoffs within each sample
padj_cut <- 1e-10
lfc_up   <- 0.50
lfc_down <- -0.50

# cap number of genes retained per sample
cap_up   <- 500
cap_down <- 500

# minimum number of target cells needed
min_target_cells <- 50

# minimum support for consensus
min_support_k <- 3


# ============================================================
# HELPERS
# ============================================================
log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

strip_prefix_to_match <- function(x) {
  sub("^.*?_", "", x)
}

detect_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) >= 1) return(hit[1])
  if (required) stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

clean_path <- function(p) {
  p <- as.character(p)
  p <- str_trim(p)
  p <- str_replace_all(p, '^"+|"+$', "")
  p <- str_replace_all(p, "^'+|'+$", "")
  p
}

read_seurat_any <- function(path) {
  path <- clean_path(path)
  if (!file.exists(path)) stop("File does not exist: ", path)
  
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  if (ext == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop("Package 'qs' is required to read .qs files.")
    }
    return(qs::qread(path))
  }
  stop("Unsupported extension for file: ", path)
}

choose_assay <- function(obj) {
  if ("SCT" %in% names(obj@assays)) return("SCT")
  if ("RNA" %in% names(obj@assays)) return("RNA")
  stop("Neither SCT nor RNA assay found.")
}

map_cancer_to_tme_cells <- function(cancer_cells, tme_cells) {
  overlap <- intersect(cancer_cells, tme_cells)
  if (length(overlap) > 0) return(unique(overlap))
  
  cancer_stripped <- strip_prefix_to_match(cancer_cells)
  tme_stripped    <- strip_prefix_to_match(tme_cells)
  
  tme_map <- setNames(tme_cells, tme_stripped)
  matched <- intersect(cancer_stripped, tme_stripped)
  mapped  <- unname(tme_map[matched])
  mapped  <- mapped[!is.na(mapped)]
  
  unique(mapped)
}

add_low_ESR1_difference <- function(obj, up_genes, down_genes, score_name = score_col) {
  DefaultAssay(obj) <- choose_assay(obj)
  
  obj <- AddModuleScore_UCell(
    obj,
    features = list(
      low_ESR1_up = up_genes,
      low_ESR1_down = down_genes
    ),
    name = NULL
  )
  
  up_col <- grep("^low_ESR1_up", colnames(obj@meta.data), value = TRUE)
  down_col <- grep("^low_ESR1_down", colnames(obj@meta.data), value = TRUE)
  
  stopifnot(length(up_col) == 1, length(down_col) == 1)
  
  obj[[score_name]] <- obj[[up_col]][, 1] - obj[[down_col]][, 1]
  obj
}

safe_write_lines <- function(x, file) {
  x <- unique(x)
  x <- x[nzchar(x)]
  writeLines(x, file)
}


# ============================================================
# READ SAMPLE TABLE
# ============================================================
tbl <- read_csv(csv_file, show_col_types = FALSE)

sample_col <- detect_col(tbl, c("sample_id", "sample", "Sample", "SampleID"), TRUE, "sample column")
cancer_col <- detect_col(tbl, c("cancer_cells_seurat", "cancer_path", "cancer_file", "cancer_seurat", "cancer"), TRUE, "cancer column")
tme_col    <- detect_col(tbl, c("TME_Seurat", "tme_path", "tme_file", "tme_seurat", "tme"), TRUE, "TME column")
frac_col   <- detect_col(tbl, c("frac_top05", "frac_top5", "frac_top5pct", "frac_low_ESR1", "top_frac", "fraction"), TRUE, "frac_top05 column")

tbl <- tbl %>%
  transmute(
    sample_id   = as.character(.data[[sample_col]]),
    cancer_path = clean_path(.data[[cancer_col]]),
    tme_path    = clean_path(.data[[tme_col]]),
    frac_top05  = as.numeric(.data[[frac_col]])
  ) %>%
  filter(!is.na(sample_id), !is.na(cancer_path), !is.na(tme_path), !is.na(frac_top05)) %>%
  arrange(desc(frac_top05))

tbl_selected <- tbl %>%
  filter(frac_top05 >= min_frac_top05)

if (!is.null(use_top_n_samples)) {
  tbl_selected <- tbl_selected %>% slice_head(n = use_top_n_samples)
}

write_csv(tbl,          file.path(out_root, "samples_all_ranked.csv"))
write_csv(tbl_selected, file.path(out_root, "samples_selected_for_training.csv"))

log_msg("Total samples in table:", nrow(tbl))
log_msg("Selected samples for training:", nrow(tbl_selected))

if (nrow(tbl_selected) == 0) {
  stop("No samples passed selection criteria.")
}


# ============================================================
# MAIN TRAINING LOOP
# ============================================================
all_up_lists <- list()
all_down_lists <- list()
run_log <- list()

for (i in seq_len(nrow(tbl_selected))) {
  sample_id   <- tbl_selected$sample_id[i]
  cancer_path <- tbl_selected$cancer_path[i]
  tme_path    <- tbl_selected$tme_path[i]
  frac_target <- tbl_selected$frac_top05[i]
  
  sample_outdir <- file.path(out_root, sample_id)
  dir.create(sample_outdir, recursive = TRUE, showWarnings = FALSE)
  
  msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sample_id, "|", ..., "\n")
  
  msg("Starting")
  msg("frac_top05 =", signif(frac_target, 4))
  msg("cancer =", cancer_path)
  msg("tme =", tme_path)
  
  cancer <- read_seurat_any(cancer_path)
  tme    <- read_seurat_any(tme_path)
  
  cancer_cells <- colnames(cancer)
  tme_cells    <- colnames(tme)
  
  mapped_cancer_in_tme <- map_cancer_to_tme_cells(cancer_cells, tme_cells)
  tme$is_cancer_from_cancerRDS <- colnames(tme) %in% mapped_cancer_in_tme
  
  msg("Cancer cells:", length(cancer_cells))
  msg("TME cells:", length(tme_cells))
  msg("Mapped cancer cells into TME:", length(mapped_cancer_in_tme))
  
  # score both objects with HCl005 seed signature
  cancer <- add_low_ESR1_difference(cancer, low_ESR1_up, low_ESR1_down)
  tme    <- add_low_ESR1_difference(tme, low_ESR1_up, low_ESR1_down)
  
  # define target cells using the per-sample low_ESR1 abundance
  x <- cancer[[score_col]][, 1]
  q <- 1 - frac_target
  cut <- as.numeric(quantile(x, probs = q, na.rm = TRUE))
  cancer$target_low_ESR1_cancer <- x >= cut
  target_cancer_cells <- colnames(cancer)[cancer$target_low_ESR1_cancer %in% TRUE]
  
  msg("Target quantile:", signif(q, 4))
  msg("Target cutoff:", signif(cut, 4))
  msg("Target cancer cells:", length(target_cancer_cells))
  
  if (length(target_cancer_cells) < min_target_cells) {
    msg("Skipping: too few target cancer cells.")
    run_log[[sample_id]] <- data.frame(
      sample_id = sample_id,
      status = "skipped_too_few_target_in_cancer",
      frac_top05 = frac_target,
      n_cancer = length(cancer_cells),
      n_tme = length(tme_cells),
      n_target_cancer = length(target_cancer_cells),
      n_target_in_tme = NA_integer_,
      n_up = NA_integer_,
      n_down = NA_integer_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # map those target cancer cells into matched TME
  target_in_tme <- map_cancer_to_tme_cells(target_cancer_cells, colnames(tme))
  tme$target_low_ESR1_cancer <- colnames(tme) %in% target_in_tme
  
  msg("Target cancer cells mapped into TME:", sum(tme$target_low_ESR1_cancer))
  
  if (sum(tme$target_low_ESR1_cancer) < min_target_cells) {
    msg("Skipping: too few target cancer cells mapped into TME.")
    run_log[[sample_id]] <- data.frame(
      sample_id = sample_id,
      status = "skipped_too_few_target_in_tme",
      frac_top05 = frac_target,
      n_cancer = length(cancer_cells),
      n_tme = length(tme_cells),
      n_target_cancer = length(target_cancer_cells),
      n_target_in_tme = sum(tme$target_low_ESR1_cancer),
      n_up = NA_integer_,
      n_down = NA_integer_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # DE: target cancer cells vs all remaining TME
  DefaultAssay(tme) <- choose_assay(tme)
  Idents(tme) <- factor(tme$target_low_ESR1_cancer, levels = c(FALSE, TRUE))
  
  de <- FindMarkers(
    tme,
    ident.1 = TRUE,
    ident.2 = FALSE,
    logfc.threshold = de_logfc_threshold,
    min.pct = de_min_pct,
    only.pos = FALSE
  )
  
  de$gene <- rownames(de)
  de <- de %>% arrange(p_val_adj, desc(avg_log2FC))
  write.csv(file = file.path(sample_outdir, "DE_targetLow_ESR1Cancer_vs_restOfTME.csv"), de, row.names = FALSE)
  
  up_sig <- de %>%
    filter(!is.na(p_val_adj), p_val_adj < padj_cut, avg_log2FC >= lfc_up) %>%
    arrange(desc(avg_log2FC)) %>%
    pull(gene) %>%
    unique() %>%
    head(cap_up)
  
  down_sig <- de %>%
    filter(!is.na(p_val_adj), p_val_adj < padj_cut, avg_log2FC <= lfc_down) %>%
    arrange(avg_log2FC) %>%
    pull(gene) %>%
    unique() %>%
    head(cap_down)
  
  safe_write_lines(up_sig,   file.path(sample_outdir, "signature_UP_low_ESR1.txt"))
  safe_write_lines(down_sig, file.path(sample_outdir, "signature_DOWN_low_ESR1.txt"))
  
  all_up_lists[[sample_id]] <- up_sig
  all_down_lists[[sample_id]] <- down_sig
  
  run_log[[sample_id]] <- data.frame(
    sample_id = sample_id,
    status = "ok",
    frac_top05 = frac_target,
    n_cancer = length(cancer_cells),
    n_tme = length(tme_cells),
    n_target_cancer = length(target_cancer_cells),
    n_target_in_tme = sum(tme$target_low_ESR1_cancer),
    n_up = length(up_sig),
    n_down = length(down_sig),
    stringsAsFactors = FALSE
  )
  
  saveRDS(tme, file.path(sample_outdir, "tme_scored_targetLow_ESR1Cells.rds"))
  msg("Done")
}

run_log_df <- bind_rows(run_log)
write.csv(run_log_df, file.path(out_root, "run_log.csv"), row.names = FALSE)


# ============================================================
# CONSENSUS SIGNATURE
# ============================================================
ok_samples <- run_log_df %>% filter(status == "ok")

if (nrow(ok_samples) == 0) {
  stop("No successful samples available for consensus.")
}

up_freq <- sort(table(unlist(all_up_lists)), decreasing = TRUE)
down_freq <- sort(table(unlist(all_down_lists)), decreasing = TRUE)

up_freq_df <- data.frame(
  gene = names(up_freq),
  n_samples = as.integer(up_freq),
  stringsAsFactors = FALSE
)

down_freq_df <- data.frame(
  gene = names(down_freq),
  n_samples = as.integer(down_freq),
  stringsAsFactors = FALSE
)

write.csv(up_freq_df,   file.path(out_root, "consensus_UP_low_ESR1_gene_frequency.csv"), row.names = FALSE)
write.csv(down_freq_df, file.path(out_root, "consensus_DOWN_low_ESR1_gene_frequency.csv"), row.names = FALSE)

cons_up <- up_freq_df %>%
  filter(n_samples >= min_support_k) %>%
  pull(gene) %>%
  head(cap_up)

cons_down <- down_freq_df %>%
  filter(n_samples >= min_support_k) %>%
  pull(gene) %>%
  head(cap_down)

safe_write_lines(cons_up,   file.path(out_root, paste0("CONSENSUS_UP_low_ESR1_supportGE", min_support_k, ".txt")))
safe_write_lines(cons_down, file.path(out_root, paste0("CONSENSUS_DOWN_low_ESR1_supportGE", min_support_k, ".txt")))

summary_df <- data.frame(
  metric = c(
    "n_total_samples_in_csv",
    "n_selected_samples",
    "n_successful_samples",
    "consensus_up_genes",
    "consensus_down_genes"
  ),
  value = c(
    nrow(tbl),
    nrow(tbl_selected),
    nrow(ok_samples),
    length(cons_up),
    length(cons_down)
  )
)
write.csv(summary_df, file.path(out_root, "training_summary.csv"), row.names = FALSE)

log_msg("Training complete")
log_msg("Total samples:", nrow(tbl))
log_msg("Selected samples:", nrow(tbl_selected))
log_msg("Successful samples:", nrow(ok_samples))
log_msg("Consensus UP genes:", length(cons_up))
log_msg("Consensus DOWN genes:", length(cons_down))
log_msg("Output directory:", out_root)