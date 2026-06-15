#!/usr/bin/env Rscript

## ============================================================
## Build single-cell breast cancer atlas from public scRNA samples
##
## For each sample:
##   1. Read 10X matrix/barcode/features triplet
##   2. Create Seurat object
##   3. Perform QC
##   4. Normalize with SCTransform
##   5. Run PCA/UMAP
##   6. Annotate cell types with SingleR / HPCA
##   7. Run inferCNV using immune cells as reference
##   8. Compute per-cell CNV burden
##   9. Extract epithelial high-CNV cancer cells
##  10. Save full object, inferCNV object, cancer-cell object, and summary
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(celldex)
  library(SingleR)
  library(GenomicFeatures)
  library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(infercnv)
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(matrixStats)
})

## ---------------------------
## Global helpers
## ---------------------------

now_stamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

safe_dir_create <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

log_line <- function(..., log_file = NULL) {
  msg <- paste0(now_stamp(), " | ", paste(..., collapse = " "))
  message(msg)
  if (!is.null(log_file)) {
    cat(msg, "\n", file = log_file, append = TRUE)
  }
}

## Seurat v4/v5-compatible assay accessor
get_assay_mat <- function(obj, assay, layer_or_slot = c("counts", "data", "scale.data")) {
  layer_or_slot <- match.arg(layer_or_slot)
  
  out <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer_or_slot),
    error = function(e1) {
      tryCatch(
        Seurat::GetAssayData(obj, assay = assay, slot = layer_or_slot),
        error = function(e2) NULL
      )
    }
  )
  
  if (is.null(out)) {
    stop("Could not retrieve ", layer_or_slot, " from assay ", assay)
  }
  
  out
}

## ---------------------------
## Find matrix/barcode/features triplets
## ---------------------------

find_filtered_triplets <- function(root_dir, sample_ids = NULL) {
  
  features_global <- file.path(root_dir, "GSE161529_features.tsv.gz")
  
  if (!file.exists(features_global)) {
    stop(
      "Expected global features file not found: ", features_global,
      "\nExpected a flat folder with GSE161529_features.tsv.gz and sample-specific matrix/barcode files."
    )
  }
  
  mat_files <- list.files(
    root_dir,
    pattern = "matrix\\.mtx\\.gz$",
    full.names = TRUE
  )
  
  if (!length(mat_files)) {
    stop("No *-matrix.mtx.gz files found in: ", root_dir)
  }
  
  if (!is.null(sample_ids)) {
    mat_files <- mat_files[
      vapply(
        mat_files,
        function(f) {
          any(vapply(sample_ids, function(id) grepl(id, basename(f), fixed = TRUE), logical(1)))
        },
        logical(1)
      )
    ]
  }
  
  if (!length(mat_files)) {
    stop("No matrix files remained after sample_ids filtering.")
  }
  
  rows <- lapply(mat_files, function(mf) {
    prefix <- sub("-matrix\\.mtx\\.gz$", "", basename(mf))
    barcode_file <- file.path(root_dir, paste0(prefix, "-barcodes.tsv.gz"))
    
    if (!file.exists(barcode_file)) {
      warning("Missing barcode file for ", prefix)
      return(NULL)
    }
    
    data.frame(
      sample_id = prefix,
      matrix    = normalizePath(mf),
      barcodes  = normalizePath(barcode_file),
      features  = normalizePath(features_global),
      stringsAsFactors = FALSE
    )
  })
  
  out <- bind_rows(rows)
  
  if (!nrow(out)) {
    stop("No complete matrix/barcode/features triplets found.")
  }
  
  out
}

## ---------------------------
## Raw counts helper
## ---------------------------

get_counts_matrix <- function(obj) {
  assay_names <- names(obj@assays)
  
  if ("RNA" %in% assay_names) {
    m <- tryCatch(get_assay_mat(obj, assay = "RNA", layer_or_slot = "counts"),
                  error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  }
  
  if ("SCT" %in% assay_names) {
    m <- tryCatch(get_assay_mat(obj, assay = "SCT", layer_or_slot = "counts"),
                  error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  }
  
  stop("No usable raw counts found in RNA or SCT assays.")
}

## ---------------------------
## Gene order for inferCNV
## ---------------------------

make_gene_order <- function(counts) {
  symbols <- rownames(counts)
  
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  keep <- which(!is.na(entrez))
  if (length(keep) < 2000) {
    warning("Few genes mapped to ENTREZ IDs; check whether rownames are gene symbols.")
  }
  
  symbols <- symbols[keep]
  entrez  <- entrez[keep]
  
  tx_genes <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  GenomeInfoDb::seqlevelsStyle(tx_genes) <- "UCSC"
  
  tx_sub <- tx_genes[as.character(S4Vectors::mcols(tx_genes)$gene_id) %in% entrez]
  
  if (length(tx_sub) == 0) {
    stop("No overlap between gene symbols/ENTREZ IDs and TxDb hg38 genes.")
  }
  
  sym_sub <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = as.character(S4Vectors::mcols(tx_sub)$gene_id),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  gene_df <- data.frame(
    gene  = unname(sym_sub),
    chr   = as.character(seqnames(tx_sub)),
    start = start(tx_sub),
    end   = end(tx_sub),
    stringsAsFactors = FALSE
  )
  
  standard_chr <- paste0("chr", c(1:22, "X"))
  
  gene_df %>%
    filter(
      !is.na(gene),
      gene %in% rownames(counts),
      chr %in% standard_chr
    ) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(factor(chr, levels = standard_chr), start)
}

## ---------------------------
## Cell-type label groups
## ---------------------------

IMMUNE_LABELS_DEFAULT <- c(
  "B_cell", "T_cells", "NK_cell", "Monocyte", "Macrophage", "DC", "CMP", "GMP",
  "Macrophages", "Dendritic_cell", "Naive_B_cell", "Memory_B_cell",
  "CD4+_T_cell", "CD8+_T_cell", "NK_cells", "Monocytes"
)

EPITHELIAL_LABELS_DEFAULT <- c(
  "Epithelial_cells", "Epithelial_cell", "Epithelial"
)

## ---------------------------
## inferCNV burden score
## ---------------------------

compute_infercnv_score <- function(infer_obj, min_variable_genes = 50, variance_floor = 5e-4) {
  
  mat <- infer_obj@expr.data
  
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    stop("inferCNV expr.data matrix is empty.")
  }
  
  vars <- matrixStats::rowVars(as.matrix(mat))
  vars_nonzero <- vars[vars > 0]
  
  if (length(vars_nonzero) == 0) {
    mat_sub <- mat
  } else {
    q80 <- as.numeric(stats::quantile(vars_nonzero, 0.80, na.rm = TRUE))
    cutoff <- max(q80, variance_floor)
    var_idx <- which(vars >= cutoff)
    
    if (length(var_idx) < min_variable_genes) {
      mat_sub <- mat
    } else {
      mat_sub <- mat[var_idx, , drop = FALSE]
    }
  }
  
  cnv_score <- Matrix::colMeans(abs(mat_sub - 1), na.rm = TRUE)
  
  if (is.null(names(cnv_score))) {
    names(cnv_score) <- colnames(mat_sub)
  }
  
  cnv_score
}

## ---------------------------
## Process one sample
## ---------------------------

process_one_sample <- function(
    sample_row,
    out_dir,
    ref_mat,
    ref_lbl,
    min_genes = 300,
    max_genes = 8000,
    max_mt = 20,
    cap_ref = 5000,
    threads = 16,
    run_umap = TRUE,
    immune_labels = IMMUNE_LABELS_DEFAULT,
    epithelial_labels = EPITHELIAL_LABELS_DEFAULT,
    epithelial_cnv_quantile = 0.60
) {
  
  sample_id <- sample_row$sample_id
  sample_out <- file.path(out_dir, sample_id)
  safe_dir_create(sample_out)
  
  log_file <- file.path(sample_out, paste0("run_", sample_id, ".log"))
  
  log_line("========== Sample:", sample_id, "==========", log_file = log_file)
  
  summary_row <- tryCatch({
    
    ## ---------------------------
    ## 1. Read sample matrix
    ## ---------------------------
    
    log_line("Reading matrix:", sample_row$matrix, log_file = log_file)
    
    counts <- Seurat::ReadMtx(
      mtx      = sample_row$matrix,
      features = sample_row$features,
      cells    = sample_row$barcodes
    )
    
    seu <- Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_id,
      assay = "RNA",
      min.cells = 3,
      min.features = 0
    )
    
    raw_cells <- ncol(seu)
    raw_genes <- nrow(seu)
    
    log_line("Loaded", raw_cells, "cells x", raw_genes, "genes", log_file = log_file)
    
    ## ---------------------------
    ## 2. QC
    ## ---------------------------
    
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
    
    keep <- (
      seu$nFeature_RNA >= min_genes &
        seu$nFeature_RNA <= max_genes &
        seu$percent.mt <= max_mt
    )
    
    kept_cells <- sum(keep)
    seu <- subset(seu, cells = colnames(seu)[keep])
    
    log_line(
      "QC kept", kept_cells, "/", raw_cells,
      "min_genes =", min_genes,
      "max_genes =", max_genes,
      "max_mt =", max_mt,
      log_file = log_file
    )
    
    if (ncol(seu) < 50) {
      stop("Too few cells after QC: ", ncol(seu))
    }
    
    ## ---------------------------
    ## 3. SCTransform / PCA / UMAP
    ## ---------------------------
    
    seu <- Seurat::SCTransform(
      seu,
      assay = "RNA",
      variable.features.n = 3000,
      conserve.memory = TRUE,
      return.only.var.genes = FALSE,
      verbose = FALSE
    )
    
    seu <- Seurat::RunPCA(seu, verbose = FALSE)
    
    if (isTRUE(run_umap)) {
      seu <- Seurat::RunUMAP(seu, dims = 1:30, verbose = FALSE)
    }
    
    log_line("SCTransform/PCA/UMAP complete", log_file = log_file)
    
    ## ---------------------------
    ## 4. SingleR annotation
    ## ---------------------------
    
    test_mat <- get_assay_mat(seu, assay = "SCT", layer_or_slot = "data")
    common <- intersect(rownames(test_mat), rownames(ref_mat))
    
    if (length(common) < 1000) {
      log_line("WARNING: only", length(common), "genes overlap with HPCA reference", log_file = log_file)
    }
    
    sr <- SingleR::SingleR(
      test   = as.matrix(test_mat[common, , drop = FALSE]),
      ref    = as.matrix(ref_mat[common, , drop = FALSE]),
      labels = ref_lbl
    )
    
    seu$SingleR_labels     <- sr$pruned.labels
    seu$SingleR_labels_raw <- sr$labels
    
    seu$SingleR_labels[is.na(seu$SingleR_labels)] <- "Unknown"
    
    write.csv(
      as.data.frame(table(seu$SingleR_labels)),
      file.path(sample_out, "singler_label_counts.csv"),
      row.names = FALSE
    )
    
    log_line("SingleR complete", log_file = log_file)
    
    ## ---------------------------
    ## 5. inferCNV using immune reference cells
    ## ---------------------------
    
    counts_raw <- get_counts_matrix(seu)
    meta <- seu@meta.data
    meta$cell <- rownames(meta)
    
    immune_present <- intersect(immune_labels, unique(meta$SingleR_labels))
    
    if (!length(immune_present)) {
      stop("No immune reference cells detected by SingleR.")
    }
    
    meta$infercnv_group <- as.character(meta$SingleR_labels)
    meta$infercnv_group[meta$infercnv_group %in% immune_present] <- "immune_cells"
    meta$infercnv_group[is.na(meta$infercnv_group)] <- "unknown"
    
    keep_mask <- rep(TRUE, nrow(meta))
    
    if (is.finite(cap_ref)) {
      imm_idx <- which(meta$infercnv_group == "immune_cells")
      
      if (length(imm_idx) > cap_ref) {
        set.seed(1)
        keep_idx <- sample(imm_idx, cap_ref)
        drop_idx <- setdiff(imm_idx, keep_idx)
        keep_mask[drop_idx] <- FALSE
      }
    }
    
    cells_use <- meta$cell[keep_mask]
    keep_cols <- intersect(colnames(counts_raw), cells_use)
    
    counts_ic <- counts_raw[, keep_cols, drop = FALSE]
    meta_ic <- meta[match(colnames(counts_ic), meta$cell), , drop = FALSE]
    
    n_immune <- sum(meta_ic$infercnv_group == "immune_cells", na.rm = TRUE)
    
    if (n_immune < 2) {
      stop("Fewer than 2 immune reference cells after filtering.")
    }
    
    gene_df <- make_gene_order(counts_ic)
    
    counts_ic <- counts_ic[intersect(rownames(counts_ic), gene_df$gene), , drop = FALSE]
    counts_ic <- counts_ic[match(gene_df$gene, rownames(counts_ic)), , drop = FALSE]
    
    out_infer <- file.path(sample_out, paste0("infercnv_out_rna_immune_ref_", sample_id))
    safe_dir_create(out_infer)
    
    anno_file <- file.path(out_infer, "annotations.tsv")
    gene_order_file <- file.path(out_infer, "gene_order_hg38.tsv")
    counts_file <- file.path(out_infer, "counts.tsv")
    
    write.table(
      data.frame(cell = colnames(counts_ic), group = meta_ic$infercnv_group),
      file = anno_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    write.table(
      gene_df,
      file = gene_order_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    write.table(
      as.matrix(counts_ic),
      file = counts_file,
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
    
    log_line("Running inferCNV with", n_immune, "immune reference cells", log_file = log_file)
    
    set.seed(1)
    
    infer_obj <- infercnv::CreateInfercnvObject(
      raw_counts_matrix = counts_file,
      annotations_file  = anno_file,
      delim = "\t",
      gene_order_file   = gene_order_file,
      ref_group_names   = "immune_cells"
    )
    
    infer_obj <- infercnv::run(
      infercnv_obj = infer_obj,
      cutoff = 0.1,
      denoise = TRUE,
      cluster_by_groups = TRUE,
      out_dir = out_infer,
      num_threads = threads,
      HMM = FALSE,
      analysis_mode = "subclusters"
    )
    
    infer_rds <- file.path(out_infer, "infercnv_object.rds")
    saveRDS(infer_obj, infer_rds)
    
    log_line("inferCNV complete", log_file = log_file)
    
    ## ---------------------------
    ## 6. Add CNV burden score
    ## ---------------------------
    
    cnv_score <- compute_infercnv_score(infer_obj)
    
    seu$infercnv_score <- NA_real_
    
    common_cells <- intersect(colnames(seu), names(cnv_score))
    seu$infercnv_score[common_cells] <- cnv_score[common_cells]
    
    log_line(
      "Mapped inferCNV score to",
      length(common_cells), "/", ncol(seu),
      "cells",
      log_file = log_file
    )
    
    ## ---------------------------
    ## 7. Extract epithelial high-CNV cancer cells
    ## ---------------------------
    
    epithelial_present <- intersect(epithelial_labels, unique(seu$SingleR_labels))
    
    if (!length(epithelial_present)) {
      stop("No epithelial cells detected by SingleR.")
    }
    
    epithelial_cells <- colnames(seu)[seu$SingleR_labels %in% epithelial_present]
    
    if (length(epithelial_cells) < 10) {
      stop("Too few epithelial cells detected: ", length(epithelial_cells))
    }
    
    epi_scores <- seu$infercnv_score[epithelial_cells]
    epi_scores_valid <- epi_scores[!is.na(epi_scores)]
    
    if (length(epi_scores_valid) < 10) {
      stop("Too few epithelial cells with valid inferCNV scores.")
    }
    
    cnv_thr <- as.numeric(stats::quantile(
      epi_scores_valid,
      probs = epithelial_cnv_quantile,
      na.rm = TRUE
    ))
    
    seu$cnv_group <- NA_character_
    seu$cnv_group[epithelial_cells] <- ifelse(
      seu$infercnv_score[epithelial_cells] >= cnv_thr,
      "Epithelial_CNV_high",
      "Epithelial_CNV_low"
    )
    
    cancer_cells <- colnames(seu)[seu$cnv_group == "Epithelial_CNV_high"]
    
    if (length(cancer_cells) < 10) {
      stop("Too few epithelial high-CNV cancer cells: ", length(cancer_cells))
    }
    
    seu_cancer <- subset(seu, cells = cancer_cells)
    
    ## ---------------------------
    ## 8. Save outputs
    ## ---------------------------
    
    full_rds <- file.path(sample_out, paste0(sample_id, "_seurat_sct_singler_infercnv.rds"))
    cancer_rds <- file.path(sample_out, paste0("seu_", sample_id, "_Epithelial_cells_cancer.RDS"))
    
    saveRDS(seu, full_rds)
    saveRDS(seu_cancer, cancer_rds)
    
    log_line("Saved full object:", full_rds, log_file = log_file)
    log_line("Saved cancer-cell object:", cancer_rds, log_file = log_file)
    
    data.frame(
      sample_id = sample_id,
      status = "success",
      error_msg = NA_character_,
      cells_raw = raw_cells,
      cells_kept_qc = kept_cells,
      n_genes_raw = raw_genes,
      pct_mt_median = round(median(seu$percent.mt, na.rm = TRUE), 2),
      n_genes_median = round(median(seu$nFeature_RNA, na.rm = TRUE)),
      n_counts_median = round(median(seu$nCount_RNA, na.rm = TRUE)),
      singler_unique_labels = length(unique(seu$SingleR_labels)),
      immune_ref_n = n_immune,
      epithelial_n = length(epithelial_cells),
      cnv_threshold = cnv_thr,
      cancer_cells_n = length(cancer_cells),
      cancer_fraction_of_epithelial = length(cancer_cells) / length(epithelial_cells),
      full_seurat_rds = normalizePath(full_rds),
      cancer_seurat_rds = normalizePath(cancer_rds),
      infercnv_rds = normalizePath(infer_rds),
      infercnv_dir = normalizePath(out_infer),
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    
    log_line("FAILED:", conditionMessage(e), log_file = log_file)
    
    data.frame(
      sample_id = sample_id,
      status = "failed",
      error_msg = conditionMessage(e),
      cells_raw = NA_integer_,
      cells_kept_qc = NA_integer_,
      n_genes_raw = NA_integer_,
      pct_mt_median = NA_real_,
      n_genes_median = NA_real_,
      n_counts_median = NA_real_,
      singler_unique_labels = NA_integer_,
      immune_ref_n = NA_integer_,
      epithelial_n = NA_integer_,
      cnv_threshold = NA_real_,
      cancer_cells_n = NA_integer_,
      cancer_fraction_of_epithelial = NA_real_,
      full_seurat_rds = NA_character_,
      cancer_seurat_rds = NA_character_,
      infercnv_rds = NA_character_,
      infercnv_dir = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  
  summary_row
}

## ---------------------------
## Main atlas function
## ---------------------------

run_atlas_pipeline <- function(
    root_dir,
    out_dir,
    sample_ids = NULL,
    min_genes = 300,
    max_genes = 8000,
    max_mt = 20,
    cap_ref = 5000,
    threads = 16,
    run_umap = TRUE,
    epithelial_cnv_quantile = 0.60
) {
  
  safe_dir_create(out_dir)
  
  main_log <- file.path(out_dir, "atlas_pipeline.log")
  log_line("Starting atlas pipeline", log_file = main_log)
  log_line("Input root:", root_dir, log_file = main_log)
  log_line("Output dir:", out_dir, log_file = main_log)
  
  log_line("Loading HPCA reference", log_file = main_log)
  ref_se <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
  ref_mat <- SummarizedExperiment::assay(ref_se, "logcounts")
  ref_lbl <- ref_se$label.main
  
  samp_df <- find_filtered_triplets(root_dir, sample_ids)
  
  log_line("Found", nrow(samp_df), "samples", log_file = main_log)
  
  all_summary <- lapply(seq_len(nrow(samp_df)), function(i) {
    process_one_sample(
      sample_row = samp_df[i, ],
      out_dir = out_dir,
      ref_mat = ref_mat,
      ref_lbl = ref_lbl,
      min_genes = min_genes,
      max_genes = max_genes,
      max_mt = max_mt,
      cap_ref = cap_ref,
      threads = threads,
      run_umap = run_umap,
      epithelial_cnv_quantile = epithelial_cnv_quantile
    )
  })
  
  summary_df <- bind_rows(all_summary)
  
  summary_file <- file.path(out_dir, "summary_per_sample.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  log_line("Wrote summary:", summary_file, log_file = main_log)
  log_line("Finished atlas pipeline", log_file = main_log)
  
  invisible(summary_df)
}

## ---------------------------
## Run
## ---------------------------

root_dir <- "/path/to/GSE161529_samples"

out_dir <- "/path/to/SC_RNA_GSE161529_atlas_processing"

## NULL = all samples with matrix.mtx.gz files in root_dir
## Or use c("ER-MH0025", "ER-MH0032") to test selected samples
sample_ids <- NULL

summary_df <- run_atlas_pipeline(
  root_dir = root_dir,
  out_dir = out_dir,
  sample_ids = sample_ids,
  min_genes = 300,
  max_genes = 8000,
  max_mt = 20,
  cap_ref = 5000,
  threads = 16,
  run_umap = TRUE,
  epithelial_cnv_quantile = 0.60
)

print(summary_df)