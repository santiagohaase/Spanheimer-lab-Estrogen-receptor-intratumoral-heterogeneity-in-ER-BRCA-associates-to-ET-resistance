


#### ESR1-associated kinase correlation and driver-gene analysis in PDX1 control cells

library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(UCell)
library(pheatmap)

###############################################################################
# 1. Inputs
###############################################################################

# Seurat object containing PDX1 control cells
PDX1_Control <- readRDS("/path/to/PDX1_Control.rds")

# Named list of signatures.
# Expected format:
# signatures[["PIM3 human kinase ARCHS4 coexpression"]] = c("GENE1","GENE2",...)
# signatures[["ESR1 human tf ARCHS4 coexpression"]]     = c(...)
if (!exists("signatures")) {
  stop("The object 'signatures' must be loaded before running this script.")
}

# Output directory
out_dir <- "/path/to/output_directory"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. Helper functions
###############################################################################

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

safe_cortest_p <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor.test(x, y, method = method)$p.value)
}

make_umap_metacells_scRNA <- function(seurat_object,
                                      reduction = "umap",
                                      cells_per_metacell = 50,
                                      seed = 1,
                                      verbose = TRUE) {
  set.seed(seed)
  
  emb <- Embeddings(seurat_object, reduction = reduction)
  if (is.null(emb) || nrow(emb) == 0) {
    stop("UMAP reduction not found or empty.")
  }
  
  remaining <- rownames(emb)
  metacell_id <- rep(NA_character_, length(remaining))
  names(metacell_id) <- remaining
  
  mc_counter <- 1
  
  while (length(remaining) > 0) {
    seed_cell <- sample(remaining, 1)
    coords_seed <- emb[seed_cell, , drop = FALSE]
    coords_remaining <- emb[remaining, , drop = FALSE]
    
    d <- sqrt(rowSums((coords_remaining - matrix(
      rep(as.numeric(coords_seed), nrow(coords_remaining)),
      nrow = nrow(coords_remaining),
      byrow = TRUE
    ))^2))
    
    take_n <- min(cells_per_metacell, length(remaining))
    nearest_cells <- names(sort(d))[seq_len(take_n)]
    
    metacell_id[nearest_cells] <- paste0("MC_", mc_counter)
    remaining <- setdiff(remaining, nearest_cells)
    
    if (verbose && mc_counter %% 50 == 0) {
      message("Assigned ", mc_counter, " metacells")
    }
    
    mc_counter <- mc_counter + 1
  }
  
  seurat_object$metacell_id <- metacell_id[colnames(seurat_object)]
  return(seurat_object)
}

aggregate_metacell_matrix <- function(seurat_object, features, assay = "SCT", layer = "data") {
  expr <- GetAssayData(seurat_object, assay = assay, layer = layer)
  features <- intersect(features, rownames(expr))
  if (length(features) == 0) {
    stop("None of the requested features were found in the expression matrix.")
  }
  
  mc_ids <- seurat_object$metacell_id
  if (is.null(mc_ids)) {
    stop("metacell_id metadata is missing.")
  }
  
  valid_cells <- !is.na(mc_ids)
  expr <- expr[features, valid_cells, drop = FALSE]
  mc_ids <- mc_ids[valid_cells]
  
  split_cells <- split(seq_along(mc_ids), mc_ids)
  
  metacell_expression <- sapply(split_cells, function(idx) {
    rowMeans(as.matrix(expr[, idx, drop = FALSE]), na.rm = TRUE)
  })
  
  if (is.vector(metacell_expression)) {
    metacell_expression <- matrix(
      metacell_expression,
      nrow = length(features),
      dimnames = list(features, names(split_cells))
    )
  }
  
  return(metacell_expression)
}

aggregate_metacell_signatures <- function(seurat_object, signature_cols) {
  meta <- seurat_object@meta.data
  
  if (!"metacell_id" %in% colnames(meta)) {
    stop("metacell_id metadata is missing.")
  }
  
  signature_cols <- intersect(signature_cols, colnames(meta))
  if (length(signature_cols) == 0) {
    stop("No requested signature columns found in metadata.")
  }
  
  sig_split <- split(meta[, signature_cols, drop = FALSE], meta$metacell_id)
  metacell_signatures <- t(sapply(sig_split, function(df) colMeans(df, na.rm = TRUE)))
  
  return(as.data.frame(metacell_signatures))
}

###############################################################################
# 3. Optional: UMAP visualization in single cells
###############################################################################

# These plots assume the metadata columns already exist in PDX1_Control
plot_vars <- c(
  "lowESR1_nonprolif_difference",
  "ESR1",
  "lowESR1_prolif_difference",
  "umap_1",
  "umap_2"
)

plot_vars <- intersect(plot_vars, c(colnames(PDX1_Control@meta.data), rownames(PDX1_Control), "umap_1", "umap_2"))

if (all(c("umap_1", "umap_2") %in% plot_vars) && "lowESR1_nonprolif_difference" %in% plot_vars) {
  umap_data <- FetchData(PDX1_Control, vars = c("lowESR1_nonprolif_difference", "umap_1", "umap_2"))
  p1 <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = lowESR1_nonprolif_difference)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_viridis_c(option = "C", limits = c(-0.1, 0.2), na.value = "lightgrey", oob = scales::squish) +
    labs(color = "lowESR1_nonprolif_difference") +
    theme_minimal()
  ggsave(file.path(out_dir, "UMAP_lowESR1_nonprolif_difference_PDX1_Control.pdf"), p1, width = 6, height = 5)
}

if (all(c("umap_1", "umap_2") %in% plot_vars)) {
  umap_data <- FetchData(PDX1_Control, vars = c("ESR1", "umap_1", "umap_2"))
  p2 <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = ESR1)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_viridis_c(option = "C", limits = c(0, 3.5), na.value = "lightgrey", oob = scales::squish) +
    labs(color = "ESR1 Expression") +
    theme_minimal()
  ggsave(file.path(out_dir, "UMAP_ESR1_PDX1_Control.pdf"), p2, width = 6, height = 5)
}

if (all(c("umap_1", "umap_2") %in% plot_vars) && "lowESR1_prolif_difference" %in% plot_vars) {
  umap_data <- FetchData(PDX1_Control, vars = c("lowESR1_prolif_difference", "umap_1", "umap_2"))
  p3 <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = lowESR1_prolif_difference)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_viridis_c(option = "C", limits = c(-0.1, 0.2), na.value = "lightgrey", oob = scales::squish) +
    labs(color = "lowESR1_prolif_difference") +
    theme_minimal()
  ggsave(file.path(out_dir, "UMAP_lowESR1_prolif_difference_PDX1_Control.pdf"), p3, width = 6, height = 5)
}

###############################################################################
# 4. ESR1 single-cell distribution diagnostics
###############################################################################

if ("ESR1" %in% rownames(GetAssayData(PDX1_Control, assay = "SCT", layer = "data"))) {
  esr1_vals <- GetAssayData(PDX1_Control, assay = "SCT", layer = "data")["ESR1", ]
  
  n_cells  <- length(esr1_vals)
  zero_cnt <- sum(esr1_vals == 0, na.rm = TRUE)
  na_cnt   <- sum(is.na(esr1_vals))
  pct_zero <- zero_cnt / n_cells * 100
  pct_na   <- na_cnt / n_cells * 100
  nonzero  <- esr1_vals[is.finite(esr1_vals) & esr1_vals != 0]
  
  esr1_summary <- data.frame(
    metric = c("n_cells", "zero_count", "zero_pct", "na_count", "na_pct", "nonzero_min", "nonzero_max"),
    value = c(
      n_cells,
      zero_cnt,
      pct_zero,
      na_cnt,
      pct_na,
      ifelse(length(nonzero) > 0, min(nonzero), NA),
      ifelse(length(nonzero) > 0, max(nonzero), NA)
    )
  )
  write.csv(esr1_summary, file.path(out_dir, "ESR1_distribution_summary_PDX1_Control.csv"), row.names = FALSE)
  
  pdf(file.path(out_dir, "ESR1_distribution_histogram_PDX1_Control.pdf"), width = 6, height = 5)
  hist(
    esr1_vals,
    breaks = 50,
    main = "ESR1 SCT expression in PDX1 control cells",
    xlab = "SCT data layer value",
    ylab = "Cell count"
  )
  dev.off()
}

###############################################################################
# 5. Generate metacells
###############################################################################

PDX1_Control_mc <- make_umap_metacells_scRNA(
  seurat_object = PDX1_Control,
  reduction = "umap",
  cells_per_metacell = 50,
  verbose = TRUE
)

write.csv(
  data.frame(cell = colnames(PDX1_Control_mc), metacell_id = PDX1_Control_mc$metacell_id),
  file.path(out_dir, "PDX1_Control_metacell_assignments.csv"),
  row.names = FALSE
)

###############################################################################
# 6. Aggregate metacell expression and signature scores
###############################################################################

# Identify kinase-signature columns in metadata
kinase_sig_cols <- grep(
  pattern = "human kinase ARCHS4 coexpression$",
  x = colnames(PDX1_Control_mc@meta.data),
  value = TRUE
)

# Identify ESR1 activity score columns if present
esr1_sig_candidates <- c(
  "ESR1 human tf ARCHS4 coexpression",
  "ESR1_human_tf_ARCHS4_coexpression",
  "ESR1_Encode"
)
esr1_sig_col <- intersect(esr1_sig_candidates, colnames(PDX1_Control_mc@meta.data))

if (length(esr1_sig_col) == 0) {
  warning("No ESR1 activity signature found in metadata. ESR1 TF activity correlation block will be skipped.")
}
if (length(esr1_sig_col) > 1) {
  esr1_sig_col <- esr1_sig_col[1]
}

signature_cols_to_aggregate <- c(kinase_sig_cols, esr1_sig_col)

metacell_signatures <- aggregate_metacell_signatures(
  seurat_object = PDX1_Control_mc,
  signature_cols = signature_cols_to_aggregate
)

metacell_expression <- aggregate_metacell_matrix(
  seurat_object = PDX1_Control_mc,
  features = "ESR1",
  assay = "SCT",
  layer = "data"
)

esr1_expr <- as.numeric(metacell_expression["ESR1", ])
names(esr1_expr) <- colnames(metacell_expression)

# Align rows
metacell_signatures <- metacell_signatures[names(esr1_expr), , drop = FALSE]

###############################################################################
# 7. Correlate metacell ESR1 expression with kinase signatures
###############################################################################

sig_df <- as.data.frame(metacell_signatures[, kinase_sig_cols, drop = FALSE])

res_ESR1_expression_kinases_corr_PDX1_Control <- data.frame(
  kinase = colnames(sig_df),
  r = NA_real_,
  p.value = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(sig_df))) {
  k <- colnames(sig_df)[i]
  vec <- sig_df[[k]]
  res_ESR1_expression_kinases_corr_PDX1_Control$r[i] <- safe_cor(esr1_expr, vec, method = "pearson")
  res_ESR1_expression_kinases_corr_PDX1_Control$p.value[i] <- safe_cortest_p(esr1_expr, vec, method = "pearson")
}

res_ESR1_expression_kinases_corr_PDX1_Control <- res_ESR1_expression_kinases_corr_PDX1_Control %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  arrange(desc(r))

write.csv(
  res_ESR1_expression_kinases_corr_PDX1_Control,
  file.path(out_dir, "ESR1_expression_kinase_correlations_PDX1_Control.csv"),
  row.names = FALSE
)

###############################################################################
# 8. Correlate ESR1 TF activity score with kinase signatures
###############################################################################

if (length(esr1_sig_col) == 1) {
  er_score <- metacell_signatures[[esr1_sig_col]]
  
  res_ESR1_tf_sign_expression_kinases_corr_PDX1_Control <- data.frame(
    kinase = kinase_sig_cols,
    r = NA_real_,
    p.value = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(kinase_sig_cols)) {
    k <- kinase_sig_cols[i]
    vec <- metacell_signatures[[k]]
    res_ESR1_tf_sign_expression_kinases_corr_PDX1_Control$r[i] <- safe_cor(er_score, vec, method = "pearson")
    res_ESR1_tf_sign_expression_kinases_corr_PDX1_Control$p.value[i] <- safe_cortest_p(er_score, vec, method = "pearson")
  }
  
  res_ESR1_tf_sign_expression_kinases_corr_PDX1_Control <- res_ESR1_tf_sign_expression_kinases_corr_PDX1_Control %>%
    mutate(padj = p.adjust(p.value, method = "BH")) %>%
    arrange(desc(r))
  
  write.csv(
    res_ESR1_tf_sign_expression_kinases_corr_PDX1_Control,
    file.path(out_dir, "ESR1_TF_activity_kinase_correlations_PDX1_Control.csv"),
    row.names = FALSE
  )
}

###############################################################################
# 9. Driver genes per negatively correlated kinase signature
###############################################################################

# Use ESR1 expression correlations as the primary ranking
corr_tbl <- res_ESR1_expression_kinases_corr_PDX1_Control %>%
  filter(is.finite(r)) %>%
  mutate(w_sig = pmax(0, -r)) %>%
  arrange(r)

# Keep only negatively correlated signatures with nonzero weight
top_neg_sigs <- corr_tbl %>%
  filter(w_sig > 0)

if (nrow(top_neg_sigs) == 0) {
  stop("No negatively correlated kinase signatures were found.")
}

# Build broader metacell expression matrix for driver analysis
all_kinase_genes <- unique(unlist(signatures[top_neg_sigs$kinase]))
all_kinase_genes <- intersect(all_kinase_genes, rownames(GetAssayData(PDX1_Control_mc, assay = "SCT", layer = "data")))

expression_data <- aggregate_metacell_matrix(
  seurat_object = PDX1_Control_mc,
  features = unique(c("ESR1", all_kinase_genes)),
  assay = "SCT",
  layer = "data"
)

esr1_vec <- as.numeric(expression_data["ESR1", ])
names(esr1_vec) <- colnames(expression_data)

compute_signature_drivers <- function(sig_name_key, w_sig) {
  genes_all <- signatures[[sig_name_key]]
  if (is.null(genes_all)) return(NULL)
  
  genes <- intersect(genes_all, rownames(expression_data))
  genes <- setdiff(genes, "ESR1")
  if (!length(genes)) return(NULL)
  
  X <- t(expression_data[genes, , drop = FALSE])   # metacells x genes
  common_mc2 <- intersect(names(esr1_vec), rownames(X))
  if (length(common_mc2) < 3) return(NULL)
  
  y2 <- esr1_vec[common_mc2]
  X <- X[common_mc2, , drop = FALSE]
  
  if (ncol(X) == 0) return(NULL)
  
  per_gene <- purrr::map_dfr(colnames(X), function(g) {
    tibble(
      gene = g,
      r_gene = safe_cor(y2, X[, g], method = "pearson"),
      p_gene = safe_cortest_p(y2, X[, g], method = "pearson")
    )
  })
  
  sig_mean <- rowMeans(X, na.rm = TRUE)
  r_sig <- safe_cor(y2, sig_mean, method = "pearson")
  
  r_loo <- purrr::map_dbl(colnames(X), function(g) {
    if (ncol(X) == 1) return(NA_real_)
    s_loo <- rowMeans(X[, setdiff(colnames(X), g), drop = FALSE], na.rm = TRUE)
    safe_cor(y2, s_loo, method = "pearson")
  })
  
  per_gene %>%
    mutate(
      signature = sig_name_key,
      r_sig = r_sig,
      delta_r = r_sig - r_loo,
      w_sig = w_sig,
      n_genes_sig = length(genes)
    )
}

drivers_list <- purrr::map2(
  .x = top_neg_sigs$kinase,
  .y = top_neg_sigs$w_sig,
  .f = compute_signature_drivers
)

drivers_all <- bind_rows(drivers_list)

if (nrow(drivers_all) == 0) {
  stop("No driver genes were computed.")
}

###############################################################################
# 10. Score driver genes across kinase signatures
###############################################################################

alpha <- 1

drivers_all_scored <- drivers_all %>%
  mutate(
    r_gene_neg = pmax(0, -r_gene),
    delta_neg  = pmax(0, -delta_r),
    contrib_corr_only    = w_sig * r_gene_neg,
    contrib_corr_plusLOO = w_sig * (r_gene_neg + alpha * delta_neg)
  )

gene_rank <- drivers_all_scored %>%
  group_by(gene) %>%
  summarise(
    score_corr_only    = sum(contrib_corr_only, na.rm = TRUE),
    score_corr_plusLOO = sum(contrib_corr_plusLOO, na.rm = TRUE),
    n_signatures_hit   = n_distinct(signature),
    mean_sig_weight    = mean(w_sig, na.rm = TRUE),
    max_sig_weight     = max(w_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(score_corr_plusLOO))

write.csv(
  drivers_all_scored,
  file.path(out_dir, "ESR1_negative_kinase_per_signature_all_genes_scored_PDX1_Control.csv"),
  row.names = FALSE
)

write.csv(
  gene_rank,
  file.path(out_dir, "ESR1_negative_kinase_DRIVERS_weighted_gene_rank_PDX1_Control.csv"),
  row.names = FALSE
)

###############################################################################
# 11. Optional heatmap of top meta-kinase driver genes
###############################################################################

top_driver_genes <- head(gene_rank$gene, 50)
top_driver_genes <- intersect(top_driver_genes, rownames(expression_data))

if (length(top_driver_genes) > 1) {
  heat_mat <- expression_data[top_driver_genes, , drop = FALSE]
  heat_mat <- t(scale(t(heat_mat)))
  
  pdf(file.path(out_dir, "Top_meta_kinase_driver_genes_heatmap_PDX1_Control.pdf"), width = 8, height = 10)
  pheatmap(
    heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = "Top meta-kinase driver genes across PDX1 metacells"
  )
  dev.off()
}

###############################################################################
# 12. Save key R objects
###############################################################################

saveRDS(PDX1_Control_mc, file.path(out_dir, "PDX1_Control_with_metacells.rds"))
saveRDS(metacell_signatures, file.path(out_dir, "PDX1_metacell_signatures.rds"))
saveRDS(expression_data, file.path(out_dir, "PDX1_metacell_expression_driver_input.rds"))
saveRDS(drivers_all_scored, file.path(out_dir, "PDX1_driver_gene_details.rds"))
saveRDS(gene_rank, file.path(out_dir, "PDX1_meta_kinase_driver_rank.rds"))









# ---------------------- HEATMAP (Top 50 genes × selected signatures) ----------------------

top_genes <- gene_rank %>%
  dplyr::arrange(dplyr::desc(score_corr_plusLOO)) %>%
  dplyr::slice_head(n = TOP_GENES) %>%
  dplyr::pull(gene)

sel_sigs <- intersect(top_neg_sigs$signature, unique(drivers_all_scored$signature))

M_long <- drivers_all_scored %>%
  dplyr::filter(gene %in% top_genes & signature %in% sel_sigs) %>%
  dplyr::select(gene, signature, contrib_corr_plusLOO)

M <- M_long %>%
  tidyr::complete(gene, signature, fill = list(contrib_corr_plusLOO = 0)) %>%
  tidyr::pivot_wider(names_from = signature, values_from = contrib_corr_plusLOO) %>%
  as.data.frame()

rownames(M) <- M$gene
M$gene <- NULL
M <- as.matrix(M)

# keep original signature names for matching annotations
orig_colnames <- colnames(M)

# Row-wise z-score for visualization
M_z <- t(scale(t(M)))
M_z[is.na(M_z)] <- 0

# ---------------------- COLUMN ANNOTATIONS ----------------------
# Build annotation table starting from the actual matrix columns
col_ann_tbl <- tibble::tibble(signature = orig_colnames) %>%
  dplyr::left_join(
    top_neg_sigs %>%
      dplyr::select(signature, signature_label, w_sig),
    by = "signature"
  ) %>%
  dplyr::left_join(
    kinase_cor_tbl %>%
      dplyr::select(signature = signature_key, r_pearson),
    by = "signature"
  )

# check what is missing
col_ann_tbl %>% dplyr::filter(is.na(w_sig) | is.na(r_pearson))

stopifnot(identical(col_ann_tbl$signature, orig_colnames))

# ---------------------- COLORS ----------------------
w_sig_col_fun <- circlize::colorRamp2(
  c(min(col_ann_tbl$w_sig, na.rm = TRUE),
    median(col_ann_tbl$w_sig, na.rm = TRUE),
    max(col_ann_tbl$w_sig, na.rm = TRUE)),
  c("#f0f0f0", "#9ecae1", "#08519c")
)

r_col_fun <- circlize::colorRamp2(
  c(min(col_ann_tbl$r_pearson, na.rm = TRUE), 0, max(col_ann_tbl$r_pearson, na.rm = TRUE)),
  c("#67000d", "#f7f7f7", "#084081")
)

# ---------------------- ANNOTATIONS ----------------------
ha_col <- ComplexHeatmap::HeatmapAnnotation(
  weight    = ComplexHeatmap::anno_simple(col_ann_tbl$w_sig, col = w_sig_col_fun, border = TRUE, gp = grid::gpar(lwd = 0.3)),
  r_pearson = ComplexHeatmap::anno_simple(col_ann_tbl$r_pearson, col = r_col_fun, border = TRUE, gp = grid::gpar(lwd = 0.3)),
  which = "column"
)

ann_legends <- list()
if (sd(col_ann_tbl$w_sig, na.rm = TRUE) > 0) {
  ann_legends <- c(ann_legends, list(
    ComplexHeatmap::Legend(title = "Signature weight (|r|)", col_fun = w_sig_col_fun)
  ))
}
if (sd(col_ann_tbl$r_pearson, na.rm = TRUE) > 0) {
  ann_legends <- c(ann_legends, list(
    ComplexHeatmap::Legend(title = "ESR1 vs signature (r)", col_fun = r_col_fun)
  ))
}

# ---------------------- DISPLAY LABELS ----------------------
shorten_sig <- function(x) {
  x <- sub("(?i)\\s*(human\\s*)?kinase\\s*archs4\\s*coexpression\\s*$", "", x, perl = TRUE)
  x <- trimws(x)
  x <- toupper(x)
  x
}

short_labels <- make.unique(shorten_sig(orig_colnames))
colnames(M_z) <- short_labels

# ---------------------- TITLE / HEATMAP COLORS ----------------------
hm_title <- sprintf(
  "Top %d ESR1-negative kinase drivers\n(values = row-Z of per-signature weighted contributions)",
  nrow(M_z)
)

heat_col_fun <- circlize::colorRamp2(
  c(-2.5, 0, 2.5),
  c("#313695", "#f7f7f7", "#a50026")
)

# ---------------------- HEATMAP ----------------------
ht <- ComplexHeatmap::Heatmap(
  M_z,
  name = "Z(contrib)",
  col = heat_col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = hm_title,
  top_annotation = ha_col,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 45,
  column_dend_height = grid::unit(25, "mm"),
  row_dend_width = grid::unit(20, "mm"),
  heatmap_legend_param = list(title = "Row-Z of contribution", at = c(-2, 0, 2))
)

ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = ann_legends
)