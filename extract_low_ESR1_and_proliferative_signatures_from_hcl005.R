#######################################
# PDX1 Signatures
#######################################

library(Seurat)
library(UCell)
library(dplyr)
library(tibble)
library(ggplot2)
library(scales)

#######################################
# Config
#######################################

pdx1_rds <- "/path/to/PDX_analysis_R/PDX1_Control_Tam.RDS"

merged_qc_out <- "/path/to/SC_RNA_ERHeterogeneity_QC/merged_qc.rds"

signature_out <- "/path/to/SC_RNA_ERHeterogeneity_QC/lowESR1_signatures.rds"

sample_keep <- "PDX1_Control"

non_prolif_clusters <- c("6", "0")
prolif_clusters     <- c("15", "7")

all_low_esr1_clusters <- c(
  "12", "18", "4", "9", "13", "16", "2", "11", "10", "15", "7"
)

rest_for_prolif <- c(
  "12", "18", "0", "6", "4", "9", "13", "16", "2", "11", "10"
)

marker_params <- list(
  p_val_adj = 1e-5,
  logfc_up = 0.25,
  logfc_down = -0.25,
  pct_diff = 0.2,
  n_up = 200,
  n_down = 100
)

#######################################
# Helper functions
#######################################

extract_signature_genes <- function(
    markers,
    direction = c("up", "down"),
    p_val_adj = 1e-5,
    logfc_up = 0.25,
    logfc_down = -0.25,
    pct_diff = 0.2,
    n = 200
) {
  direction <- match.arg(direction)
  
  marker_df <- markers %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene")
  
  if (direction == "up") {
    marker_df %>%
      filter(
        p_val_adj < !!p_val_adj,
        avg_log2FC > !!logfc_up,
        pct.1 - pct.2 > !!pct_diff
      ) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = n) %>%
      pull(gene)
  } else {
    marker_df %>%
      filter(
        p_val_adj < !!p_val_adj,
        avg_log2FC < !!logfc_down,
        pct.2 - pct.1 > !!pct_diff
      ) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = n) %>%
      pull(gene)
  }
}


get_score_col <- function(seurat_obj, signature_name) {
  # UCell column naming can vary by version/settings.
  # This allows either "signature" or "signature_UCell".
  possible_cols <- c(signature_name, paste0(signature_name, "_UCell"))
  
  score_col <- possible_cols[possible_cols %in% colnames(seurat_obj@meta.data)]
  
  if (length(score_col) == 0) {
    stop(
      "Could not find UCell score column for signature: ",
      signature_name,
      "\nTried: ",
      paste(possible_cols, collapse = ", ")
    )
  }
  
  score_col[1]
}


add_lowESR1_scores <- function(seurat_obj, signatures) {
  seurat_obj <- AddModuleScore_UCell(
    seurat_obj,
    features = signatures,
    name = NULL
  )
  
  up_np   <- get_score_col(seurat_obj, "lowESR1_nonprolif_up")
  down_np <- get_score_col(seurat_obj, "lowESR1_nonprolif_dw")
  up_pr   <- get_score_col(seurat_obj, "lowESR1_prolif_up")
  down_pr <- get_score_col(seurat_obj, "lowESR1_prolif_dw")
  
  seurat_obj$lowESR1_nonprolif_difference <-
    seurat_obj@meta.data[[up_np]] - seurat_obj@meta.data[[down_np]]
  
  seurat_obj$lowESR1_prolif_difference <-
    seurat_obj@meta.data[[up_pr]] - seurat_obj@meta.data[[down_pr]]
  
  seurat_obj
}


plot_umap_score <- function(
    seurat_obj,
    score_col,
    limits = NULL,
    point_size = 1,
    alpha = 0.5
) {
  umap_data <- FetchData(
    seurat_obj,
    vars = c(score_col, "umap_1", "umap_2")
  )
  
  ggplot(
    umap_data,
    aes(
      x = umap_1,
      y = umap_2,
      color = .data[[score_col]]
    )
  ) +
    geom_point(alpha = alpha, size = point_size) +
    scale_color_viridis_c(
      option = "C",
      limits = limits,
      na.value = "lightgrey",
      oob = scales::squish
    ) +
    labs(color = score_col) +
    theme_minimal()
}

#######################################
# Load and prepare PDX1 control object
#######################################

PDX1_Control_Tam <- readRDS(pdx1_rds)

Idents(PDX1_Control_Tam) <- "sampleID"

PDX1_Control <- subset(
  PDX1_Control_Tam,
  subset = sampleID == sample_keep
)

rm(PDX1_Control_Tam)

DefaultAssay(PDX1_Control) <- "SCT"
Idents(PDX1_Control) <- "seurat_clusters"

#######################################
# Find markers
#######################################

da_genes_non_prolif_vs_all_rest <- FindMarkers(
  object = PDX1_Control,
  ident.1 = non_prolif_clusters,
  ident.2 = all_low_esr1_clusters,
  assay = "SCT"
)

da_genes_prolif_vs_all_rest <- FindMarkers(
  object = PDX1_Control,
  ident.1 = prolif_clusters,
  ident.2 = rest_for_prolif,
  assay = "SCT"
)

#######################################
# Extract signatures
#######################################

signatures <- list(
  lowESR1_nonprolif_up = extract_signature_genes(
    da_genes_non_prolif_vs_all_rest,
    direction = "up",
    p_val_adj = marker_params$p_val_adj,
    logfc_up = marker_params$logfc_up,
    pct_diff = marker_params$pct_diff,
    n = marker_params$n_up
  ),
  
  lowESR1_nonprolif_dw = extract_signature_genes(
    da_genes_non_prolif_vs_all_rest,
    direction = "down",
    p_val_adj = marker_params$p_val_adj,
    logfc_down = marker_params$logfc_down,
    pct_diff = marker_params$pct_diff,
    n = marker_params$n_down
  ),
  
  lowESR1_prolif_up = extract_signature_genes(
    da_genes_prolif_vs_all_rest,
    direction = "up",
    p_val_adj = marker_params$p_val_adj,
    logfc_up = marker_params$logfc_up,
    pct_diff = marker_params$pct_diff,
    n = marker_params$n_up
  ),
  
  lowESR1_prolif_dw = extract_signature_genes(
    da_genes_prolif_vs_all_rest,
    direction = "down",
    p_val_adj = marker_params$p_val_adj,
    logfc_down = marker_params$logfc_down,
    pct_diff = marker_params$pct_diff,
    n = marker_params$n_down
  )
)

saveRDS(signatures, signature_out)

#######################################
# Score PDX1 control object
#######################################

PDX1_Control <- add_lowESR1_scores(
  seurat_obj = PDX1_Control,
  signatures = signatures
)

#######################################
# Example plots on PDX1 control
#######################################

plot_umap_score(
  PDX1_Control,
  score_col = get_score_col(PDX1_Control, "lowESR1_prolif_up"),
  limits = c(-0.1, 0.4)
)

plot_umap_score(
  PDX1_Control,
  score_col = "lowESR1_nonprolif_difference",
  limits = c(-0.2, 0.4)
)

plot_umap_score(
  PDX1_Control,
  score_col = get_score_col(PDX1_Control, "lowESR1_nonprolif_up"),
  limits = c(0, 0.6)
)

#######################################
# Apply signatures to merged_qc
#######################################

merged_qc <- add_lowESR1_scores(
  seurat_obj = merged_qc,
  signatures = signatures
)

saveRDS(merged_qc, merged_qc_out)