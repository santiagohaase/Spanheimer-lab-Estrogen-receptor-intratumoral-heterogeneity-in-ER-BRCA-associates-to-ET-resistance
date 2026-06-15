

###################### Cancer Cell stemness states


library(Seurat)
library(dplyr)
library(UCell)
library(ggplot2)
library(purrr)
library(tibble)

# ---------------------------
# Choose assay and layer
# ---------------------------
score_assay <- "SCT"   # or "RNA"
score_layer <- "data"

# ---------------------------
# Example signature list
# Each signature can have:
#   $positive = c(...)
#   $negative = c(...)
# ---------------------------
signature_list <- list(
  
  Luminal_Progenitor = list(
    positive = c("KIT", "ELF5", "SOX4", "ALDH1A3", "JAG1", "PROM1"),
    negative = c("ESR1", "PGR", "TFF1", "FOXA1", "GATA3")
  ),
  
  BCSC_ALDH = list(
    positive = c("ALDH1A1", "ALDH1A3", "PROM1", "CXCR4", "ITGA6"),
    negative = c("ESR1", "PGR", "KRT18")
  ),
  
  CD44hi_CD24low = list(
    positive = c("CD44", "ITGA6", "CXCL8", "VIM"),
    negative = c("CD24", "EPCAM", "KRT8", "KRT18")
  ),
  
  EMT_Plasticity = list(
    positive = c("VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1"),
    negative = c("EPCAM", "CDH1", "KRT8", "KRT18", "MUC1")
  ),
  
  Quiescence_Persister = list(
    positive = c("CDKN1A", "CDKN1B", "ATF3", "JUN", "FOS", "KLF4", "TXNIP"),
    negative = c("MKI67", "TOP2A", "PCNA", "UBE2C", "CENPF")
  ),
  
  WNT_Stem = list(
    positive = c("TCF7L2", "AXIN2", "LGR5", "MYC", "LEF1"),
    negative = c("ESR1", "PGR")
  )
)


# genes available in chosen assay
available_genes <- rownames(PDX1_Control[[score_assay]])

clean_signature_list <- lapply(signature_list, function(sig) {
  pos <- intersect(sig$positive, available_genes)
  neg <- intersect(sig$negative, available_genes)
  
  list(
    positive = unique(pos),
    negative = unique(neg)
  )
})

# inspect how many genes survived per signature
signature_qc <- tibble(
  signature = names(clean_signature_list),
  n_positive = sapply(clean_signature_list, function(x) length(x$positive)),
  n_negative = sapply(clean_signature_list, function(x) length(x$negative))
)


print(signature_qc)

# Build separate UCell signatures for positive and negative components
ucell_features <- list()

for (sig_name in names(clean_signature_list)) {
  sig <- clean_signature_list[[sig_name]]
  
  if (length(sig$positive) > 0) {
    ucell_features[[paste0(sig_name, "_pos")]] <- sig$positive
  }
  
  if (length(sig$negative) > 0) {
    ucell_features[[paste0(sig_name, "_neg")]] <- sig$negative
  }
}


print(names(ucell_features))




PDX1_Control <- AddModuleScore_UCell(
  obj = PDX1_Control,
  features = ucell_features,
  assay = score_assay,
  name = NULL
)






for (sig_name in names(clean_signature_list)) {
  pos_col <- paste0(sig_name, "_pos")
  neg_col <- paste0(sig_name, "_neg")
  signed_col <- paste0(sig_name, "_UCell_signed")
  
  pos_exists <- pos_col %in% colnames(PDX1_Control@meta.data)
  neg_exists <- neg_col %in% colnames(PDX1_Control@meta.data)
  
  if (pos_exists && neg_exists) {
    PDX1_Control[[signed_col]] <- PDX1_Control@meta.data[[pos_col]] -
      PDX1_Control@meta.data[[neg_col]]
  } else if (pos_exists && !neg_exists) {
    PDX1_Control[[signed_col]] <- PDX1_Control@meta.data[[pos_col]]
  } else if (!pos_exists && neg_exists) {
    PDX1_Control[[signed_col]] <- -PDX1_Control@meta.data[[neg_col]]
  }
}



PDX1_Control$intrinsic_tolerant_group <- ifelse(
  PDX1_Control$seurat_clusters %in% c("0", "6"),
  "clusters_0_6",
  "other_clusters"
)

PDX1_Control$intrinsic_tolerant_group <- factor(
  PDX1_Control$intrinsic_tolerant_group,
  levels = c("other_clusters", "clusters_0_6")
)

table(PDX1_Control$intrinsic_tolerant_group)





score_cols <- c(
  grep("_UCell_signed$", colnames(PDX1_Control@meta.data), value = TRUE),
  grep("_AMS_signed$", colnames(PDX1_Control@meta.data), value = TRUE)
)

score_test_results <- lapply(score_cols, function(sc) {
  df <- PDX1_Control@meta.data |>
    dplyr::select(intrinsic_tolerant_group, all_of(sc)) |>
    dplyr::rename(score = all_of(sc)) |>
    dplyr::filter(!is.na(score))
  
  if (length(unique(df$intrinsic_tolerant_group)) < 2) return(NULL)
  
  wt <- wilcox.test(score ~ intrinsic_tolerant_group, data = df)
  
  summary_df <- df |>
    dplyr::group_by(intrinsic_tolerant_group) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(score, na.rm = TRUE),
      median = median(score, na.rm = TRUE),
      .groups = "drop"
    )
  
  mean_other <- summary_df$mean[summary_df$intrinsic_tolerant_group == "other_clusters"]
  mean_06    <- summary_df$mean[summary_df$intrinsic_tolerant_group == "clusters_0_6"]
  med_other  <- summary_df$median[summary_df$intrinsic_tolerant_group == "other_clusters"]
  med_06     <- summary_df$median[summary_df$intrinsic_tolerant_group == "clusters_0_6"]
  
  tibble(
    score = sc,
    mean_clusters_0_6 = mean_06,
    mean_other = mean_other,
    median_clusters_0_6 = med_06,
    median_other = med_other,
    delta_mean = mean_06 - mean_other,
    delta_median = med_06 - med_other,
    p_value = wt$p.value
  )
}) |>
  bind_rows() |>
  mutate(FDR = p.adjust(p_value, method = "fdr")) |>
  arrange(FDR, desc(delta_mean))

print(score_test_results)






top_scores <- score_test_results$score[1:min(8, nrow(score_test_results))]

for (sc in top_scores) {
  p <- VlnPlot(
    object = PDX1_Control,
    features = sc,
    group.by = "intrinsic_tolerant_group",
    pt.size = 0
  ) +
    ggtitle(sc) +
    theme_bw()
  
  print(p)
}
