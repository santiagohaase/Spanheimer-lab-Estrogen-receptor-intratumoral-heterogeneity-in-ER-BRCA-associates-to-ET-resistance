library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

#==================================================
# 1. Subset to Control cells and define populations
#==================================================
pdx_control <- subset(PDX1_Control_Tam, idents = "PDX1_Control")

pdx_control$ER_state_pop <- NA_character_
pdx_control$ER_state_pop[pdx_control$seurat_clusters %in% c(6, 0)] <- "low_ESR1_non_prolif"
pdx_control$ER_state_pop[pdx_control$seurat_clusters %in% c(7)]    <- "low_ESR1_prolif"
pdx_control$ER_state_pop[pdx_control$seurat_clusters %in% c(16, 2, 11, 9, 13)] <- "high_ESR1"

pdx_pop_sr_Control <- subset(pdx_control, subset = !is.na(ER_state_pop))

pdx_pop_sr_Control$ER_state_pop <- factor(
  pdx_pop_sr_Control$ER_state_pop,
  levels = c("high_ESR1", "low_ESR1_non_prolif", "low_ESR1_prolif")
)

#==================================================
# 2. Define cell-cycle genes to exclude
#==================================================
cc_genes <- unique(c(
  cc.genes.updated.2019$s.genes,
  cc.genes.updated.2019$g2m.genes
))

#==================================================
# 3. Extract RNA expression matrix
#    RNA data is preferred for Shannon entropy
#==================================================
expr <- GetAssayData(pdx_pop_sr_Control, assay = "RNA", layer = "data")

genes_use <- setdiff(rownames(expr), cc_genes)
expr <- expr[genes_use, , drop = FALSE]

#==================================================
# 4. Function to calculate Shannon entropy per cell
#==================================================
calc_shannon_entropy <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x[x < 0] <- 0
  
  total_expr <- sum(x)
  if (total_expr <= 0) return(NA_real_)
  
  p <- x / total_expr
  p <- p[p > 0]
  
  -sum(p * log2(p))
}

#==================================================
# 5. Calculate per-cell entropy
#==================================================
shannon_entropy <- apply(expr, 2, calc_shannon_entropy)

pdx_pop_sr_Control$Shannon_entropy <- shannon_entropy[colnames(pdx_pop_sr_Control)]

#==================================================
# 6. Residualize entropy for technical covariates
#==================================================
meta_df <- pdx_pop_sr_Control@meta.data

entropy_model <- lm(
  Shannon_entropy ~ log10(nCount_RNA + 1) + log10(nFeature_RNA + 1),
  data = meta_df
)

pdx_pop_sr_Control$Shannon_entropy_resid <- resid(entropy_model)
pdx_pop_sr_Control$Shannon_entropy_fitted <- fitted(entropy_model)

summary(entropy_model)

#==================================================
# 7. Plot raw entropy
#==================================================
p_raw <- ggplot(
  pdx_pop_sr_Control@meta.data,
  aes(x = ER_state_pop, y = Shannon_entropy, fill = ER_state_pop)
) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 13) +
  labs(
    title = "Shannon entropy",
    subtitle = "RNA data, cell-cycle genes removed",
    x = NULL,
    y = "Shannon entropy"
  ) +
  theme(legend.position = "none")

p_raw

#==================================================
# 8. Plot residualized entropy
#==================================================
p_resid <- ggplot(
  pdx_pop_sr_Control@meta.data,
  aes(x = ER_state_pop, y = Shannon_entropy_resid, fill = ER_state_pop)
) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 13) +
  labs(
    title = "Residualized Shannon entropy",
    subtitle = "Adjusted for nCount_RNA and nFeature_RNA",
    x = NULL,
    y = "Residual Shannon entropy"
  ) +
  theme(legend.position = "none")

p_resid

#==================================================
# 9. Statistical tests
#==================================================
pairwise.wilcox.test(
  x = pdx_pop_sr_Control$Shannon_entropy_resid,
  g = pdx_pop_sr_Control$ER_state_pop,
  p.adjust.method = "BH"
)

kruskal.test(
  Shannon_entropy_resid ~ ER_state_pop,
  data = pdx_pop_sr_Control@meta.data
)