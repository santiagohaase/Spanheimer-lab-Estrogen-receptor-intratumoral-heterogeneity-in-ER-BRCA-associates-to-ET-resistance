################## Trajectory PDX3



##### trajectory analysis PDX3_control to PDX3_tam
R.version.string

Sys.which("R")
R.version.string
.libPaths()

save.image("/path/to/all.RData")
library(remotes)
install_github("r-spatial/sf")

Sys.setenv(CXXFLAGS="-std=c++17")
install.packages("units", dependencies = TRUE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages('devtools')
install.packages('usethis')

devtools::install_github('cole-trapnell-lab/monocle3')

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))


BiocManager::install('Rsamtools')
install.packages("Signac")
install.packages("Signac")
Sys.setenv(GDAL_CONFIG = "/path/to/miniconda3/bin/gdal-config")
install.packages("terra", type = "source")

Sys.setenv(PKG_CONFIG_PATH = "/path/to/miniconda3/lib/pkgconfig")


Sys.setenv(GDAL_CONFIG = "/path/to/miniconda3/bin/gdal-config")
Sys.setenv(LD_LIBRARY_PATH = "/path/to/miniconda3/lib")
Sys.setenv(PROJ_LIB = "/path/to/miniconda3/share/proj")
install.packages("terra", type = "source")

install.packages('terra')

library(Seurat)
library(Signac)

library(monocle3)
library(SeuratWrappers) # Contains as.cell_data_set() converter
remotes::install_github('satijalab/seurat-wrappers')


####################### Subset Seurat

PDX_combined_PDX3 <- readRDS('/path/to/PDX_analysis_R/PDX_combined_PDX3.rds')



DimPlot(PDX_combined_PDX3)


levels(Idents(PDX_combined_PDX3))


PDX3_Control_Tam <- subset(x = PDX_combined_PDX3, idents = c("PDX3_Control", "PDX3_Tam"))


###############
library(ggplot2)

DefaultAssay(PDX3_Control_Tam) <- "SCT"

# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(PDX3_Control_Tam, 
                       vars = c("ESR1", "umap_1", "umap_2"))

# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = ESR1)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(0, 3.5),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "ESR1 Expression") +
  theme_minimal()


DimPlot(PDX3_Control_Tam, reduction = 'umap')




library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(UCell)


#################


cds <- as.cell_data_set(PDX3_Control_Tam)
# Suppose you've already set up or transferred your UMAP from Seurat
reducedDims(cds)$UMAP <- Embeddings(PDX3_Control_Tam, "umap")

# Optionally cluster
cds <- cluster_cells(cds)  # uses the existing "UMAP" dims

# Learn graph (no "reduction_method" argument needed)
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# If you do want pseudotime
cds <- order_cells(cds)

plot_cells(cds, 
           color_cells_by = "pseudotime", 
           cell_size = 1,
           show_trajectory_graph = TRUE)

?plot_cells

###############





library(dplyr)
library(tibble)
library(magrittr)


PDX1_Control_Tam <- readRDS(
  "/path/to/PDX_analysis_R/PDX1_Control_Tam.RDS"
)

Idents(PDX1_Control_Tam) <- "sampleID"

PDX1_Control <- subset(PDX1_Control_Tam, subset = sampleID == "PDX1_Control")

rm(PDX1_Control_Tam)

DefaultAssay(PDX1_Control) <- 'SCT'

Idents(PDX1_Control) <- 'seurat_clusters'

da_genes_non_prolif_vs_all_rest <- FindMarkers(
  object  = PDX1_Control,
  ident.1 = c('6', '0'),
  ident.2 = c('12', '18', '4', '9', '13', '16', '2', '11', '10', '15', '7'),
  assay   = 'SCT'
)


da_genes_high_ESR1_all_rest <- FindMarkers(
  object  = PDX1_Control,
  ident.1 = c('16', '2', '9', '13', '11'),
  ident.2 = c('0', '6', '4', '0', '18', '12', '15', '7'),
  assay   = 'SCT'
)


da_genes_prolif_vs_all_rest <- FindMarkers(
  object  = PDX1_Control,
  ident.1 = c('15', '7'),
  ident.2 = c('12', '18', '0', '6', '4', '9', '13', '16', '2', '11', '10'),
  assay   = 'SCT'
)


prolif_up <- da_genes_prolif_vs_all_rest %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 1e-5, avg_log2FC > 0.25, (pct.1 - pct.2) > 0.2) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%   # strongest UP first
  dplyr::slice_head(n = 200) %>%
  dplyr::pull(gene)


prolif_down <- da_genes_prolif_vs_all_rest %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 1e-5, avg_log2FC < -0.25, (pct.2 - pct.1) > 0.2) %>%
  dplyr::arrange(avg_log2FC) %>%        # most negative first (strongest down)
  dplyr::slice_head(n = 100) %>%
  dplyr::pull(gene)


non_prolif_up <- da_genes_non_prolif_vs_all_rest %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 1e-5, avg_log2FC > 0.25, (pct.1 - pct.2) > 0.2) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%   # strongest UP first
  dplyr::slice_head(n = 200) %>%
  dplyr::pull(gene)

non_prolif_down <- da_genes_non_prolif_vs_all_rest %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 1e-5, avg_log2FC < -0.25, (pct.2 - pct.1) > 0.2) %>%
  dplyr::arrange(avg_log2FC) %>%        # most negative first (strongest down)
  dplyr::slice_head(n = 100) %>%
  dplyr::pull(gene)



High_ESR1_down <- da_genes_high_ESR1_all_rest %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 1e-5, avg_log2FC < -0.25, (pct.2 - pct.1) > 0.2) %>%
  dplyr::arrange(avg_log2FC) %>%        # most negative first (strongest down)
  dplyr::slice_head(n = 100) %>%
  dplyr::pull(gene)



High_ESR1_up <- da_genes_high_ESR1_all_rest %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 1e-5, avg_log2FC > 0.25, (pct.1 - pct.2) > 0.2) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%   # strongest UP first
  dplyr::slice_head(n = 200) %>%
  dplyr::pull(gene)






PDX3_Control_Tam <- AddModuleScore_UCell(PDX3_Control_Tam, features = list(
  lowESR1_nonprolif_up = non_prolif_up,
  lowESR1_nonprolif_dw = non_prolif_down,
  lowESR1_prolif_up = prolif_up,
  lowESR1_prolif_dw = prolif_down
), name = NULL)





# column names created by UCell (usually "<name>_UCell")
up_np   <- "lowESR1_nonprolif_up"
down_np <- "lowESR1_nonprolif_dw"
up_pr   <- "lowESR1_prolif_up"
down_pr <- "lowESR1_prolif_dw"



# add difference metadata
PDX3_Control_Tam$lowESR1_nonprolif_difference <- PDX3_Control_Tam[[up_np]] - PDX3_Control_Tam[[down_np]]
PDX3_Control_Tam$lowESR1_prolif_difference    <- PDX3_Control_Tam[[up_pr]] - PDX3_Control_Tam[[down_pr]]



# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(PDX3_Control_Tam, 
                       vars = c("lowESR1_nonprolif_difference", "umap_1", "umap_2"))


# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = lowESR1_nonprolif_difference)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(-0.1,0.4),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "lowESR1_nonprolif_difference") +
  theme_minimal()





# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(PDX3_Control_Tam, 
                       vars = c("lowESR1_nonprolif_up", "umap_1", "umap_2"))


# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = lowESR1_nonprolif_up)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(-0.1,0.4),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "lowESR1_nonprolif_up") +
  theme_minimal()






















library(RColorBrewer)
library(ggplot2)

FeaturePlot(PDX3_Control_Tam, features = "ESR1", reduction = "umap") +
  scale_colour_gradientn(limits = c(0.5,3.2), colours = rev(brewer.pal(n = 11, name = "RdBu")), oob = scales::squish)


##################
## Annotate clusters

# 6) See which neighbor-graphs (SNN, WNN, etc.) have been generated
names(PDX3_Control_Tam@graphs)

Reductions(PDX3_Control_Tam)

options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB


PDX3_Control_Tam <- FindClusters(
  object = PDX3_Control_Tam,
  graph.name = "wknn",    # <-- Important!
  resolution = 2          # Adjust as needed
)



DimPlot(PDX3_Control_Tam, 
        reduction = "umap", 
        label = TRUE, 
        repel = TRUE)




DimPlot(
  PDX3_Control_Tam, 
  reduction = "umap",
  label = TRUE,
  cells.highlight = CellsByIdentities(PDX3_Control_Tam, idents = 14),
  sizes.highlight = 0.05,
  cols.highlight = "red"
)



DimPlot(
  PDX3_Control_and_Vand, 
  reduction = "umap",
  label = TRUE,
  cells.highlight = CellsByIdentities(PDX3_Control_and_Vand, idents = 3),
  sizes.highlight = 0.05,
  cols.highlight = "red"
)



#####################


library(dplyr)
library(Seurat)
library(ggplot2)

# Convert the Seurat metadata to a data frame with "cell_id" from row names:
meta_data <- data.frame(
  cell_id = rownames(PDX3_Control_Tam@meta.data),
  PDX3_Control_Tam@meta.data,
  stringsAsFactors = FALSE
)

# Create a new group column, labeling only cluster 2 or 3 for PDX3_Vand vs PDX3_Control
meta_data <- meta_data %>%
  mutate(compareGroup = case_when(
    sampleID == "PDX3_Tam"    & seurat_clusters %in% c("7", "8", "14") ~ "PDX3_Tam_7_8_14",
    sampleID == "PDX3_Control" & seurat_clusters %in% c("14") ~ "PDX3_Control_14",
    TRUE                                                           ~ NA_character_
  ))

meta_data_sub <- meta_data %>%
  filter(!is.na(compareGroup))

kinase_cols <- grep("kinase", colnames(meta_data_sub), value = TRUE)
kinase_cols

results_list <- list()

for (k in kinase_cols) {
  tam_vals <- meta_data_sub %>% 
    filter(compareGroup == "PDX3_Tam_7_8_14") %>% 
    pull(k)
  
  ctrl_vals <- meta_data_sub %>%
    filter(compareGroup == "PDX3_Control_14") %>%
    pull(k)
  
  # Perform the t-test
  ttest_res <- t.test(tam_vals, ctrl_vals)
  
  # Store the results
  results_list[[k]] <- data.frame(
    signature     = k,
    mean_vand     = mean(tam_vals, na.rm = TRUE),
    mean_control  = mean(ctrl_vals, na.rm = TRUE),
    p_value       = ttest_res$p.value
  )
}

# Combine into a single results data frame
results_df <- do.call(rbind, results_list)
results_df


write.csv(results_df, "results_kinase_signatures_PDX3_Tam_clusterss_8_11_12_vs_PDX3_control_cluster_5.csv")



##################
library(Seurat)
library(dplyr)

# 1. Convert the Seurat metadata to a data frame including cell IDs
meta_data <- data.frame(
  cell_id = rownames(PDX3_Control_Tam@meta.data),
  PDX3_Control_Tam@meta.data,
  stringsAsFactors = FALSE
)

# 2. Create a new group column that labels cells:
#    - "PDX3_Tam_7_8_14" for cells from sampleID "PDX3_Tam" in clusters 8, 11, or 12
#    - "PDX3_Control_14" for cells from sampleID "PDX3_Control" in cluster 5
meta_data <- meta_data %>%
  mutate(compareGroup = case_when(
    sampleID == "PDX3_Tam"    & seurat_clusters %in% c("7", "8", "14") ~ "PDX3_Tam_7_8_14",
    sampleID == "PDX3_Control" & seurat_clusters %in% c("14") ~ "PDX3_Control_14",
    TRUE ~ NA_character_
  ))

# 3. Update the Seurat object's metadata with the new compareGroup column.
#    Ensure that the rownames of meta_data match the cell names in the Seurat object.
rownames(meta_data) <- meta_data$cell_id
PDX3_Control_Tam <- AddMetaData(PDX3_Control_Tam, metadata = meta_data["compareGroup"])

# 4. Subset the Seurat object to include only cells that have been assigned a compareGroup
cells_to_keep <- rownames(PDX3_Control_Tam@meta.data)[!is.na(PDX3_Control_Tam@meta.data$compareGroup)]
filtered_obj <- subset(PDX3_Control_Tam, cells = cells_to_keep)

# 5. Set the cell identities to the compareGroup column so that the DE analysis knows the grouping.
Idents(filtered_obj) <- filtered_obj@meta.data$compareGroup

# Optional: Check the number of cells in each group
table(Idents(filtered_obj))

filtered_obj <- PrepSCTFindMarkers(filtered_obj)

# 6. Run differential expression analysis on the "SCT" assay
de_results <- FindMarkers(
  object = filtered_obj,
  ident.1 = "PDX3_Tam_7_8_14",   # Group with clusters 8, 11, and 12 from PDX3_Tam
  ident.2 = "PDX3_Control_14",      # Group with cluster 5 from PDX3_Control
  assay = "SCT"                  # Specify that the analysis should be done on the SCT assay
)

# 7. View the top differentially expressed genes
head(de_results)

write.csv(de_results, "de_results_PDX3_Tam_7_8_14_vs_PDX3_Control_14.csv")


################################## kinase PCA



# 1. Identify metadata columns with "kinase"
kinase_cols <- grep("kinase", colnames(PDX3_Control_Tam@meta.data), value = TRUE)

# 2. Extract the data (ensure these columns are numeric)
kinase_data <- PDX3_Control_Tam@meta.data[, kinase_cols]

# Optionally, check the data type and convert if needed:
# kinase_data <- sapply(kinase_data, as.numeric)

# 3. Scale the data (centering and scaling is common before PCA)
kinase_data_scaled <- scale(kinase_data)

# 4. Run PCA
pca_results <- prcomp(kinase_data_scaled, center = TRUE, scale. = TRUE)

# View a summary of the PCA
summary(pca_results)

# 5. Visualize the PCA results (for example, a scatter plot of PC1 vs. PC2)
plot(pca_results$x[,1:2], 
     xlab = "PC1", 
     ylab = "PC2", 
     main = "PCA of Kinase Signature Scores",
     pch = 10, col = "blue")




# Create a DimReduc object
kinase_pca <- CreateDimReducObject(embeddings = pca_results$x,
                                   key = "PCK_",
                                   assay = DefaultAssay(PDX3_Control_Tam))

# Add it to the Seurat object
PDX3_Control_Tam[["kinasePCA"]] <- kinase_pca

Reductions(PDX3_Control_Tam)

head(PDX3_Control_Tam)

Idents(PDX3_Control_Tam)<- "sampleID"

DimPlot(PDX3_Control_Tam, reduction = "kinasePCA" )

# Now, you can use Seurat's DimPlot or other functions:
FeaturePlot(PDX3_Control_Tam, features="ESR1" ,reduction = "kinasePCA")

library(RColorBrewer)
library(ggplot2)

# Now, you can use Seurat's DimPlot or other functions:
DimPlot(PDX3_Control_Tam, reduction = "kinasePCA")


FeaturePlot(PDX3_Control_Tam, features = "ESR1", reduction = "kinasePCA") +
  scale_colour_gradientn(limits = c(0.5,3.2), colours = rev(brewer.pal(n = 11, name = "RdBu")), oob = scales::squish)



# Example: boxplot comparison for PC1 scores between groups
pca_scores <- pca_results$x  # Matrix with PCA scores for each cell
boxplot(PC1 ~ PDX3_Control_and_Vand@meta.data$sampleID,
        main = "PC1 Scores by Treatment Group",
        xlab = "Treatment", ylab = "PC1 Score")


boxplot(pca_scores[, 1] ~ PDX3_Control_and_Vand@meta.data$sampleID,
        main = "PC1 Scores by Treatment Group",
        xlab = "Treatment", ylab = "PC1 Score")

###
# Example: boxplot comparison for PC1 scores between groups
pca_scores <- pca_results$x  # Matrix with PCA scores for each cell

boxplot(pca_scores[, 2] ~ PDX3_Control_and_Vand@meta.data$sampleID,
        main = "PC2 Scores by Treatment Group",
        xlab = "Treatment", ylab = "PC2 Score")

boxplot(pca_scores[, 3] ~ PDX3_Control_and_Vand@meta.data$sampleID,
        main = "PC3 Scores by Treatment Group",
        xlab = "Treatment", ylab = "PC3 Score")

boxplot(pca_scores[, 4] ~ PDX3_Control_and_Vand@meta.data$sampleID,
        main = "PC4 Scores by Treatment Group",
        xlab = "Treatment", ylab = "PC4 Score")

boxplot(pca_scores[, 5] ~ PDX3_Control_and_Vand@meta.data$sampleID,
        main = "PC5 Scores by Treatment Group",
        xlab = "Treatment", ylab = "PC5 Score")

boxplot(pca_scores[, 6] ~ PDX3_Control_and_Vand@meta.data$sampleID,
        main = "PC6 Scores by Treatment Group",
        xlab = "Treatment", ylab = "PC6 Score")
###


# Extract the loadings from the PCA result
pca_loadings <- pca_results$rotation

# For example, inspect the top contributing kinases for PC1:
top_PC1 <- sort(abs(pca_loadings[, "PC1"]), decreasing = TRUE)
top_contributors_PC1 <- head(top_PC1, 10)
print(top_contributors_PC1)



# Option 1: Using the PCA scores matrix
num_PCs <- ncol(pca_results$x)

# Number of PCs to inspect
num_PCs <- 6  # Change this to how many PCs you want to inspect

for (i in 1:num_PCs) {
  pc_name <- paste0("PC", i)
  cat("Top contributors for", pc_name, ":\n")
  top_contributors <- head(sort(abs(pca_loadings[, pc_name]), decreasing = TRUE), 50)
  print(top_contributors)
  cat("\n")
}

# Calculate the variance explained by each PC
var_explained <- (pca_results$sdev)^2
percent_var_explained <- var_explained / sum(var_explained) * 100

# Print the percentage of variance explained by each PC
print(percent_var_explained)



##################### DA peaks

DefaultAssay(filtered_obj)<- "peaks"

da_peaks <- FindMarkers(
  object = filtered_obj,
  ident.1 = 'PDX3_Tam_7_8_14',
  ident.2 = 'PDX3_Control_14',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

da_peaks_control <- FindMarkers(
  object = filtered_obj,
  ident.1 = 'PDX3_Control_14',
  ident.2 = 'PDX3_Tam_7_8_14',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)


# get top differentially accessible peaks
top.da.peak_upr_Tam <- rownames(da_peaks[da_peaks$p_val_adj < 0.1 & da_peaks$pct.1 > 0.05 & da_peaks$avg_log2FC > 0.2, ])

top.da.peak_dw_Tam <- rownames(da_peaks_control[da_peaks_control$p_val_adj < 0.1 & da_peaks_control$pct.1 > 0.05 & da_peaks_control$avg_log2FC > 0.2, ])


# test enrichment
enriched.motifs_dw_Tam <- FindMotifs(
  object = filtered_obj,
  features = top.da.peak_dw_Tam
)

# test enrichment
enriched.motifs_up_Tam <- FindMotifs(
  object = filtered_obj,
  features = top.da.peak_upr_Tam
)

#############################


# Duplicate the RNA assay into a new assay
filtered_obj[["RNA_transformed"]] <- CreateAssayObject(counts = filtered_obj[["RNA"]]@layers[["counts"]])

filtered_obj[["RNA_transformed"]]@data <- filtered_obj[["RNA"]]@data

# Extract the sparse count matrix from the RNA assay
counts_mat <- filtered_obj@assays[["RNA"]]@layers[["counts"]]
# rownames = the features in the RNA assay
rownames(counts_mat) <- rownames(filtered_obj@assays[["RNA"]])

# colnames = the cell names in the overall Seurat object
colnames(counts_mat) <- colnames(filtered_obj)

nrow(counts_mat) == length(rownames(filtered_obj@assays[["RNA"]]))
ncol(counts_mat) == length(colnames(filtered_obj))


new_assay <- CreateAssayObject(counts = counts_mat)

filtered_obj[["RNA_transformed"]] <- new_assay






# If you want a one-time log-normalized matrix:
# 1) extract raw counts
cmat <- GetAssayData(object = filtered_obj, assay = "RNA", layer = "counts") 
# 2) log-normalize on your own or using `NormalizeData()`:
library(Matrix)
# Convert sparse to dense if needed, or keep as sparse if implementing carefully
cmat <- as(as.matrix(cmat), "dgCMatrix") 

# A quick manual log-normalization (pseudo-code):
# sum per cell, scale factor, log1p
sf <- colSums(cmat)
sf[sf == 0] <- 1  # avoid division by zero
norm_mat <- sweep(cmat, 2, sf, "/") * 1e4  # scale factor 10k
lognorm_mat <- log1p(norm_mat)

# Now you have log-normalized expression, reminiscent of what `NormalizeData()` would produce.

# Provide row/col names (if needed):
rownames(lognorm_mat) <- rownames(cmat)
colnames(lognorm_mat) <- colnames(cmat)

# Then store it in a new assay, e.g. "RNA_transformed"
filtered_obj[["RNA_transformed"]] <- CreateAssayObject(counts = cmat)
filtered_obj[["RNA_transformed"]]@data <- lognorm_mat



#################33



# 2) Identify your control cells
ctrl_cells <- WhichCells(filtered_obj, expression = sampleID == "PDX3_Control")


# read.csv for CSV files
df <- read.csv("tam_vector.csv", stringsAsFactors = FALSE)
head(df)


logFC_vector <- setNames(df$Log2FC, df$gene)



####################





##################

control_obj <- subset(filtered_obj, subset = sampleID == "PDX3_Control")

control_obj <- RenameCells(
  object = control_obj,
  new.names = paste0(colnames(control_obj), "_TRANSFORMED")
)


# Step C: Assign them a unique label
control_obj$treatment_label <- "Control_Transformed"

mat <- GetAssayData(
  object = control_obj,
  assay = "RNA",
  layer = "counts"  # The raw counts layer
)



common_genes <- intersect(rownames(mat), names(logFC_vector))

submat <- mat[common_genes, , drop = FALSE]
submat_untr <- mat[common_genes, , drop = FALSE]

str(submat_untr)
gene_index <- which(rownames(submat_untr) == "ESR1")
ESR1_values_before <- submat_untr[gene_index, ]



shift <- 2 ^ logFC_vector[common_genes]

submat <- sweep(submat, MARGIN = 1, STATS = shift, FUN = "*")
submat <- round(submat)
# Optionally clamp negative values if you do some advanced offset 
# (but normally 2^(log2FC) won't be negative).

str(submat)
gene_index <- which(rownames(submat) == "ESR1")
ESR1_values_after_transf <- submat[gene_index, ]


mat[common_genes, ] <- submat

gene_index <- which(mat@Dimnames[[1]] == "ESR1")

ESR1_values_in_mat <- mat[gene_index, ]



# Step A: Create a classic v3/v4-style Assay object
new_assay <- CreateAssayObject(counts = mat)

# Step B: Store it in 'combined_obj' under a new name, e.g. "MERGED"
control_obj[["RNA_trans"]] <- new_assay



combined_obj <- merge(
  x = filtered_obj,
  y = control_obj
)


#############

head(combined_obj)

# 1) Extract partial matrices
mat1 <- GetAssayData(combined_obj, assay = "RNA", layer = "counts.1")
mat2 <- GetAssayData(combined_obj, assay = "RNA_trans", layer = "counts")

# 2) Check row names match. They should be the same gene set, in the same order:
stopifnot(identical(rownames(mat), rownames(mat2)))

# 3) cbind them. The resulting matrix has (genes) x (sum of cells)
merged_counts <- cbind(mat1, mat2)

library(Seurat)

# Step A: Create a classic v3/v4-style Assay object
new_assay <- CreateAssayObject(counts = merged_counts)

# Step B: Store it in 'combined_obj' under a new name, e.g. "MERGED"
combined_obj[["MERGED"]] <- new_assay


combined_obj <- SCTransform(
  combined_obj,
  assay = "MERGED",
  verbose = TRUE,
  new.assay.name = "SCT_MERGED"
)



DefaultAssay(combined_obj) <- "SCT_MERGED"

# Then do PCA, UMAP:
combined_obj <- RunPCA(combined_obj, dims=20)
combined_obj <- RunUMAP(combined_obj, dims = 1:20)

Idents(combined_obj) <- "treatment_label"



DimPlot(combined_obj, reduction = "umap")
DimPlot(combined_obj, reduction = "pca")

Reductions(combined_obj)



combined_obj[["MERGED"]] <- CreateAssayObject(counts = merged_counts)
ctrl_cells <- WhichCells(combined_obj, expression = treatment_label == "Control")
trans_cells <- WhichCells(combined_obj, expression = treatment_label == "Control_Transformed")


exp_mat <- GetAssayData(combined_obj, assay = "MERGED", slot = "counts")
# If you performed 'NormalizeData(MERGED)', you might want slot="data" for log1p expression.



# Check a few example cells:
head( exp_mat["HSPA5", ctrl_cells] )
head( exp_mat["HSPA5", trans_cells] )

# Summaries (mean count, median, etc.)
mean_control <- mean(exp_mat["HSPA5", ctrl_cells])
mean_trans   <- mean(exp_mat["HSPA5", trans_cells])

median_control <- median(exp_mat["HSPA5", ctrl_cells])
median_trans   <- median(exp_mat["HSPA5", trans_cells])

mean_control
mean_trans

head(control_obj)

tail(filtered_obj)
head(control_obj)









########################################33


control_obj <- subset(filtered_obj, subset = sampleID == "PDX3_Control")
control_obj <- RenameCells(
  object = control_obj,
  new.names = paste0(colnames(control_obj), "_TRANSFORMED")
)
control_obj$treatment_label <- "Control_Transformed"

mat <- GetAssayData(
  object = control_obj,
  assay = "RNA",
  layer = "counts"  # The raw counts
)

common_genes <- intersect(rownames(mat), names(logFC_vector))
submat <- mat[common_genes, , drop = FALSE]

# Multiply by 2^(log2FC), round:
shift <- 2 ^ logFC_vector[common_genes]
submat <- sweep(submat, 1, shift, "*")
submat <- round(submat)

# Store back
mat[common_genes, ] <- submat
control_obj[["RNA2"]]@layers[["counts"]] <- mat

# Now 'control_obj' has updated raw counts




exp_mat <- GetAssayData(control_obj, assay = "RNA", slot = "counts")
# Check a few example cells:
head(exp_mat["HSPA5",  ] )
head(exp_mat["HSPA5",  ] )






































suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(FNN)
  library(patchwork)
})

###############################################################################
## SETTINGS
###############################################################################


outdir <- "/path/to/PDX_analysis_R/PDX3_harmony_soft_integration"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


theta_values <- c(0.1, 0.25, 0.5, 1, 2, 4)
npcs_max <- 50
umap_neighbors <- 30
cluster_resolution <- 0.4

sample_col <- "sampleID"
score_col  <- "lowESR1_nonprolif_up"

control_label <- "PDX3_Control"
treated_label <- "PDX3_Tam"

# define "PDX3 cells most similar to the PDX1 low_ESR1_nonprolif program"
top_fraction_query <- 0.10   # top 10% within control
fallback_dims <- 20
k_neighbors <- 30

###############################################################################
## HELPERS
###############################################################################

`%||%` <- function(a, b) if (!is.null(a)) a else b

choose_dims_from_elbow <- function(stdev, max_dims = 30, min_dims = 10) {
  if (length(stdev) < 3) return(min(length(stdev), max_dims))
  
  pct_var <- (stdev^2) / sum(stdev^2)
  drops <- diff(pct_var)
  elbow <- which(abs(drops) < 0.002)[1]
  
  if (is.na(elbow)) elbow <- fallback_dims
  elbow <- max(min_dims, elbow)
  elbow <- min(max_dims, elbow, length(stdev))
  elbow
}

safe_featureplot <- function(obj, feature, reduction_use, title_txt = NULL) {
  if (!feature %in% colnames(obj@meta.data) && !feature %in% rownames(obj)) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = paste("Feature not found:", feature), size = 6) +
        theme_void() +
        ggtitle(title_txt %||% feature)
    )
  }
  
  FeaturePlot(
    obj,
    reduction = reduction_use,
    features = feature,
    raster = FALSE
  ) + ggtitle(title_txt %||% feature)
}

centroid_distance <- function(mat, cells_a, cells_b) {
  cells_a <- intersect(cells_a, rownames(mat))
  cells_b <- intersect(cells_b, rownames(mat))
  if (length(cells_a) < 2 || length(cells_b) < 2) return(NA_real_)
  
  cen_a <- colMeans(mat[cells_a, , drop = FALSE])
  cen_b <- colMeans(mat[cells_b, , drop = FALSE])
  sqrt(sum((cen_a - cen_b)^2))
}

mean_cell_to_group_distance <- function(mat, cells_from, cells_to) {
  cells_from <- intersect(cells_from, rownames(mat))
  cells_to   <- intersect(cells_to, rownames(mat))
  if (length(cells_from) < 2 || length(cells_to) < 2) return(NA_real_)
  
  sub_from <- mat[cells_from, , drop = FALSE]
  sub_to   <- mat[cells_to, , drop = FALSE]
  
  dmat <- as.matrix(dist(rbind(sub_from, sub_to)))
  n1 <- nrow(sub_from)
  n2 <- nrow(sub_to)
  block <- dmat[seq_len(n1), n1 + seq_len(n2), drop = FALSE]
  mean(apply(block, 1, min))
}

nn_composition <- function(mat, meta, query_cells, k = 30, label_col = "sampleID") {
  query_cells <- intersect(query_cells, rownames(mat))
  if (length(query_cells) < 2) return(NULL)
  
  nn <- FNN::get.knnx(data = mat, query = mat[query_cells, , drop = FALSE], k = k + 1)
  idx <- nn$nn.index
  
  nn_cells <- apply(idx, 1, function(x) rownames(mat)[x])
  nn_cells <- lapply(seq_along(query_cells), function(i) setdiff(nn_cells[, i], query_cells[i])[1:k])
  
  out <- lapply(seq_along(query_cells), function(i) {
    neigh <- nn_cells[[i]]
    labs <- meta[neigh, label_col, drop = TRUE]
    tab <- table(labs)
    data.frame(
      cell = query_cells[i],
      label = names(tab),
      n = as.numeric(tab),
      frac = as.numeric(tab) / sum(tab),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(out)
}

###############################################################################
## LOAD OBJECT
###############################################################################

cat("Loading object...\n")
obj <- PDX3_Control_Tam
DefaultAssay(obj) <- "SCT"

if (!sample_col %in% colnames(obj@meta.data)) {
  stop("Metadata column not found: ", sample_col)
}
if (!score_col %in% colnames(obj@meta.data)) {
  stop("Metadata column not found: ", score_col)
}

obj <- subset(obj, subset = .data[[sample_col]] %in% c(control_label, treated_label))

cat("Cells after subset:\n")
print(table(obj@meta.data[[sample_col]]))

###############################################################################
## DEFINE QUERY GROUP FROM SIGNATURE
###############################################################################

meta <- obj@meta.data

control_cells <- rownames(meta)[meta[[sample_col]] == control_label]
treated_cells <- rownames(meta)[meta[[sample_col]] == treated_label]

control_scores <- meta[control_cells, score_col, drop = TRUE]
score_cutoff <- as.numeric(stats::quantile(control_scores, probs = 1 - top_fraction_query, na.rm = TRUE))

obj$query_group_lowESR1_nonprolif_like <- "other"
obj$query_group_lowESR1_nonprolif_like[
  rownames(obj@meta.data) %in% control_cells &
    obj@meta.data[[score_col]] >= score_cutoff
] <- "top_score_control"

obj$query_group_lowESR1_nonprolif_like[
  rownames(obj@meta.data) %in% treated_cells
] <- "treated"

cells_query <- rownames(obj@meta.data)[
  obj$query_group_lowESR1_nonprolif_like == "top_score_control"
]

cells_other_ctrl <- rownames(obj@meta.data)[
  obj$query_group_lowESR1_nonprolif_like == "other" &
    obj@meta.data[[sample_col]] == control_label
]

cells_treated <- rownames(obj@meta.data)[
  obj@meta.data[[sample_col]] == treated_label
]

cat("Query definition based on ", score_col, "\n", sep = "")
cat("Top fraction in control = ", top_fraction_query, "\n", sep = "")
cat("Score cutoff = ", signif(score_cutoff, 4), "\n", sep = "")
cat("n query cells = ", length(cells_query), "\n", sep = "")
cat("n other control = ", length(cells_other_ctrl), "\n", sep = "")
cat("n treated = ", length(cells_treated), "\n", sep = "")

fwrite(
  data.frame(
    top_fraction_query = top_fraction_query,
    score_cutoff = score_cutoff,
    n_query = length(cells_query),
    n_other_control = length(cells_other_ctrl),
    n_treated = length(cells_treated)
  ),
  file = file.path(outdir, "00_query_definition.tsv"),
  sep = "\t"
)

###############################################################################
## PCA + ELBOW
###############################################################################

cat("Running PCA...\n")
obj <- RunPCA(obj, npcs = npcs_max, verbose = FALSE)

pdf(file.path(outdir, "01_ElbowPlot_PDX3_Control_Tam.pdf"), width = 7, height = 5)
print(ElbowPlot(obj, ndims = npcs_max))
dev.off()

pca_stdev <- Stdev(obj, reduction = "pca")
dims_use_n <- choose_dims_from_elbow(pca_stdev, max_dims = 30, min_dims = 10)
dims_use <- 1:dims_use_n

writeLines(
  paste0(
    "Chosen dims from elbow heuristic: ", dims_use_n, "\n",
    "Dims used: ", paste(dims_use, collapse = ","), "\n"
  ),
  con = file.path(outdir, "01_dims_used.txt")
)

###############################################################################
## PRE-HARMONY UMAP
###############################################################################

obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap_preHarmony",
  reduction.key = "UMAPPRE_",
  verbose = FALSE
)

p_pre1 <- DimPlot(
  obj, reduction = "umap_preHarmony", group.by = sample_col, raster = FALSE
) + ggtitle("Pre-Harmony: sampleID")

p_pre2 <- safe_featureplot(
  obj, "ESR1", "umap_preHarmony", "Pre-Harmony: ESR1"
)

p_pre3 <- safe_featureplot(
  obj, score_col, "umap_preHarmony", paste0("Pre-Harmony: ", score_col)
)

p_pre4 <- DimPlot(
  obj, reduction = "umap_preHarmony", group.by = "query_group_lowESR1_nonprolif_like", raster = FALSE
) + ggtitle("Pre-Harmony: query group")

pdf(file.path(outdir, "02_PreHarmony_UMAPs.pdf"), width = 15, height = 10)
print((p_pre1 + p_pre2) / (p_pre3 + p_pre4))
dev.off()

###############################################################################
## HARMONY RUNS
###############################################################################

distance_summary <- list()
nn_summary <- list()

for (theta in theta_values) {
  cat("\n============================\n")
  cat("Running Harmony, theta = ", theta, "\n", sep = "")
  cat("============================\n")
  
  obj_h <- obj
  
  harm_name <- paste0("harmony_theta_", gsub("\\.", "_", as.character(theta)))
  umap_name <- paste0("umap_harmony_theta_", gsub("\\.", "_", as.character(theta)))
  
  obj_h <- RunHarmony(
    object = obj_h,
    group.by.vars = sample_col,
    reduction = "pca",
    assay.use = "SCT",
    dims.use = dims_use,
    theta = theta,
    lambda = 1,
    plot_convergence = TRUE,
    reduction.save = harm_name,
    verbose = TRUE
  )
  
  obj_h <- RunUMAP(
    obj_h,
    reduction = harm_name,
    dims = dims_use,
    reduction.name = umap_name,
    reduction.key = paste0("UMAPTH", gsub("\\.", "", as.character(theta)), "_"),
    n.neighbors = umap_neighbors,
    verbose = FALSE
  )
  
  obj_h <- FindNeighbors(
    obj_h,
    reduction = harm_name,
    dims = dims_use,
    verbose = FALSE
  )
  
  obj_h <- FindClusters(
    obj_h,
    resolution = cluster_resolution,
    verbose = FALSE
  )
  
  saveRDS(
    obj_h,
    file = file.path(outdir, paste0("PDX3_Control_Tam_", harm_name, ".rds"))
  )
  
  p1 <- DimPlot(
    obj_h,
    reduction = umap_name,
    group.by = sample_col,
    raster = FALSE
  ) + ggtitle(paste0("sampleID | theta=", theta))
  
  p2 <- safe_featureplot(
    obj_h, "ESR1", umap_name, paste0("ESR1 | theta=", theta)
  )
  
  p3 <- safe_featureplot(
    obj_h, score_col, umap_name, paste0(score_col, " | theta=", theta)
  )
  
  p4 <- DimPlot(
    obj_h,
    reduction = umap_name,
    group.by = "query_group_lowESR1_nonprolif_like",
    raster = FALSE
  ) + ggtitle(paste0("query group | theta=", theta))
  
  pdf(file.path(outdir, paste0("UMAP_theta_", gsub("\\.", "_", theta), ".pdf")), width = 15, height = 10)
  print((p1 + p2) / (p3 + p4))
  dev.off()
  
  emb <- Embeddings(obj_h, reduction = harm_name)[, dims_use, drop = FALSE]
  
  d_centroid_query_treated <- centroid_distance(emb, cells_query, cells_treated)
  d_centroid_otherctrl_treated <- centroid_distance(emb, cells_other_ctrl, cells_treated)
  
  d_min_query_treated <- mean_cell_to_group_distance(emb, cells_query, cells_treated)
  d_min_otherctrl_treated <- mean_cell_to_group_distance(emb, cells_other_ctrl, cells_treated)
  
  distance_summary[[as.character(theta)]] <- data.frame(
    theta = theta,
    n_dims = length(dims_use),
    n_query = length(intersect(cells_query, rownames(emb))),
    n_other_ctrl = length(intersect(cells_other_ctrl, rownames(emb))),
    n_treated = length(intersect(cells_treated, rownames(emb))),
    centroid_dist_query_to_treated = d_centroid_query_treated,
    centroid_dist_otherctrl_to_treated = d_centroid_otherctrl_treated,
    mean_minDist_query_to_treated = d_min_query_treated,
    mean_minDist_otherctrl_to_treated = d_min_otherctrl_treated,
    delta_centroid = d_centroid_query_treated - d_centroid_otherctrl_treated,
    delta_mean_minDist = d_min_query_treated - d_min_otherctrl_treated,
    stringsAsFactors = FALSE
  )
  
  nn_df <- nn_composition(
    mat = emb,
    meta = obj_h@meta.data,
    query_cells = cells_query,
    k = k_neighbors,
    label_col = sample_col
  )
  
  if (!is.null(nn_df)) {
    nn_df$theta <- theta
    fwrite(
      nn_df,
      file = file.path(outdir, paste0("NN_composition_theta_", gsub("\\.", "_", theta), ".tsv")),
      sep = "\t"
    )
    
    nn_summary[[as.character(theta)]] <- nn_df %>%
      group_by(theta, label) %>%
      summarise(
        mean_frac = mean(frac, na.rm = TRUE),
        mean_n = mean(n, na.rm = TRUE),
        .groups = "drop"
      )
  }
}

###############################################################################
## SAVE SUMMARIES
###############################################################################

distance_summary_df <- bind_rows(distance_summary)
fwrite(
  distance_summary_df,
  file = file.path(outdir, "distance_summary_across_thetas.tsv"),
  sep = "\t"
)

if (length(nn_summary) > 0) {
  nn_summary_df <- bind_rows(nn_summary)
  fwrite(
    nn_summary_df,
    file = file.path(outdir, "NN_summary_across_thetas.tsv"),
    sep = "\t"
  )
}

writeLines(
  c(
    "Interpretation guide:",
    "",
    "Query cells = top-scoring PDX3_Control cells for lowESR1_nonprolif_up.",
    "",
    "If delta_centroid < 0, the high-score control cells are closer to treated cells",
    "than the rest of control cells in Harmony space.",
    "",
    "If delta_mean_minDist < 0, the high-score control cells have treated cells as",
    "closer neighbors than the rest of control cells.",
    "",
    "In NN_summary_across_thetas.tsv, higher mean_frac for PDX3_Tam among neighbors",
    "of query cells supports that these cells sit closer to treated cells.",
    "",
    "Compare UMAPs across theta values and prefer the smallest theta that improves",
    "mixing without erasing structure."
  ),
  con = file.path(outdir, "README_interpretation.txt")
)

cat("\nDone.\nOutputs written to:\n", outdir, "\n", sep = "")


































suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(FNN)
  library(data.table)
  library(patchwork)
})

###############################################################################
## SETTINGS
###############################################################################

DefaultAssay(PDX3_Control_Tam) <- "SCT"

dims_use <- 1:10   # restricted PCA usually needs fewer PCs; start with 10
score_col <- "lowESR1_nonprolif_difference"
sample_col <- "sampleID"

###############################################################################
## 1. BUILD RESTRICTED GENE SET
###############################################################################

genes_state_up <- unique(c(
  High_ESR1_up,
  non_prolif_up,
  prolif_up
))

# keep only genes present in PDX3 object
genes_state_up <- intersect(genes_state_up, rownames(PDX3_Control_Tam))

length(genes_state_up)
head(genes_state_up)

###############################################################################
## 2. SCALE ONLY THOSE GENES AND RUN PCA
###############################################################################

# keep object as-is, just create a new PCA reduction from restricted genes
PDX3_Control_Tam <- ScaleData(
  PDX3_Control_Tam,
  features = genes_state_up,
  verbose = FALSE
)

PDX3_Control_Tam <- RunPCA(
  PDX3_Control_Tam,
  features = genes_state_up,
  npcs = min(30, length(genes_state_up) - 1),
  reduction.name = "pca_state_up",
  verbose = FALSE
)

ElbowPlot(PDX3_Control_Tam, reduction = "pca_state_up", ndims = 30)

###############################################################################
## 3. BUILD UMAP FROM THE RESTRICTED PCA
###############################################################################

PDX3_Control_Tam <- FindNeighbors(
  PDX3_Control_Tam,
  reduction = "pca_state_up",
  dims = dims_use,
  graph.name = "stateup_nn",
  verbose = FALSE
)

PDX3_Control_Tam <- FindClusters(
  PDX3_Control_Tam,
  graph.name = "stateup_nn",
  resolution = 0.8,
  algorithm = 1,
  verbose = FALSE
)

PDX3_Control_Tam <- RunUMAP(
  PDX3_Control_Tam,
  reduction = "pca_state_up",
  dims = dims_use,
  reduction.name = "umap_state_up",
  reduction.key = "UMAPSTATEUP_",
  verbose = FALSE
)

###############################################################################
## 4. VISUALIZE
###############################################################################

p1 <- DimPlot(
  PDX3_Control_Tam,
  reduction = "umap_state_up",
  group.by = "sampleID",
  raster = FALSE
) + ggtitle("Restricted PCA/UMAP | sampleID")

p2 <- FeaturePlot(
  PDX3_Control_Tam,
  reduction = "umap_state_up",
  features = "ESR1",
  raster = FALSE
) + ggtitle("Restricted PCA/UMAP | ESR1")

p3 <- FeaturePlot(
  PDX3_Control_Tam,
  reduction = "umap_state_up",
  features = "lowESR1_nonprolif_difference",
  raster = FALSE
) + ggtitle("Restricted PCA/UMAP | lowESR1_nonprolif_difference")

p4 <- FeaturePlot(
  PDX3_Control_Tam,
  reduction = "umap_state_up",
  features = "lowESR1_prolif_difference",
  raster = FALSE
) + ggtitle("Restricted PCA/UMAP | lowESR1_prolif_difference")

print((p1 + p2) / (p3 + p4))



colnames(Embeddings(PDX3_Control_Tam, "umap_state_up"))
Reductions(PDX3_Control_Tam)

# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(PDX3_Control_Tam, 
                       vars = c("lowESR1_nonprolif_difference", "UMAPSTATEUP_1", "UMAPSTATEUP_2"))


# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = UMAPSTATEUP_1, y = UMAPSTATEUP_2, color = lowESR1_nonprolif_difference)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(-0.1,0.4),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "lowESR1_nonprolif_difference") +
  theme_minimal()






# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(PDX3_Control_Tam, 
                       vars = c("ESR1", "UMAPSTATEUP_1", "UMAPSTATEUP_2"))


# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = UMAPSTATEUP_1, y = UMAPSTATEUP_2, color = ESR1)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(-0.1,3.5),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "ESR1") +
  theme_minimal()




# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(PDX3_Control_Tam, 
                       vars = c("lowESR1_prolif_difference", "UMAPSTATEUP_1", "UMAPSTATEUP_2"))


# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = UMAPSTATEUP_1, y = UMAPSTATEUP_2, color = lowESR1_prolif_difference)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(-0.1,0.20),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "lowESR1_prolif_difference") +
  theme_minimal()








cds <- as.cell_data_set(PDX3_Control_Tam)

# Suppose you've already set up or transferred your UMAP from Seurat
reducedDims(cds)$UMAP <- Embeddings(PDX3_Control_Tam, "umap_state_up")

# Optionally cluster
cds <- cluster_cells(cds)  # uses the existing "UMAP" dims

# Learn graph (no "reduction_method" argument needed)
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# If you do want pseudotime
cds <- order_cells(cds)

plot_cells(cds, 
           color_cells_by = "pseudotime", 
           cell_size = 1,
           show_trajectory_graph = TRUE)

?plot_cells
