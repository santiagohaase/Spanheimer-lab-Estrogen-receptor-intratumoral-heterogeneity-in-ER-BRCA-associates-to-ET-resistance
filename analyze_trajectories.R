#### Analyze trajectories

###################### Trajectory analysis in control and treated cells

library(Seurat)
library(monocle3)
library(ggplot2)
library(BiocParallel)
library(UCell)

###############################################################################
# 1. Load merged object and subset samples of interest
###############################################################################

PDX_combined <- readRDS("/path/to/combined_multiome_object.rds")

DimPlot(PDX_combined)
levels(Idents(PDX_combined))

# Assumes identities in PDX_combined correspond to sample identities
PDX3_Control_Tam <- subset(PDX_combined, idents = c("PDX3_Control", "PDX3_Tam"))

###############################################################################
# 2. Clean metadata if needed
###############################################################################

meta_cols <- colnames(PDX3_Control_Tam@meta.data)

# Remove metadata columns where "kinase" appears twice in the column name
duplicate_cols <- grep(".*kinase.*kinase.*", meta_cols, value = TRUE)

if (length(duplicate_cols) > 0) {
  PDX3_Control_Tam@meta.data <- PDX3_Control_Tam@meta.data[
    , !colnames(PDX3_Control_Tam@meta.data) %in% duplicate_cols, drop = FALSE
  ]
}

# Remove selected UCell metadata columns before recalculating, only if present
cols_to_remove <- c(
  "ESR1_human_tf_ARCHS4_coexpression_UCell",
  "GATA3_human_tf_ARCHS4_coexpression_UCell",
  "RET_human_kinase_ARCHS4_coexpression_UCell"
)

existing_remove <- intersect(cols_to_remove, colnames(PDX3_Control_Tam@meta.data))
if (length(existing_remove) > 0) {
  PDX3_Control_Tam@meta.data[, existing_remove] <- NULL
}

###############################################################################
# 3. Add missing UCell signatures
###############################################################################
# The object "signatures" must exist in the environment before this step.
# It should be a named list where each element is a vector of gene symbols.

if (!exists("signatures")) {
  stop("The object 'signatures' is not defined. Load or define it before running AddModuleScore_UCell().")
}

bp <- MulticoreParam(workers = 32)

PDX3_Control_Tam <- AddModuleScore_UCell(
  obj = PDX3_Control_Tam,
  features = signatures,
  assay = "SCT",
  slot = "counts",
  BPPARAM = bp,
  name = NULL
)

bpstop(bp)

saveRDS(
  PDX3_Control_Tam,
  "/path/to/output/PDX3_Control_Tam.RDS"
)

###############################################################################
# 4. Dimensional reduction, clustering, and UMAP
###############################################################################

DefaultAssay(PDX3_Control_Tam) <- "SCT"

PDX3_Control_Tam <- RunPCA(
  PDX3_Control_Tam,
  features = VariableFeatures(PDX3_Control_Tam)
)

ElbowPlot(PDX3_Control_Tam)

PDX3_Control_Tam <- FindNeighbors(PDX3_Control_Tam, dims = 1:13)
PDX3_Control_Tam <- FindClusters(PDX3_Control_Tam, resolution = 1)
PDX3_Control_Tam <- RunUMAP(PDX3_Control_Tam, dims = 1:13)

DimPlot(PDX3_Control_Tam, reduction = "umap", label = TRUE)

###############################################################################
# 5. Control-only subset for signature/expression visualization
###############################################################################

# Subset by metadata rather than identities, because identities may change
PDX3_Control <- subset(PDX3_Control_Tam, subset = sampleID == "PDX3_Control")

Idents(PDX3_Control_Tam) <- "seurat_clusters"
DimPlot(PDX3_Control_Tam, reduction = "umap", label = TRUE)

Idents(PDX3_Control) <- "seurat_clusters"
DimPlot(PDX3_Control, reduction = "umap", label = TRUE)

Idents(PDX3_Control_Tam) <- "sampleID"
DimPlot(PDX3_Control_Tam, reduction = "umap", label = TRUE)

###############################################################################
# 6. UMAP plots colored by metadata/signatures in control cells
###############################################################################

umap_data <- FetchData(
  PDX3_Control,
  vars = c("lowESR1_nonprolif_difference", "umap_1", "umap_2")
)

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = lowESR1_nonprolif_difference)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_viridis_c(
    option = "C",
    limits = c(-0.1, 0.2),
    na.value = "lightgrey",
    oob = scales::squish
  ) +
  labs(color = "lowESR1_nonprolif_difference") +
  theme_minimal()

umap_data <- FetchData(
  PDX3_Control,
  vars = c("ESR1", "umap_1", "umap_2")
)

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = ESR1)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_viridis_c(
    option = "C",
    limits = c(0, 3.5),
    na.value = "lightgrey",
    oob = scales::squish
  ) +
  labs(color = "ESR1 Expression") +
  theme_minimal()

umap_data <- FetchData(
  PDX3_Control,
  vars = c("lowESR1_prolif_difference", "umap_1", "umap_2")
)

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = lowESR1_prolif_difference)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_viridis_c(
    option = "C",
    limits = c(-0.1, 0.2),
    na.value = "lightgrey",
    oob = scales::squish
  ) +
  labs(color = "lowESR1_prolif_difference") +
  theme_minimal()

###############################################################################
# 7. Monocle3 trajectory and pseudotime
###############################################################################

cds <- as.cell_data_set(PDX3_Control_Tam)

# Transfer Seurat UMAP coordinates to Monocle3
reducedDims(cds)$UMAP <- Embeddings(PDX3_Control_Tam, "umap")

# Cluster cells in Monocle3 using the transferred UMAP
cds <- cluster_cells(cds)

# Learn principal graph
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# Order cells to compute pseudotime
# IMPORTANT:
# When this command opens the interactive plot, manually select the root node
# corresponding to the biologically earliest cell population, as done in the
# original analysis.
cds <- order_cells(cds)

plot_cells(
  cds,
  color_cells_by = "pseudotime",
  cell_size = 1,
  trajectory_graph_segment_size = 2,
  show_trajectory_graph = TRUE
)