## Recurrent samples + PDX analysis 

## =========================================================
## Replace these placeholders with local paths before running.
project_dir <- "/path/to/project_root"
data_dir <- file.path(project_dir, "data")
pdx_dir <- file.path(data_dir, "PDX_analysis_R")
multiome_dir <- file.path(data_dir, "Multiome_tumors")
fragment_dir <- file.path(data_dir, "fragments")
requant_dir <- file.path(project_dir, "results", "requant")
uncover_dir <- file.path(project_dir, "results", "UNCover_recurrent_vs_PDX_lowESR1")
chip_catalog_bed <- file.path(project_dir, "resources", "remap_tf_4col.bed")
conda_sh <- "/path/to/miniconda3/etc/profile.d/conda.sh"
python_bin <- "/path/to/miniconda3/envs/myenv/bin/python"

dir.create(requant_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(uncover_dir, recursive = TRUE, showWarnings = FALSE)

library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)

PDX3_Control_Tam <- readRDS(
  file.path(pdx_dir, "PDX3_Control_Tam.RDS")
)

PDX1_Control_Tam <- readRDS(
  file.path(pdx_dir, "PDX1_Control_Tam.RDS")
)

Idents(PDX3_Control_Tam) <- 'sampleID'
PDX3_Control <- subset(PDX3_Control_Tam, idents = 'PDX3_Control')

Idents(PDX1_Control_Tam) <- 'sampleID'
PDX1_Control <- subset(PDX1_Control_Tam, idents = 'PDX1_Control')


Tumor_5572CL_epithelial <- readRDS(file.path(multiome_dir, "Tumor_5572CL_epithelial.RDS"))
Tumor_52BC3L_epithelial <- readRDS(file.path(multiome_dir, "Tumor_52BC3L_epithelial.RDS"))
Tumor_644D9L_epithelial <- readRDS(file.path(multiome_dir, "Tumor_644D9L_epithelial.RDS"))
Tumor_73607L_epithelial <- readRDS(file.path(multiome_dir, "Tumor_73607L_epithelial.RDS"))


## =========================================================
## 1. Helpers
## =========================================================

safe_assays <- function(obj) {
  tryCatch(names(obj@assays), error = function(e) character(0))
}

safe_reductions <- function(obj) {
  tryCatch(names(obj@reductions), error = function(e) character(0))
}

rename_assay_if_present <- function(obj, old_name, new_name) {
  assays_now <- safe_assays(obj)
  if (old_name %in% assays_now && !(new_name %in% assays_now)) {
    obj[[new_name]] <- obj[[old_name]]
    obj[[old_name]] <- NULL
  }
  obj
}

standardize_multiome_object <- function(obj,
                                        sample_name,
                                        source_group = c("PDX", "Recurrent"),
                                        treatment = "Untreated") {
  source_group <- match.arg(source_group)
  
  message("---- Standardizing ", sample_name, " ----")
  message("Assays before: ", paste(safe_assays(obj), collapse = ", "))
  message("Reductions before: ", paste(safe_reductions(obj), collapse = ", "))
  
  ## rename chromatin assay to ATAC if needed
  obj <- rename_assay_if_present(obj, old_name = "peaks", new_name = "ATAC")
  
  assays_now <- safe_assays(obj)
  
  ## minimal checks
  if (!("RNA" %in% assays_now) && !("SCT" %in% assays_now)) {
    stop(sample_name, ": no RNA or SCT assay found.")
  }
  if (!("ATAC" %in% assays_now)) {
    warning(sample_name, ": no ATAC assay found after standardization.")
  }
  
  ## metadata
  obj$sampleID <- sample_name
  obj$dataset  <- source_group
  obj$treatment_group <- treatment
  
  obj$tumor_class <- ifelse(source_group == "PDX", "PDX_untreated", "Recurrent_tumor")
  obj$model_group <- ifelse(grepl("^PDX", sample_name), "PDX", "Tumor")
  
  ## prefix cell names so nothing collides
  obj <- RenameCells(obj, add.cell.id = sample_name)
  
  ## choose a sane default assay
  if ("SCT" %in% safe_assays(obj)) {
    DefaultAssay(obj) <- "SCT"
  } else {
    DefaultAssay(obj) <- "RNA"
  }
  
  message("Assays after: ", paste(safe_assays(obj), collapse = ", "))
  message("Cells: ", ncol(obj), " | Features: ", nrow(obj))
  obj
}

## =========================================================
## 2. Standardize all objects
## =========================================================

PDX1_Control <- standardize_multiome_object(
  PDX1_Control,
  sample_name  = "PDX1_Control",
  source_group = "PDX",
  treatment    = "Untreated"
)

PDX3_Control <- standardize_multiome_object(
  PDX3_Control,
  sample_name  = "PDX3_Control",
  source_group = "PDX",
  treatment    = "Untreated"
)

Tumor_5572CL_epithelial <- standardize_multiome_object(
  Tumor_5572CL_epithelial,
  sample_name  = "Tumor_5572CL",
  source_group = "Recurrent",
  treatment    = "Recurrent"
)

Tumor_52BC3L_epithelial <- standardize_multiome_object(
  Tumor_52BC3L_epithelial,
  sample_name  = "Tumor_52BC3L",
  source_group = "Recurrent",
  treatment    = "Recurrent"
)

Tumor_644D9L_epithelial <- standardize_multiome_object(
  Tumor_644D9L_epithelial,
  sample_name  = "Tumor_644D9L",
  source_group = "Recurrent",
  treatment    = "Recurrent"
)

Tumor_73607L_epithelial <- standardize_multiome_object(
  Tumor_73607L_epithelial,
  sample_name  = "Tumor_73607L",
  source_group = "Recurrent",
  treatment    = "Recurrent"
)

obj_list <- list(
  PDX1_Control,
  PDX3_Control,
  Tumor_5572CL_epithelial,
  Tumor_52BC3L_epithelial,
  Tumor_644D9L_epithelial,
  Tumor_73607L_epithelial
)

names(obj_list) <- c(
  "PDX1_Control",
  "PDX3_Control",
  "Tumor_5572CL",
  "Tumor_52BC3L",
  "Tumor_644D9L",
  "Tumor_73607L"
)

## quick assay summary
assay_summary <- lapply(names(obj_list), function(nm) {
  x <- obj_list[[nm]]
  data.frame(
    sample = nm,
    n_cells = ncol(x),
    assays = paste(safe_assays(x), collapse = ";"),
    reductions = paste(safe_reductions(x), collapse = ";"),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

print(assay_summary)


suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

safe_assays <- function(obj) {
  tryCatch(names(obj@assays), error = function(e) character(0))
}

safe_reductions <- function(obj) {
  tryCatch(names(obj@reductions), error = function(e) character(0))
}

strip_fragments_from_chrom_assay <- function(obj, assay_name = "ATAC") {
  if (!(assay_name %in% safe_assays(obj))) return(obj)
  
  message("Stripping fragments from assay: ", assay_name)
  
  ## safest route: clear fragments slot directly
  try({
    slot(obj[[assay_name]], "fragments") <- list()
  }, silent = TRUE)
  
  ## also clear misc fragment refs if present
  try({
    obj[[assay_name]]@misc$fragments <- NULL
  }, silent = TRUE)
  
  obj
}

drop_redundant_assays <- function(obj) {
  assays_now <- safe_assays(obj)
  
  ## if both ATAC and peaks exist, keep ATAC only
  if ("ATAC" %in% assays_now && "peaks" %in% assays_now) {
    message("Dropping redundant assay: peaks")
    obj[["peaks"]] <- NULL
  }
  
  obj
}

standardize_for_merge <- function(obj, sample_name) {
  message("---- Preparing ", sample_name, " ----")
  message("Assays before: ", paste(safe_assays(obj), collapse = ", "))
  
  ## if object has peaks but no ATAC, rename peaks -> ATAC
  assays_now <- safe_assays(obj)
  if ("peaks" %in% assays_now && !("ATAC" %in% assays_now)) {
    obj[["ATAC"]] <- obj[["peaks"]]
    obj[["peaks"]] <- NULL
  }
  
  ## remove redundant peaks assay if both exist
  obj <- drop_redundant_assays(obj)
  
  ## strip fragments from chromatin assays that may trigger merge failure
  if ("ATAC" %in% safe_assays(obj)) {
    obj <- strip_fragments_from_chrom_assay(obj, "ATAC")
  }
  if ("peaks_union" %in% safe_assays(obj)) {
    obj <- strip_fragments_from_chrom_assay(obj, "peaks_union")
  }
  
  ## metadata
  obj$sampleID <- sample_name
  
  message("Assays after: ", paste(safe_assays(obj), collapse = ", "))
  obj
}

## --------------------------------------------------
## Apply to your objects
## --------------------------------------------------

PDX1_Control <- standardize_for_merge(PDX1_Control, "PDX1_Control")
PDX3_Control <- standardize_for_merge(PDX3_Control, "PDX3_Control")
Tumor_5572CL_epithelial <- standardize_for_merge(Tumor_5572CL_epithelial, "Tumor_5572CL")
Tumor_52BC3L_epithelial <- standardize_for_merge(Tumor_52BC3L_epithelial, "Tumor_52BC3L")
Tumor_644D9L_epithelial <- standardize_for_merge(Tumor_644D9L_epithelial, "Tumor_644D9L")
Tumor_73607L_epithelial <- standardize_for_merge(Tumor_73607L_epithelial, "Tumor_73607L")

obj_list <- list(
  PDX1_Control,
  PDX3_Control,
  Tumor_5572CL_epithelial,
  Tumor_52BC3L_epithelial,
  Tumor_644D9L_epithelial,
  Tumor_73607L_epithelial
)
names(obj_list) <- c(
  "PDX1_Control",
  "PDX3_Control",
  "Tumor_5572CL",
  "Tumor_52BC3L",
  "Tumor_644D9L",
  "Tumor_73607L"
)

## retry merge
merged_multiome <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  merge.data = TRUE
)

table(merged_multiome$sampleID)


## =========================================================
## 4. RNA-only diagnostic workflow WITHOUT integration
## =========================================================

## for initial diagnostics, use RNA or SCT only
if ("SCT" %in% safe_assays(merged_multiome)) {
  DefaultAssay(merged_multiome) <- "SCT"
} else {
  DefaultAssay(merged_multiome) <- "RNA"
}

## if variable features are missing after merge, recompute on RNA
if (length(VariableFeatures(merged_multiome)) < 200) {
  if ("RNA" %in% safe_assays(merged_multiome)) {
    DefaultAssay(merged_multiome) <- "RNA"
    merged_multiome <- FindVariableFeatures(
      merged_multiome,
      selection.method = "vst",
      nfeatures = 3000
    )
    merged_multiome <- ScaleData(merged_multiome, verbose = FALSE)
    merged_multiome <- RunPCA(merged_multiome, npcs = 50, verbose = FALSE)
  } else {
    merged_multiome <- RunPCA(merged_multiome, npcs = 50, verbose = FALSE)
  }
} else {
  merged_multiome <- RunPCA(merged_multiome, npcs = 50, verbose = FALSE)
}


DefaultAssay(merged_multiome) <- "RNA"

## Seurat v5: merge leaves multiple layers; join them first
merged_multiome[["RNA"]] <- JoinLayers(merged_multiome[["RNA"]])

## now create a unified normalized data layer
merged_multiome <- NormalizeData(merged_multiome, verbose = FALSE)

## variable features + scaling + PCA
merged_multiome <- FindVariableFeatures(
  merged_multiome,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

merged_union <- ScaleData(
  merged_union,
  features = VariableFeatures(merged_union),
  verbose = FALSE
)

merged_union <- RunPCA(
  merged_union,
  features = VariableFeatures(merged_union),
  npcs = 50,
  verbose = FALSE
)

merged_multiome <- RunUMAP(
  merged_multiome,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "umap_rna_unintegrated",
  verbose = FALSE
)

merged_multiome <- FindNeighbors(
  merged_multiome,
  dims = 1:30,
  reduction = "pca",
  verbose = FALSE
)

merged_multiome <- FindClusters(
  merged_multiome,
  resolution = 0.4,
  verbose = FALSE
)


merged_multiome <- RunUMAP(
  merged_multiome,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "umap_rna_unintegrated",
  verbose = FALSE
)

merged_multiome <- FindNeighbors(
  merged_multiome,
  dims = 1:30,
  reduction = "pca",
  verbose = FALSE
)

merged_multiome <- FindClusters(
  merged_multiome,
  resolution = 0.4,
  verbose = FALSE
)


## =========================================================
## 5. Diagnostics to decide whether integration is needed
## =========================================================

p1 <- DimPlot(
  merged_multiome,
  reduction = "umap_rna_unintegrated",
  group.by = "sampleID",
  label = TRUE,
  repel = TRUE
) + ggtitle("Unintegrated RNA UMAP by sample")

p2 <- DimPlot(
  merged_multiome,
  reduction = "umap_rna_unintegrated",
  group.by = "tumor_class",
  label = TRUE,
  repel = TRUE
) + ggtitle("Unintegrated RNA UMAP by tumor class")

p3 <- FeaturePlot(
  merged_multiome,
  reduction = "umap_rna_unintegrated",
  features = "ESR1",
  order = TRUE
) + ggtitle("ESR1 on unintegrated RNA UMAP")

print(p1 + p2)
print(p3)


## sample composition across clusters
cluster_sample_table <- table(
  cluster = Idents(merged_multiome),
  sample = merged_multiome$sampleID
)
print(cluster_sample_table)
print(prop.table(cluster_sample_table, margin = 1))

## ESR1 per sample
VlnPlot(
  merged_multiome,
  features = "ESR1",
  group.by = "sampleID",
  pt.size = 0
) + RotatedAxis()

## =====


DefaultAssay(merged_multiome) <- "RNA"
merged_multiome[["RNA"]] <- JoinLayers(merged_multiome[["RNA"]])

merged_multiome <- NormalizeData(merged_multiome, verbose = FALSE)
merged_multiome <- FindVariableFeatures(
  merged_multiome,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)
merged_union <- ScaleData(
  merged_union,
  features = VariableFeatures(merged_union),
  verbose = FALSE
)
merged_union <- RunPCA(
  merged_union,
  features = VariableFeatures(merged_union),
  npcs = 50,
  verbose = FALSE
)

library(harmony)

merged_multiome <- RunHarmony(
  object = merged_multiome,
  group.by.vars = "sampleID",
  reduction.use = "pca",
  dims.use = 1:30,
  project.dim = FALSE
)

merged_multiome <- RunUMAP(
  merged_multiome,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony_rna",
  verbose = FALSE
)


merged_multiome <- FindNeighbors(
  merged_multiome,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

merged_multiome <- FindClusters(
  merged_multiome,
  resolution = 0.4,
  verbose = FALSE
)


DimPlot(merged_multiome, reduction = "umap_harmony_rna", group.by = "sampleID", label = TRUE)
DimPlot(merged_multiome, reduction = "umap_harmony_rna", group.by = "tumor_class", label = TRUE)
FeaturePlot(merged_multiome, reduction = "umap_harmony_rna", features = "ESR1", order = TRUE)

tab_harmony <- table(Idents(merged_multiome), merged_multiome$sampleID)
prop.table(tab_harmony, 1)


library(future)

options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB

DefaultAssay(merged_multiome) <- "RNA"
merged_multiome[["RNA"]] <- JoinLayers(merged_multiome[["RNA"]])

merged_multiome <- SCTransform(
  merged_multiome,
  assay = "RNA",
  new.assay.name = "SCT_recalc",
  verbose = TRUE
)


PDX1_Control_Tam <- readRDS(
  file.path(pdx_dir, "PDX1_Control_Tam.RDS")
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


DefaultAssay(merged_multiome) <- 'SCT'

merged_multiome$lowESR1_nonprolif_difference <- NULL
merged_multiome$lowESR1_nonprolif_up <- NULL
merged_multiome$lowESR1_nonprolif_dw <- NULL
merged_multiome$lowESR1_prolif_difference <- NULL
merged_multiome$lowESR1_prolif_dw <- NULL
merged_multiome$lowESR1_prolif_up <- NULL


merged_multiome$lowESR1_nonprolif_difference

merged_multiome <- AddModuleScore_UCell(merged_multiome, features = list(
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
merged_multiome$lowESR1_nonprolif_difference <- merged_multiome[[up_np]] - merged_multiome[[down_np]]
merged_multiome$lowESR1_prolif_difference    <- merged_multiome[[up_pr]] - merged_multiome[[down_pr]]


# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(merged_multiome, 
                       vars = c("lowESR1_nonprolif_difference", "umapharmonyrna_1", "umapharmonyrna_2"))


# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = umapharmonyrna_1, y = umapharmonyrna_2, color = lowESR1_nonprolif_difference)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(-0.1,0.4),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "lowESR1_nonprolif_difference") +
  theme_minimal()
 
DefaultAssay(merged_multiome)<- 'SCT'
# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(merged_multiome, 
                       vars = c("MUC6", "umapharmonyrna_1", "umapharmonyrna_2"))


# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = umapharmonyrna_1, y = umapharmonyrna_2, color = MUC6)) +
  geom_point(alpha = 0.5, size = 1) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(0.1,3),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "MUC6") +
  theme_minimal()


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
})

## --------------------------------------------------
## 1. Keep samples of interest
## --------------------------------------------------

cells_keep <- colnames(merged_multiome)[
  merged_multiome$sampleID %in% c(
    "PDX1_Control",
    "PDX3_Control",
    "Tumor_52BC3L",
    "Tumor_644D9L",
    "Tumor_73607L"
  )
]

obj_sub <- subset(merged_multiome, cells = cells_keep)

obj_sub$group_for_DE <- ifelse(
  obj_sub$sampleID %in% c("PDX1_Control", "PDX3_Control"),
  "PDX",
  "Recurrent"
)

table(obj_sub$sampleID)
table(obj_sub$group_for_DE)

## --------------------------------------------------
## 2. Select top 20% lowESR1_nonprolif cells WITHIN each sample
## --------------------------------------------------

obj_sub$top20_lowESR1_nonprolif <- FALSE

for (s in unique(obj_sub$sampleID)) {
  idx <- which(obj_sub$sampleID == s)
  scores <- obj_sub$lowESR1_nonprolif_difference[idx]
  
  thr <- quantile(scores, probs = 0.80, na.rm = TRUE)
  obj_sub$top20_lowESR1_nonprolif[idx] <- scores >= thr
}

table(obj_sub$sampleID, obj_sub$top20_lowESR1_nonprolif)

obj_top20 <- subset(obj_sub, subset = top20_lowESR1_nonprolif)

table(obj_top20$sampleID)
table(obj_top20$group_for_DE)

## optional sanity plots
VlnPlot(
  obj_sub,
  features = "lowESR1_nonprolif_difference",
  group.by = "sampleID",
  pt.size = 0
) + RotatedAxis()

DimPlot(
  obj_top20,
  reduction = "umap_harmony_rna",
  group.by = "sampleID",
  label = TRUE
)

## --------------------------------------------------
## 3. Pseudobulk counts per sample
## --------------------------------------------------

DefaultAssay(obj_top20) <- "RNA"
obj_top20[["RNA"]] <- JoinLayers(obj_top20[["RNA"]])

rna_counts <- GetAssayData(obj_top20, assay = "RNA", layer = "counts")

samples <- unique(obj_top20$sampleID)

pb_counts <- sapply(samples, function(s) {
  cells_s <- colnames(obj_top20)[obj_top20$sampleID == s]
  Matrix::rowSums(rna_counts[, cells_s, drop = FALSE])
})

pb_counts <- as.matrix(pb_counts)
colnames(pb_counts) <- samples

meta_pb <- data.frame(
  sampleID = samples,
  condition = ifelse(samples %in% c("PDX1_Control", "PDX3_Control"), "PDX", "Recurrent"),
  row.names = samples,
  stringsAsFactors = FALSE
)

meta_pb$condition <- factor(meta_pb$condition, levels = c("PDX", "Recurrent"))
meta_pb
colSums(pb_counts)

## --------------------------------------------------
## 4. DESeq2 pseudobulk DE
## --------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(pb_counts),
  colData = meta_pb,
  design = ~ condition
)

## filter low-count genes
keep_genes <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep_genes, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Recurrent", "PDX"))
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

head(res_df, 50)

write.csv(
  res_df,
  file.path(multiome_dir, "DE_top20_lowESR1_nonprolif_Recurrent_vs_PDX_pseudobulk.csv"),
  row.names = FALSE
)

## --------------------------------------------------
## 5. Quick volcano
## --------------------------------------------------

res_df$signif <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                        "yes", "no")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point(size = 1, alpha = 0.8) +
  theme_minimal()

## --------------------------------------------------
## 6. Save selected object
## --------------------------------------------------

saveRDS(
  obj_top20,
  file.path(multiome_dir, "top20_lowESR1_nonprolif_cells_5samples.rds")
)


frag_paths <- c(
  PDX1_Control = file.path(fragment_dir, "PDX1_Control", "atac_fragments.tsv.gz"),
  PDX3_Control = file.path(fragment_dir, "PDX3_Control", "atac_fragments.tsv.gz"),
  Tumor_5572CL = file.path(fragment_dir, "Tumor_5572CL", "atac_fragments.tsv.gz"),
  Tumor_52BC3L = file.path(fragment_dir, "Tumor_52BC3L", "atac_fragments.tsv.gz"),
  Tumor_644D9L = file.path(fragment_dir, "Tumor_644D9L", "atac_fragments.tsv.gz"),
  Tumor_73607L = file.path(fragment_dir, "Tumor_73607L", "atac_fragments.tsv.gz")
)

debug_fragment_file <- function(frag_path, sample_name, n = 20) {
  out <- list(
    sample = sample_name,
    path = frag_path,
    frag_exists = file.exists(frag_path),
    tbi_exists = file.exists(paste0(frag_path, ".tbi"))
  )
  
  if (!out$frag_exists) {
    out$status <- "fragment_missing"
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  
  con <- gzfile(frag_path, open = "rt")
  txt <- tryCatch(readLines(con, n = n), error = function(e) e)
  close(con)
  
  if (inherits(txt, "error")) {
    out$status <- paste0("read_error: ", txt$message)
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  
  out$n_lines_read <- length(txt)
  out$first_lines <- paste(txt[1:min(3, length(txt))], collapse = " || ")
  
  if (length(txt) == 0) {
    out$status <- "empty_file"
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  
  split_counts <- sapply(strsplit(txt, "\t"), length)
  out$min_cols <- min(split_counts)
  out$max_cols <- max(split_counts)
  
  if (out$max_cols < 4) {
    out$status <- "fewer_than_4_columns"
  } else {
    out$status <- "looks_ok"
  }
  
  as.data.frame(out, stringsAsFactors = FALSE)
}


read_fragment_barcodes_safe <- function(frag_path, n = 50000) {
  con <- gzfile(frag_path, open = "rt")
  on.exit(close(con))
  
  txt <- readLines(con, n = n)
  
  ## remove comment/header lines from cellranger-arc fragments
  txt <- txt[!grepl("^#", txt)]
  
  if (length(txt) == 0) {
    stop("No non-comment lines found in fragment file")
  }
  
  split_txt <- strsplit(txt, "\t")
  split_len <- sapply(split_txt, length)
  
  if (max(split_len) < 4) {
    stop("Fragment file did not parse into >=4 tab-separated columns after removing comment lines")
  }
  
  split_txt <- split_txt[split_len >= 4]
  bc <- vapply(split_txt, function(x) x[4], character(1))
  unique(bc)
}

check_barcode_match_safe <- function(obj, frag_path, sample_name) {
  tryCatch({
    frag_barcodes <- read_fragment_barcodes_safe(frag_path, n = 100000)
    obj_barcodes  <- Cells(obj)
    
    exact_match_n <- sum(obj_barcodes %in% frag_barcodes)
    
    ## fallback 1: remove sample prefix up to first underscore
    obj_stripped1 <- sub("^[^_]+_", "", obj_barcodes)
    stripped1_match_n <- sum(obj_stripped1 %in% frag_barcodes)
    
    ## fallback 2: remove everything up to the LAST underscore
    obj_stripped2 <- sub("^.*_", "", obj_barcodes)
    stripped2_match_n <- sum(obj_stripped2 %in% frag_barcodes)
    
    data.frame(
      sample = sample_name,
      n_obj_cells = length(obj_barcodes),
      n_frag_barcodes_sampled = length(frag_barcodes),
      exact_match_n = exact_match_n,
      exact_match_frac = exact_match_n / length(obj_barcodes),
      stripped1_match_n = stripped1_match_n,
      stripped1_match_frac = stripped1_match_n / length(obj_barcodes),
      stripped2_match_n = stripped2_match_n,
      stripped2_match_frac = stripped2_match_n / length(obj_barcodes),
      status = "ok"
    )
  }, error = function(e) {
    data.frame(
      sample = sample_name,
      n_obj_cells = length(Cells(obj)),
      n_frag_barcodes_sampled = NA,
      exact_match_n = NA,
      exact_match_frac = NA,
      stripped1_match_n = NA,
      stripped1_match_frac = NA,
      stripped2_match_n = NA,
      stripped2_match_frac = NA,
      status = paste0("error: ", e$message)
    )
  })
}

frag_debug <- dplyr::bind_rows(lapply(names(frag_paths), function(nm) {
  debug_fragment_file(frag_paths[[nm]], nm)
}))


frag_debug


barcode_check <- dplyr::bind_rows(
  check_barcode_match_safe(PDX1_Control, frag_paths["PDX1_Control"], "PDX1_Control"),
  check_barcode_match_safe(PDX3_Control, frag_paths["PDX3_Control"], "PDX3_Control"),
  check_barcode_match_safe(Tumor_5572CL_epithelial, frag_paths["Tumor_5572CL"], "Tumor_5572CL"),
  check_barcode_match_safe(Tumor_52BC3L_epithelial, frag_paths["Tumor_52BC3L"], "Tumor_52BC3L"),
  check_barcode_match_safe(Tumor_644D9L_epithelial, frag_paths["Tumor_644D9L"], "Tumor_644D9L"),
  check_barcode_match_safe(Tumor_73607L_epithelial, frag_paths["Tumor_73607L"], "Tumor_73607L")
)

barcode_check


library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(Matrix)

obj_list <- list(
  PDX1_Control = PDX1_Control,
  PDX3_Control = PDX3_Control,
  Tumor_5572CL = Tumor_5572CL_epithelial,
  Tumor_52BC3L = Tumor_52BC3L_epithelial,
  Tumor_644D9L = Tumor_644D9L_epithelial,
  Tumor_73607L = Tumor_73607L_epithelial
)

## collect peak ranges
peak_gr_list <- lapply(obj_list, function(obj) {
  granges(obj[["ATAC"]])
})

## combine peaks
all_peaks <- do.call(c, unname(peak_gr_list))

## optional but recommended: keep standard chromosomes
all_peaks <- keepStandardChromosomes(all_peaks, pruning.mode = "coarse")


union_peaks <- GenomicRanges::reduce(all_peaks, ignore.strand = TRUE)

## sanity checks
union_peaks
length(union_peaks)
seqlevels(union_peaks)


read_fragment_barcodes_safe <- function(frag_path, n = 100000) {
  con <- gzfile(frag_path, open = "rt")
  on.exit(close(con))
  
  txt <- readLines(con, n = n)
  txt <- txt[!grepl("^#", txt)]
  
  if (length(txt) == 0) stop("No non-comment lines found in fragment file")
  
  split_txt <- strsplit(txt, "\t")
  split_len <- sapply(split_txt, length)
  split_txt <- split_txt[split_len >= 4]
  
  if (length(split_txt) == 0) stop("No valid fragment rows with >=4 columns found")
  
  unique(vapply(split_txt, function(x) x[4], character(1)))
}


requantify_atac_on_union <- function(obj, frag_path, sample_name, union_peaks) {
  obj_cells <- Cells(obj)
  frag_barcodes <- read_fragment_barcodes_safe(frag_path, n = 100000)
  
  stripped <- sub("^.*_", "", obj_cells)
  match_frac <- sum(stripped %in% frag_barcodes) / length(obj_cells)
  
  message(sample_name, " barcode match frac after stripping = ", round(match_frac, 4))
  
  if (match_frac < 0.5) {
    stop("Barcode mismatch for sample: ", sample_name)
  }
  
  cells_use <- stripped
  names(cells_use) <- obj_cells
  
  frag_obj <- CreateFragmentObject(
    path = frag_path,
    cells = unique(unname(cells_use)),
    validate.fragments = TRUE
  )
  
  mat <- FeatureMatrix(
    fragments = frag_obj,
    features = union_peaks,
    cells = unique(unname(cells_use))
  )
  
  ## restore Seurat cell names
  barcode_map <- data.frame(
    obj_cell = names(cells_use),
    frag_cell = unname(cells_use),
    stringsAsFactors = FALSE
  )
  
  idx <- match(colnames(mat), barcode_map$frag_cell)
  colnames(mat) <- barcode_map$obj_cell[idx]
  
  ## keep only original object cells and same order
  keep_cells <- obj_cells[obj_cells %in% colnames(mat)]
  mat <- mat[, keep_cells, drop = FALSE]
  
  new_atac <- CreateChromatinAssay(
    counts = mat,
    ranges = union_peaks,
    fragments = frag_obj
  )
  
  try({
    Annotation(new_atac) <- Annotation(obj[["ATAC"]])
  }, silent = TRUE)
  
  obj[["ATAC_union"]] <- new_atac
  obj
}


PDX1_Control_test <- requantify_atac_on_union(
  obj = PDX1_Control,
  frag_path = frag_paths["PDX1_Control"],
  sample_name = "PDX1_Control",
  union_peaks = union_peaks
)

dim(GetAssayData(PDX1_Control_test, assay = "ATAC_union", layer = "counts"))
sum(colnames(GetAssayData(PDX1_Control_test, assay = "ATAC_union", layer = "counts")) == Cells(PDX1_Control_test))


dim(GetAssayData(PDX1_Control_test, assay = "ATAC_union", layer = "counts"))

sum(
  colnames(GetAssayData(PDX1_Control_test, assay = "ATAC_union", layer = "counts")) ==
    Cells(PDX1_Control_test)[Cells(PDX1_Control_test) %in%
                               colnames(GetAssayData(PDX1_Control_test, assay = "ATAC_union", layer = "counts"))]
)

ncol(GetAssayData(PDX1_Control_test, assay = "ATAC_union", layer = "counts"))
length(Cells(PDX1_Control_test))


PDX3_Control_requant <- requantify_atac_on_union(
  obj = PDX3_Control,
  frag_path = frag_paths["PDX3_Control"],
  sample_name = "PDX3_Control",
  union_peaks = union_peaks
)

Tumor_52BC3L_requant <- requantify_atac_on_union(
  obj = Tumor_52BC3L_epithelial,
  frag_path = frag_paths["Tumor_52BC3L"],
  sample_name = "Tumor_52BC3L",
  union_peaks = union_peaks
)

Tumor_644D9L_requant <- requantify_atac_on_union(
  obj = Tumor_644D9L_epithelial,
  frag_path = frag_paths["Tumor_644D9L"],
  sample_name = "Tumor_644D9L",
  union_peaks = union_peaks
)

Tumor_73607L_requant <- requantify_atac_on_union(
  obj = Tumor_73607L_epithelial,
  frag_path = frag_paths["Tumor_73607L"],
  sample_name = "Tumor_73607L",
  union_peaks = union_peaks
)


saveRDS(PDX3_Control_requant, file.path(requant_dir, "PDX3_Control_requant.rds"))
saveRDS(PDX1_Control_test, file.path(requant_dir, "PDX1_Control_test.rds"))
saveRDS(Tumor_52BC3L_requant, file.path(requant_dir, "Tumor_52BC3L_requant.rds"))
saveRDS(Tumor_644D9L_requant, file.path(requant_dir, "Tumor_644D9L_requant.rds"))
saveRDS(Tumor_73607L_requant, file.path(requant_dir, "Tumor_73607L_requant.rds"))


requant_summary <- data.frame(
  sample = c("PDX1_Control", "PDX3_Control", "Tumor_52BC3L", "Tumor_644D9L", "Tumor_73607L"),
  original_cells = c(
    ncol(PDX1_Control),
    ncol(PDX3_Control),
    ncol(Tumor_52BC3L_epithelial),
    ncol(Tumor_644D9L_epithelial),
    ncol(Tumor_73607L_epithelial)
  ),
  requant_cells = c(
    ncol(GetAssayData(PDX1_Control_test, assay = "ATAC_union", layer = "counts")),
    ncol(GetAssayData(PDX3_Control_requant, assay = "ATAC_union", layer = "counts")),
    ncol(GetAssayData(Tumor_52BC3L_requant, assay = "ATAC_union", layer = "counts")),
    ncol(GetAssayData(Tumor_644D9L_requant, assay = "ATAC_union", layer = "counts")),
    ncol(GetAssayData(Tumor_73607L_requant, assay = "ATAC_union", layer = "counts"))
  )
)

requant_summary$retention_frac <- requant_summary$requant_cells / requant_summary$original_cells
requant_summary


PDX1_Control_requant <- PDX1_Control_test

obj_list_requant <- list(
  PDX1_Control = PDX1_Control_requant,
  PDX3_Control = PDX3_Control_requant,
  Tumor_52BC3L = Tumor_52BC3L_requant,
  Tumor_644D9L = Tumor_644D9L_requant,
  Tumor_73607L = Tumor_73607L_requant
)

for (nm in names(obj_list_requant)) {
  obj_list_requant[[nm]]$sampleID <- nm
  obj_list_requant[[nm]]$group_for_DE <- ifelse(
    nm %in% c("PDX1_Control", "PDX3_Control"),
    "PDX",
    "Recurrent"
  )
}


merged_union <- merge(
  x = obj_list_requant[[1]],
  y = obj_list_requant[2:length(obj_list_requant)],
  merge.data = TRUE
)


table(merged_union$sampleID)
names(merged_union@assays)


saveRDS(merged_union, file.path(requant_dir, "merged_union.rds"))

merged_union <- readRDS(file.path(requant_dir, "merged_union.rds"))


Assays(merged_union)


head(obj_top20)


library(future)

plan(sequential)  # simplest fix: disable parallel workers
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB, raise if needed

merged_union <- SCTransform(
  merged_union,
  assay = "RNA",
  new.assay.name = "SCT_recalc",
  verbose = TRUE
)


DefaultAssay(merged_union) <- 'SCT_recalc'

merged_union$lowESR1_nonprolif_difference <- NULL
merged_union$lowESR1_nonprolif_up <- NULL
merged_union$lowESR1_nonprolif_dw <- NULL
merged_union$lowESR1_prolif_difference <- NULL
merged_union$lowESR1_prolif_dw <- NULL
merged_union$lowESR1_prolif_up <- NULL


merged_union <- AddModuleScore_UCell(merged_union, features = list(
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
merged_union$lowESR1_nonprolif_difference <- merged_union[[up_np]] - merged_union[[down_np]]
merged_union$lowESR1_prolif_difference    <- merged_union[[up_pr]] - merged_union[[down_pr]]


DefaultAssay(merged_union) <- "RNA"

## Seurat v5: merge leaves multiple layers; join them first
merged_union[["RNA"]] <- JoinLayers(merged_union[["RNA"]])

## now create a unified normalized data layer
merged_union <- NormalizeData(merged_union, verbose = FALSE)

## variable features + scaling + PCA
merged_union <- FindVariableFeatures(
  merged_union,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

merged_union <- ScaleData(
  merged_union,
  features = VariableFeatures(merged_union),
  verbose = FALSE
)

merged_union <- RunPCA(
  merged_union,
  features = VariableFeatures(merged_union),
  npcs = 50,
  verbose = FALSE
)

merged_union <- RunUMAP(
  merged_union,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "umap_rna_unintegrated",
  verbose = FALSE
)

merged_union <- FindNeighbors(
  merged_union,
  dims = 1:30,
  reduction = "pca",
  verbose = FALSE
)

merged_union <- FindClusters(
  merged_union,
  resolution = 0.4,
  verbose = FALSE
)


merged_union <- RunUMAP(
  merged_union,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "umap_rna_unintegrated",
  verbose = FALSE
)

merged_union <- FindNeighbors(
  merged_union,
  dims = 1:30,
  reduction = "pca",
  verbose = FALSE
)

merged_union <- FindClusters(
  merged_union,
  resolution = 0.4,
  verbose = FALSE
)


## =========================================================
## 5. Diagnostics to decide whether integration is needed
## =========================================================

p1 <- DimPlot(
  merged_union,
  reduction = "umap_rna_unintegrated",
  group.by = "sampleID",
  label = TRUE,
  repel = TRUE
) + ggtitle("Unintegrated RNA UMAP by sample")

p2 <- DimPlot(
  merged_union,
  reduction = "umap_rna_unintegrated",
  group.by = "tumor_class",
  label = TRUE,
  repel = TRUE
) + ggtitle("Unintegrated RNA UMAP by tumor class")

p3 <- FeaturePlot(
  merged_union,
  reduction = "umap_rna_unintegrated",
  features = "ESR1",
  order = TRUE
) + ggtitle("ESR1 on unintegrated RNA UMAP")

print(p1 + p2)
print(p3)

## sample composition across clusters
cluster_sample_table <- table(
  cluster = Idents(merged_union),
  sample = merged_union$sampleID
)
print(cluster_sample_table)
print(prop.table(cluster_sample_table, margin = 1))

## ESR1 per sample
VlnPlot(
  merged_union,
  features = "ESR1",
  group.by = "sampleID",
  pt.size = 0
) + RotatedAxis()

## =====


DefaultAssay(merged_union) <- "RNA"
merged_union[["RNA"]] <- JoinLayers(merged_union[["RNA"]])

merged_union <- NormalizeData(merged_union, verbose = FALSE)
merged_union <- FindVariableFeatures(
  merged_union,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)
merged_union <- ScaleData(
  merged_union,
  features = VariableFeatures(merged_union),
  verbose = FALSE
)
merged_union <- RunPCA(
  merged_union,
  features = VariableFeatures(merged_union),
  npcs = 50,
  verbose = FALSE
)

library(harmony)

merged_union <- RunHarmony(
  object = merged_union,
  group.by.vars = "sampleID",
  reduction.use = "pca",
  dims.use = 1:30,
  project.dim = FALSE
)

merged_union <- RunUMAP(
  merged_union,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony_rna",
  verbose = FALSE
)


merged_union <- FindNeighbors(
  merged_union,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

merged_union <- FindClusters(
  merged_union,
  resolution = 0.4,
  verbose = FALSE
)


DimPlot(merged_union, reduction = "umap_harmony_rna", group.by = "sampleID", label = TRUE)
DimPlot(merged_union, reduction = "umap_harmony_rna", group.by = "tumor_class", label = TRUE)
FeaturePlot(merged_union, reduction = "umap_harmony_rna", features = "ESR1", order = TRUE)

tab_harmony <- table(Idents(merged_union), merged_union$sampleID)
prop.table(tab_harmony, 1)


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
})

## --------------------------------------------------
## 1. Keep samples of interest
## --------------------------------------------------

cells_keep <- colnames(merged_union)[
  merged_union$sampleID %in% c(
    "PDX1_Control",
    "PDX3_Control",
    "Tumor_52BC3L",
    "Tumor_644D9L",
    "Tumor_73607L"
  )
]

obj_sub <- subset(merged_union, cells = cells_keep)

obj_sub$group_for_DE <- ifelse(
  obj_sub$sampleID %in% c("PDX1_Control", "PDX3_Control"),
  "PDX",
  "Recurrent"
)

table(obj_sub$sampleID)
table(obj_sub$group_for_DE)

## --------------------------------------------------
## 2. Select top 20% lowESR1_nonprolif cells WITHIN each sample
## --------------------------------------------------

obj_sub$top20_lowESR1_nonprolif <- FALSE

for (s in unique(obj_sub$sampleID)) {
  idx <- which(obj_sub$sampleID == s)
  scores <- obj_sub$lowESR1_nonprolif_difference[idx]
  
  thr <- quantile(scores, probs = 0.80, na.rm = TRUE)
  obj_sub$top20_lowESR1_nonprolif[idx] <- scores >= thr
}

table(obj_sub$sampleID, obj_sub$top20_lowESR1_nonprolif)

obj_top20 <- subset(obj_sub, subset = top20_lowESR1_nonprolif)

table(obj_top20$sampleID)
table(obj_top20$group_for_DE)

## optional sanity plots
VlnPlot(
  obj_sub,
  features = "lowESR1_nonprolif_difference",
  group.by = "sampleID",
  pt.size = 0
) + RotatedAxis()


Reductions(obj_top20)

DimPlot(
  obj_top20,
  reduction = "umap_harmony_rna",
  group.by = "sampleID",
  label = TRUE
)


Idents(obj_top20) <- 'tumor_class'

Idents(obj_top20)

library(UNCover)

DefaultAssay(obj_top20) <- 'ds'

Reductions(obj_top20)


res_recurrent_vs_PDX_lowESR1 <- uncover_run(
  object = obj_top20,
  group_col = "tumor_class",
  group1 = c("Recurrent_tumor"),
  group2 = c("PDX_untreated"),
  rna_assay = "SCT_recalc",
  peak_assay = "ATAC",
  candidate_r_threshold = 0.30,
  metacell_reduction = "umap_harmony_rna",
  metacell_dims = 1:2,
  ml_use_metacells = TRUE,
  ml_cells_per_metacell = 50,
  ml_metacell_reduction = "umap_harmony_rna",
  ml_metacell_dims = 1:2,
  n_perm_fast = 100,
  n_random_pairs = 100,
  n_perm_random = 20,
  peak_gene_p_method = "zscore",
  peak_gene_fdr_cutoff = 0.05,
  network_fdr_cor_max = 0.05,
  tf_alpha = 0.05,
  tf_pct_threshold = 0.050,
  metacell_method = "fixed_size_umap",
  cells_per_metacell = 50,
  ml_mode = "slurm",
  submit = TRUE,
  collect_results = TRUE,
  wait_for_results = TRUE,
  conda_env = "myenv",
  conda_sh = conda_sh,
  python = python_bin,
  ml_output_dir = uncover_dir,
  sbatch_dir = uncover_dir,
  partition = "allnodes",
  cpus_per_task = 4,
  mem_per_cpu = "2G",
  time = "08:00:00",
  peaks_per_job = 10,
  tf_peak_I_min = 0.025,
  tf_peak_low = 0.3,
  tf_peak_high = 0.8,
  chip_catalog_bed = chip_catalog_bed,
  run_chip_intersect = TRUE,
  intersect_mode = "local",
  tf_peak_require_chip = TRUE,
  checkpoint_dir = file.path(uncover_dir, "checkpoints"),
  checkpoint_file = file.path(uncover_dir, "checkpoints", "UNCover_pipeline_state.rds"),
  resume_from_checkpoint = TRUE,
  save_final_rds = TRUE,
  final_rds_file = file.path(uncover_dir, "UNCover_result_final.rds"),
  run_tf_dual_profile_analysis = TRUE,
  tf_dual_profile_top_n = 20,
  tf_dual_profile_rank_by = "score",
  tf_dual_profile_min_peaks = 1,
  tf_dual_profile_min_triplets = 1,
  tf_dual_profile_use_weighted = TRUE,
  save_network_html = TRUE
)


"umap_harmony_rna" %in% Reductions(obj_top20)


extract_effector_subnetwork <- function(
    result,
    target_genes,
    min_I = 0,
    max_fdr_pg = 0.05,
    min_abs_cor_pg = 0,
    out_dir = NULL,
    prefix = NULL,
    make_html = TRUE,
    per_gene = FALSE,
    include_all_peak_to_gene_edges_for_selected_peaks = FALSE,
    ignore_missing_genes = TRUE,
    run_tf_dual_profile_analysis = TRUE,
    tf_dual_profile_top_n = 20,
    tf_dual_profile_rank_by = "score",
    tf_dual_profile_min_peaks = 1,
    tf_dual_profile_min_triplets = 1,
    tf_dual_profile_use_weighted = TRUE,
    save_tf_dual_profile_plot = TRUE,
    print_tf_dual_profile_plot = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(visNetwork)
    library(scales)
    library(purrr)
  })
  
  stopifnot(!is.null(result))
  stopifnot(all(c("edges_tf_peak", "edges_peak_gene", "nodes") %in% names(result)))
  
  if (missing(target_genes) || length(target_genes) == 0) {
    stop("Please provide at least one gene in `target_genes`.")
  }
  
  target_genes <- unique(as.character(target_genes))
  target_genes <- target_genes[!is.na(target_genes) & nzchar(target_genes)]
  
  if (length(target_genes) == 0) {
    stop("`target_genes` is empty after cleaning.")
  }
  
  edges_tf_peak   <- result$edges_tf_peak
  edges_peak_gene <- result$edges_peak_gene
  nodes           <- result$nodes
  triplets_full   <- if ("triplets" %in% names(result)) result$triplets else NULL
  
  num <- function(x) suppressWarnings(as.numeric(x))
  
  safe_pick_col <- function(df, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  make_color_fun <- function(x, palette, min_span = 1) {
    x2 <- x[is.finite(x)]
    if (!length(x2)) {
      dom <- c(-1, 1)
    } else {
      a <- stats::quantile(abs(x2), probs = 0.98, na.rm = TRUE, names = FALSE)
      a <- max(a, min_span / 2)
      dom <- c(-a, a)
    }
    function(v) {
      vv <- suppressWarnings(as.numeric(v))
      vv[is.finite(vv)] <- pmin(dom[2], pmax(dom[1], vv[is.finite(vv)]))
      pal <- scales::col_numeric(palette = palette, domain = dom, na.color = "#BDBDBD")
      pal(vv)
    }
  }
  
  make_color_fun_pos <- function(x, palette) {
    x2 <- x[is.finite(x)]
    xmax <- if (length(x2)) max(x2, na.rm = TRUE) else 1
    if (!is.finite(xmax) || xmax <= 0) xmax <- 1
    function(v) {
      vv <- suppressWarnings(as.numeric(v))
      vv[is.finite(vv)] <- pmin(xmax, pmax(0, vv[is.finite(vv)]))
      pal <- scales::col_numeric(palette = palette, domain = c(0, xmax), na.color = "#9ECAE1")
      pal(vv)
    }
  }
  
  build_tf_summary_from_triplets <- function(triplets_sub) {
    out <- triplets_sub %>%
      group_by(.data$tf) %>%
      summarise(
        n_genes = n_distinct(.data$gene),
        n_peaks = n_distinct(.data$peak),
        sum_contrib = sum(.data$peak_contrib, na.rm = TRUE),
        mean_abs_influence = mean(abs(.data$peak_contrib), na.rm = TRUE),
        sign_consistency = dplyr::if_else(
          sum(abs(.data$peak_contrib), na.rm = TRUE) > 0,
          abs(sum(.data$peak_contrib, na.rm = TRUE)) / sum(abs(.data$peak_contrib), na.rm = TRUE),
          0
        ),
        tf_pct_change = dplyr::first(.data$tf_pct_change),
        tf_p_adj = if ("tf_p_adj" %in% colnames(triplets_sub)) dplyr::first(.data$tf_p_adj) else NA_real_,
        hub = if ("hub" %in% colnames(triplets_sub)) dplyr::first(.data$hub) else NA_real_,
        pr = if ("pr" %in% colnames(triplets_sub)) dplyr::first(.data$pr) else NA_real_,
        score = mean(.data$I_weighted, na.rm = TRUE) * dplyr::n(),
        n_peaks_informative = n_distinct(.data$peak[is.finite(.data$peak_contrib) & abs(.data$peak_contrib) > 0]),
        n_peaks_total = n_distinct(.data$peak),
        n_triplets = dplyr::n(),
        n_triplets_informative = sum(is.finite(.data$peak_contrib) & abs(.data$peak_contrib) > 0, na.rm = TRUE),
        .groups = "drop"
      )
    
    out
  }
  
  build_one_subnetwork <- function(gene_set, label_for_files = "subnetwork") {
    
    gene_set <- unique(as.character(gene_set))
    
    epg_seed <- edges_peak_gene %>%
      transmute(
        peak   = .data$from,
        gene   = .data$to,
        cor_pg = num(.data$cor_pg),
        fdr_pg = num(.data$fdr_pg)
      ) %>%
      dplyr::filter(
        .data$gene %in% gene_set,
        is.finite(.data$cor_pg),
        abs(.data$cor_pg) >= min_abs_cor_pg
      )
    
    if (!is.na(max_fdr_pg)) {
      epg_seed <- epg_seed %>%
        dplyr::filter(is.finite(.data$fdr_pg), .data$fdr_pg <= max_fdr_pg)
    }
    
    found_genes <- sort(unique(epg_seed$gene))
    missing_genes <- setdiff(gene_set, found_genes)
    
    if (length(found_genes) == 0) {
      stop(
        "No Peak→Gene edges found for requested genes after filters. Requested genes: ",
        paste(gene_set, collapse = ", ")
      )
    }
    
    peaks_target <- unique(epg_seed$peak)
    
    etf <- edges_tf_peak %>%
      transmute(
        tf         = .data$from,
        peak       = .data$to,
        I          = if ("I" %in% colnames(edges_tf_peak)) num(.data$I) else NA_real_,
        I_weighted = if ("I_weighted" %in% colnames(edges_tf_peak)) num(.data$I_weighted) else NA_real_,
        chip_hits  = if ("chip_hits" %in% colnames(edges_tf_peak)) suppressWarnings(as.integer(.data$chip_hits)) else NA_integer_,
        tf_peak_cor = if ("tf_peak_cor" %in% colnames(edges_tf_peak)) num(.data$tf_peak_cor) else NA_real_
      ) %>%
      mutate(
        I_filter = dplyr::case_when(
          is.finite(.data$I_weighted) ~ .data$I_weighted,
          is.finite(.data$I) ~ .data$I,
          TRUE ~ NA_real_
        )
      ) %>%
      dplyr::filter(
        .data$peak %in% peaks_target,
        is.finite(.data$I_filter),
        .data$I_filter >= min_I
      )
    
    if (nrow(etf) == 0) {
      stop("No TF→Peak edges found for the selected peaks after I filtering.")
    }
    
    epg <- if (isTRUE(include_all_peak_to_gene_edges_for_selected_peaks)) {
      tmp <- edges_peak_gene %>%
        transmute(
          peak   = .data$from,
          gene   = .data$to,
          cor_pg = num(.data$cor_pg),
          fdr_pg = num(.data$fdr_pg)
        ) %>%
        dplyr::filter(
          .data$peak %in% peaks_target,
          is.finite(.data$cor_pg),
          abs(.data$cor_pg) >= min_abs_cor_pg
        )
      
      if (!is.na(max_fdr_pg)) {
        tmp <- tmp %>%
          dplyr::filter(is.finite(.data$fdr_pg), .data$fdr_pg <= max_fdr_pg)
      }
      tmp
    } else {
      epg_seed
    }
    
    if (nrow(epg) == 0) {
      stop("No Peak→Gene edges remain after expansion/filtering.")
    }
    
    triplets_sub <- NULL
    
    if (!is.null(triplets_full) && nrow(triplets_full) > 0) {
      if (all(c("tf", "peak", "gene") %in% colnames(triplets_full))) {
        triplets_sub <- triplets_full %>%
          dplyr::filter(
            .data$gene %in% unique(epg$gene),
            .data$peak %in% unique(etf$peak),
            .data$tf %in% unique(etf$tf)
          )
        
        if (!is.na(max_fdr_pg) && "fdr_pg" %in% colnames(triplets_sub)) {
          triplets_sub <- triplets_sub %>%
            dplyr::filter(is.finite(.data$fdr_pg), .data$fdr_pg <= max_fdr_pg)
        }
        
        if ("cor_pg" %in% colnames(triplets_sub)) {
          triplets_sub <- triplets_sub %>%
            dplyr::filter(is.finite(.data$cor_pg), abs(.data$cor_pg) >= min_abs_cor_pg)
        }
        
        if ("I_weighted" %in% colnames(triplets_sub)) {
          triplets_sub <- triplets_sub %>%
            dplyr::filter(is.finite(.data$I_weighted), .data$I_weighted >= min_I)
        } else if ("I" %in% colnames(triplets_sub)) {
          triplets_sub <- triplets_sub %>%
            dplyr::filter(is.finite(.data$I), .data$I >= min_I)
        }
      }
    }
    
    if (is.null(triplets_sub) || nrow(triplets_sub) == 0) {
      triplets_sub <- etf %>%
        transmute(
          tf,
          peak,
          I = .data$I,
          I_weighted = dplyr::case_when(
            is.finite(.data$I_weighted) ~ .data$I_weighted,
            is.finite(.data$I) ~ .data$I,
            TRUE ~ NA_real_
          ),
          chip_hits,
          tf_peak_cor
        ) %>%
        inner_join(epg, by = "peak", relationship = "many-to-many") %>%
        arrange(desc(.data$I_weighted), desc(abs(.data$cor_pg)), .data$tf, .data$peak, .data$gene)
    }
    
    if (nrow(triplets_sub) == 0) {
      stop("No TF→Peak→Gene triplets could be built.")
    }
    
    # force restriction to target genes for final subnetwork
    triplets_sub <- triplets_sub %>%
      dplyr::filter(.data$gene %in% unique(epg$gene))
    
    sub_names <- unique(c(triplets_sub$tf, triplets_sub$peak, triplets_sub$gene))
    
    nodes_sub <- nodes %>%
      dplyr::filter(.data$name %in% sub_names) %>%
      mutate(
        tf_pct_change = if ("tf_pct_change" %in% colnames(nodes)) num(.data$tf_pct_change) else NA_real_,
        peak_log2FC   = if ("peak_log2FC" %in% colnames(nodes)) num(.data$peak_log2FC) else NA_real_,
        gene_log2FC   = if ("gene_log2FC" %in% colnames(nodes)) num(.data$gene_log2FC) else NA_real_,
        tf_p_adj      = if ("tf_p_adj" %in% colnames(nodes)) num(.data$tf_p_adj) else NA_real_,
        hub           = if ("hub" %in% colnames(nodes)) num(.data$hub) else NA_real_,
        pr            = if ("pr" %in% colnames(nodes)) num(.data$pr) else NA_real_
      ) %>%
      distinct(.data$name, .keep_all = TRUE)
    
    tf_node_tbl <- nodes_sub %>%
      dplyr::filter(.data$type == "TF") %>%
      dplyr::transmute(
        tf = .data$name,
        tf_pct_change = .data$tf_pct_change,
        tf_p_adj = .data$tf_p_adj,
        hub = .data$hub,
        pr = .data$pr
      )
    
    peak_node_tbl <- nodes_sub %>%
      dplyr::filter(.data$type == "Peak") %>%
      dplyr::transmute(peak = .data$name, peak_log2FC = .data$peak_log2FC)
    
    gene_node_tbl <- nodes_sub %>%
      dplyr::filter(.data$type == "Gene") %>%
      dplyr::transmute(gene = .data$name, gene_log2FC = .data$gene_log2FC)
    
    if (!"tf_pct_change" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>% left_join(tf_node_tbl %>% select(tf, tf_pct_change), by = "tf")
    }
    if (!"tf_p_adj" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>% left_join(tf_node_tbl %>% select(tf, tf_p_adj), by = "tf")
    }
    if (!"hub" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>% left_join(tf_node_tbl %>% select(tf, hub), by = "tf")
    }
    if (!"pr" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>% left_join(tf_node_tbl %>% select(tf, pr), by = "tf")
    }
    if (!"peak_log2FC" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>% left_join(peak_node_tbl, by = "peak")
    }
    if (!"gene_log2FC" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>% left_join(gene_node_tbl, by = "gene")
    }
    
    if (!"I_weighted" %in% colnames(triplets_sub)) {
      if ("I" %in% colnames(triplets_sub)) {
        triplets_sub <- triplets_sub %>% mutate(I_weighted = num(.data$I))
      } else {
        triplets_sub <- triplets_sub %>% mutate(I_weighted = NA_real_)
      }
    }
    
    if (!"tf_peak_cor" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>% mutate(tf_peak_cor = NA_real_)
    }
    
    if (!"peak_contrib" %in% colnames(triplets_sub)) {
      triplets_sub <- triplets_sub %>%
        mutate(
          I_weighted  = ifelse(is.finite(.data$I_weighted), .data$I_weighted, 0),
          cor_pg      = ifelse(is.finite(.data$cor_pg), .data$cor_pg, 0),
          peak_log2FC = ifelse(is.finite(.data$peak_log2FC), .data$peak_log2FC, 0),
          peak_contrib = .data$I_weighted * .data$cor_pg * .data$peak_log2FC
        )
    }
    
    edges_tf_peak_sub <- triplets_sub %>%
      transmute(
        from      = .data$tf,
        to        = .data$peak,
        edge_type = "TF_to_Peak",
        I         = dplyr::case_when(
          "I" %in% colnames(triplets_sub) & is.finite(.data$I) ~ .data$I,
          TRUE ~ .data$I_weighted
        ),
        I_weighted = .data$I_weighted,
        chip_hits = if ("chip_hits" %in% colnames(triplets_sub)) .data$chip_hits else NA_integer_,
        cor_pg    = NA_real_,
        fdr_pg    = NA_real_
      ) %>%
      distinct()
    
    edges_peak_gene_sub <- triplets_sub %>%
      transmute(
        from      = .data$peak,
        to        = .data$gene,
        edge_type = "Peak_to_Gene",
        I         = NA_real_,
        I_weighted = NA_real_,
        chip_hits = NA_integer_,
        cor_pg    = .data$cor_pg,
        fdr_pg    = .data$fdr_pg
      ) %>%
      distinct()
    
    edges_sub <- bind_rows(edges_tf_peak_sub, edges_peak_gene_sub)
    
    gene_summary <- triplets_sub %>%
      group_by(.data$gene) %>%
      summarise(
        n_peaks = n_distinct(.data$peak),
        n_tfs = n_distinct(.data$tf),
        n_triplets = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(.data$n_triplets), desc(.data$n_peaks), desc(.data$n_tfs), .data$gene)
    
    tf_summary_sub <- build_tf_summary_from_triplets(triplets_sub)
    
    peak_summary <- triplets_sub %>%
      group_by(.data$peak) %>%
      summarise(
        n_tfs = n_distinct(.data$tf),
        n_genes = n_distinct(.data$gene),
        n_triplets = n(),
        max_abs_cor_pg = max(abs(.data$cor_pg), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(.data$n_triplets), desc(.data$max_abs_cor_pg), .data$peak)
    
    html_file <- NULL
    
    if (isTRUE(make_html)) {
      tf_pal <- make_color_fun(nodes_sub$tf_pct_change, c("#2166AC", "#FFFFFF", "#B2182B"), min_span = 0.1)
      pk_pal <- make_color_fun(nodes_sub$peak_log2FC, c("#313695", "#FFFFFF", "#A50026"), min_span = 0.5)
      gn_pal <- make_color_fun(nodes_sub$gene_log2FC, c("#2C7BB6", "#FFFFFF", "#D7191C"), min_span = 0.5)
      I_pal  <- make_color_fun_pos(edges_tf_peak_sub$I_weighted, c("#DEEBF7", "#3182BD"))
      r_pal  <- make_color_fun(edges_peak_gene_sub$cor_pg, c("#313695", "#FFFFFF", "#A50026"), min_span = 0.2)
      
      shape_map <- c(TF = "triangle", Peak = "diamond", Gene = "ellipse")
      target_gene_set <- unique(gene_set)
      
      nodes_v <- nodes_sub %>%
        mutate(
          id    = paste0(.data$type, ":", .data$name),
          label = .data$name,
          group = .data$type,
          shape = shape_map[.data$type],
          color.background = case_when(
            .data$type == "TF"   ~ tf_pal(.data$tf_pct_change),
            .data$type == "Peak" ~ pk_pal(.data$peak_log2FC),
            .data$type == "Gene" ~ gn_pal(.data$gene_log2FC),
            TRUE ~ "#BDBDBD"
          ),
          color.border = case_when(
            .data$type == "Gene" & .data$name %in% target_gene_set ~ "#000000",
            TRUE ~ "#333333"
          ),
          borderWidth = case_when(
            .data$type == "Gene" & .data$name %in% target_gene_set ~ 4,
            TRUE ~ 1
          ),
          title = case_when(
            .data$type == "TF" ~ paste0(
              .data$name,
              " | TF %Δ = ",
              ifelse(is.finite(.data$tf_pct_change), sprintf("%+.1f%%", 100 * .data$tf_pct_change), "NA")
            ),
            .data$type == "Peak" ~ paste0(
              .data$name,
              " | Peak log2FC = ",
              ifelse(is.finite(.data$peak_log2FC), sprintf("%+.2f", .data$peak_log2FC), "NA")
            ),
            .data$type == "Gene" ~ paste0(
              .data$name,
              " | Gene log2FC = ",
              ifelse(is.finite(.data$gene_log2FC), sprintf("%+.2f", .data$gene_log2FC), "NA"),
              ifelse(.data$name %in% target_gene_set, " | target gene", "")
            )
          )
        ) %>%
        transmute(
          id, label, group, shape,
          color.background, color.border, borderWidth, title
        )
      
      edges_v <- edges_sub %>%
        mutate(
          from_id = paste0(ifelse(.data$edge_type == "TF_to_Peak", "TF", "Peak"), ":", .data$from),
          to_id   = paste0(ifelse(.data$edge_type == "TF_to_Peak", "Peak", "Gene"), ":", .data$to),
          value   = case_when(
            .data$edge_type == "TF_to_Peak" ~ pmax(1, 12 * .data$I_weighted),
            TRUE ~ pmax(1, 12 * abs(.data$cor_pg))
          ),
          color   = case_when(
            .data$edge_type == "TF_to_Peak" ~ I_pal(.data$I_weighted),
            TRUE ~ r_pal(.data$cor_pg)
          ),
          title   = case_when(
            .data$edge_type == "TF_to_Peak" ~ paste0(
              .data$from, " → ", .data$to,
              " | I_weighted=", round(.data$I_weighted, 3),
              ifelse(is.finite(.data$chip_hits), paste0(" | ChIP hits=", .data$chip_hits), "")
            ),
            TRUE ~ paste0(
              .data$from, " → ", .data$to,
              " | r=", round(.data$cor_pg, 3),
              ifelse(is.finite(.data$fdr_pg), paste0(" | FDR=", signif(.data$fdr_pg, 3)), "")
            )
          ),
          arrows = "to"
        ) %>%
        transmute(from = .data$from_id, to = .data$to_id, value, color, title, arrows)
      
      g_vis <- visNetwork(nodes_v, edges_v, width = "100%", height = "850px") %>%
        visLegend(
          addNodes = data.frame(
            label = c("TF", "Peak", "Gene", "Target gene"),
            shape = c("triangle", "diamond", "ellipse", "ellipse"),
            color = c("#BDBDBD", "#BDBDBD", "#BDBDBD", "#BDBDBD"),
            borderWidth = c(1, 1, 1, 4),
            stringsAsFactors = FALSE
          ),
          useGroups = FALSE
        ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visPhysics(stabilization = TRUE)
      
      if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        html_file <- file.path(out_dir, paste0(label_for_files, ".html"))
        visSave(g_vis, html_file)
      }
    }
    
    tf_dual_profiles <- NULL
    tf_dual_profile_summary <- NULL
    tf_dual_profile_plot <- NULL
    tf_dual_profile_summary_file <- NULL
    tf_dual_profile_plot_file <- NULL
    
    if (isTRUE(run_tf_dual_profile_analysis) &&
        nrow(triplets_sub) > 0 &&
        nrow(tf_summary_sub) > 0 &&
        exists("summarize_tf_dual_profiles", mode = "function")) {
      
      req_cols <- c(
        "tf", "peak", "gene", "tf_peak_cor", "tf_pct_change",
        "peak_log2FC", "gene_log2FC", "I_weighted", "peak_contrib"
      )
      missing_req <- setdiff(req_cols, colnames(triplets_sub))
      
      if (length(missing_req) == 0) {
        message("[subnetwork] Running TF dual-profile post-analysis for ", label_for_files)
        
        tf_dual_profiles <- tryCatch(
          summarize_tf_dual_profiles(
            list(
              triplets = triplets_sub,
              tf_summary = tf_summary_sub
            )
          ),
          error = function(e) {
            warning("TF dual-profile analysis failed for ", label_for_files, ": ", conditionMessage(e))
            NULL
          }
        )
        
        if (!is.null(tf_dual_profiles) &&
            !is.null(tf_dual_profiles$tf_dual_profile_summary)) {
          tf_dual_profile_summary <- tf_dual_profiles$tf_dual_profile_summary
          
          if (!is.null(out_dir)) {
            tf_dual_profile_summary_file <- file.path(
              out_dir, paste0(label_for_files, "_tf_dual_profile_summary.csv")
            )
            readr::write_csv(tf_dual_profile_summary, tf_dual_profile_summary_file)
          }
          
          if (exists("plot_tf_dual_profiles", mode = "function")) {
            tf_dual_profile_plot <- tryCatch(
              plot_tf_dual_profiles(
                tf_dual_profile_summary,
                top_n = tf_dual_profile_top_n,
                rank_by = tf_dual_profile_rank_by,
                min_peaks = tf_dual_profile_min_peaks,
                min_triplets = tf_dual_profile_min_triplets,
                use_weighted = tf_dual_profile_use_weighted
              ),
              error = function(e) {
                warning("TF dual-profile plot failed for ", label_for_files, ": ", conditionMessage(e))
                NULL
              }
            )
          }
          
          if (!is.null(tf_dual_profile_plot) && isTRUE(print_tf_dual_profile_plot)) {
            print(tf_dual_profile_plot)
          }
          
          if (!is.null(tf_dual_profile_plot) &&
              !is.null(out_dir) &&
              isTRUE(save_tf_dual_profile_plot) &&
              requireNamespace("ggplot2", quietly = TRUE)) {
            tf_dual_profile_plot_file <- file.path(
              out_dir, paste0(label_for_files, "_tf_dual_profiles.pdf")
            )
            ggplot2::ggsave(
              filename = tf_dual_profile_plot_file,
              plot = tf_dual_profile_plot,
              width = 12,
              height = 8,
              units = "in"
            )
          }
        }
      } else {
        warning(
          "Skipping TF dual-profile analysis for ", label_for_files,
          " because triplets_sub is missing required columns: ",
          paste(missing_req, collapse = ", ")
        )
      }
    }
    
    if (!is.null(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      readr::write_csv(triplets_sub,        file.path(out_dir, paste0(label_for_files, "_triplets.csv")))
      readr::write_csv(edges_tf_peak_sub,   file.path(out_dir, paste0(label_for_files, "_edges_tf_peak.csv")))
      readr::write_csv(edges_peak_gene_sub, file.path(out_dir, paste0(label_for_files, "_edges_peak_gene.csv")))
      readr::write_csv(edges_sub,           file.path(out_dir, paste0(label_for_files, "_edges_all.csv")))
      readr::write_csv(nodes_sub,           file.path(out_dir, paste0(label_for_files, "_nodes.csv")))
      readr::write_csv(gene_summary,        file.path(out_dir, paste0(label_for_files, "_gene_summary.csv")))
      readr::write_csv(tf_summary_sub,      file.path(out_dir, paste0(label_for_files, "_tf_summary.csv")))
      readr::write_csv(peak_summary,        file.path(out_dir, paste0(label_for_files, "_peak_summary.csv")))
    }
    
    list(
      target_genes_requested = gene_set,
      target_genes_found = found_genes,
      target_genes_missing = missing_genes,
      peaks = sort(unique(triplets_sub$peak)),
      tfs = sort(unique(triplets_sub$tf)),
      genes = sort(unique(triplets_sub$gene)),
      triplets = triplets_sub,
      edges_tf_peak = edges_tf_peak_sub,
      edges_peak_gene = edges_peak_gene_sub,
      edges = edges_sub,
      nodes = nodes_sub,
      gene_summary = gene_summary,
      tf_summary = tf_summary_sub,
      peak_summary = peak_summary,
      tf_dual_profiles = tf_dual_profiles,
      tf_dual_profile_summary = tf_dual_profile_summary,
      tf_dual_profile_plot = tf_dual_profile_plot,
      tf_dual_profile_summary_file = tf_dual_profile_summary_file,
      tf_dual_profile_plot_file = tf_dual_profile_plot_file,
      html_file = html_file
    )
  }
  
  genes_present_anywhere <- unique(as.character(edges_peak_gene$to))
  genes_found_input <- intersect(target_genes, genes_present_anywhere)
  genes_missing_input <- setdiff(target_genes, genes_present_anywhere)
  
  if (length(genes_missing_input) > 0 && !isTRUE(ignore_missing_genes)) {
    stop(
      "These genes are not present in result$edges_peak_gene: ",
      paste(genes_missing_input, collapse = ", ")
    )
  }
  
  if (length(genes_found_input) == 0) {
    stop("None of the requested genes are present in result$edges_peak_gene.")
  }
  
  if (is.null(prefix)) {
    if (length(genes_found_input) == 1) {
      prefix <- paste0("subnetwork_", genes_found_input)
    } else {
      prefix <- "subnetwork_gene_set"
    }
  }
  
  combined <- build_one_subnetwork(
    gene_set = genes_found_input,
    label_for_files = prefix
  )
  
  per_gene_results <- NULL
  if (isTRUE(per_gene)) {
    per_gene_results <- purrr::map(
      genes_found_input,
      function(g) {
        gene_prefix <- paste0(prefix, "_", g)
        tryCatch(
          build_one_subnetwork(
            gene_set = g,
            label_for_files = gene_prefix
          ),
          error = function(e) {
            list(
              target_genes_requested = g,
              error = conditionMessage(e)
            )
          }
        )
      }
    )
    names(per_gene_results) <- genes_found_input
  }
  
  list(
    requested_genes = target_genes,
    found_genes = genes_found_input,
    missing_genes = genes_missing_input,
    combined = combined,
    per_gene = per_gene_results
  )
}


PIM3_kinase <- c(
  "AP1B1", "ACOD1", "EML3", "PML", "CAPN15", "CD4", "NINJ1", "TNFSF9", "ZNF710", "ATF3", "AP5Z1", "STXBP2", "NOD2", "SZRD1", "STX11", "FZR1", "NIPAL4",
  "STAT6", "SLAMF7", "FBRS", "GBP2", "N4BP1", "KLF13", "PRKCD", "GADD45B", "PLK3", "HLA-A", "VPS37C", "HLA-B", "KLF16", "EHD1", "VASP", "FOSL2", "PARP12", 
  "PARP10", "GRB2", "FERMT3", "HLA-DRB1", "IL1RN", "SPI1", "PPP1R18", "HDAC10", "PDGFB", "ELK1", "IFI30", "ICAM1", "TRPM2", "IL4I1", "HELZ2", "ZC3H12A", 
  "RASSF5", "FTH1", "FFAR3", "KCNN4", "LGALS9", "JUNB", "ZBTB7B", "IER3", "STAT5A", "DUSP5", "LYN", "CCR1", "UPF1", "IL10", "DUSP2", "SPHK1", "IL10RA",
  "LILRB1", "LILRB2", "AMPD3", "LILRB3", "RHOG", "LILRB4", "TYK2", "GMIP", "VAV1", "TLR2", "IL1B", "SIGLEC1", "MMP19", "PLEKHM2", "ARHGEF2", "RAPGEF1",
  "SIGLEC9", "CH507-9B2.1", "STK40", "ZNF598", "DOT1L", "CSRNP1", "RAB5C", "GRAMD1A", "CD40", "ELL", "FLII", "CTSZ", "AKAP17A", "LITAF", "SOCS3", "SBNO2", 
  "SOCS1", "MED15", "CTSL", "FLVCR2", "TFE3", "SNX8", "CCR5", "CTSD", "SERPINB2", "RFTN1", "GAA", "CREG1", "LAPTM5", "UBALD2", "SCAF1", "MYO9B", "HCK",
  "XAB2", "GAK", "SPSB1", "PLEKHO1", "KCNQ1", "SHKBP1", "NUCB1", "PKN1", "PLEKHO2", "PNPLA6", "TRIP10", "C1QA", "GRN", "SLC43A3", "CTDP1", "IL27", "OGFR",
  "SLC43A2", "GNAI2", "PTGS2", "IFIH1", "FCGRT", "GPR132", "SDF4", "MGAT1", "PTPN1", "TGFB1", "FBXW5", "PPAN-P2RY11", "TAP2", "TAP1", "TYROBP", "ARHGAP27",
  "TICAM1", "AP5B1", "BATF", "MOB3C", "OSCAR", "TMEM259", "MAPKAPK2", "RAB35", "EHBP1L1", "PPIF", "CHMP4B", "ZNF316", "ADAM8", "CD68", "PTPN6", "GRINA",
  "TYMP", "RELA", "CD83", "CRTC2", "NCF1C", "TNF", "METRNL", "RELB", "CXCL16", "RAB20", "TMEM127", "LAMP1", "SH3BP1", "RNF19B", "CCRL2", "RELT", "ZNF787",
  "SH2B2", "SLC15A3", "MAP2K3", "ARID5A", "TRAF1", "TRAF3", "CYBA", "NFKB1", "MESDC1", "RGCC", "GNB2", "MAFG", "MAP3K11", "CLASRP", "MAFK", "FTL", "PFKFB3",
  "CD274", "C5AR1", "INSIG1", "GPR84", "SLC2A6", "C17ORF96", "ETS2", "ARHGAP4", "HCAR2", "FPGS", "SIRPA", "TOM1", "EFHD2", "TNFRSF4", "MSC", "REXO1", 
  "ATP13A1", "PILRA", "AXIN1", "ARAP1", "TNIP1", "DNM2", "TNFRSF1B", "CLCN7", "TNIP2", "IRF5", "TNIP3", "ADORA2A", "BCL3", "MAN2B1", "IRF7", "LCP2",
  "PPP1R12C", "NFKBIE", "CXCL8", "RP1-47A17.1", "KDM5C", "CCL3L3", "BCL2A1", "TRADD", "PIK3R5", "TNFAIP3", "DTX2", "CXCL2", "TNFAIP2", "CXCL3", "GRK2", 
  "APH1A", "C15ORF39", "DPP7", "DRAM1", "MFSD7", "RHOT2", "IRAK2", "PSAP", "PPP6R1", "IRAK1", "EMILIN2", "ZNF865", "TNFRSF14", "MCL1", "KDM6B", "IL15RA",
  "KHNYN", "ANXA11", "PLAUR", "USF2", "SOD2", "PFKL", "CDC37", "FMNL1", "SH3D21", "PIK3AP1", "MYO1G", "ANAPC2", "B4GALT5", "NRROS", "VAC14", "ITGB2",
  "CEBPB", "TRABD", "SRC", "CEBPD", "UNC93B1", "PLEK", "FURIN", "USP12", "TCIRG1", "SQRDL", "SART1", "ARHGDIA", "GNA15", "HLX", "RPS6KA1", "ABHD16B",
  "TRIM25", "ITGAX", "CSK", "ABCD1", "CYTH1"
)


PIM3_kinase_net <- extract_effector_subnetwork(
  result = res_recurrent_vs_PDX_lowESR1,
  target_genes = PIM3_kinase,
  min_I = 0.05,
  max_fdr_pg = 0.05,
  min_abs_cor_pg = 0,
  out_dir = file.path(uncover_dir, "final_outputs", "subnetwork_pim3"),
  prefix = "PIM3_kinase",
  make_html = TRUE,
  per_gene = FALSE,
  run_tf_dual_profile_analysis = TRUE,
  tf_dual_profile_top_n = 20,
  tf_dual_profile_rank_by = "score",
  tf_dual_profile_min_peaks = 1,
  tf_dual_profile_min_triplets = 1,
  tf_dual_profile_use_weighted = TRUE
)

