##### Preparation of merged multiome Seurat object


# 1. Define input files
# 2. Load libraries
# 3. Create per-sample multiome objects
# 4. Merge samples
# 5. Compute QC metrics
# 6. Filter cells
# 7. Call peaks
# 8. Requantify counts on unified peak set
# 9. RNA preprocessing
# 10. ATAC preprocessing
# 11. Multi-modal integration
# 12. Save outputs


suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
})

###############################################################################
# 1. User-defined inputs
###############################################################################

# Named list of sample input directories and fragment files
sample_info <- data.frame(
  sample_id = c("Sample1", "Sample2", "Sample3", "Sample4"),
  count_dir = c(
    "/path/to/sample1/outs/filtered_feature_bc_matrix/",
    "/path/to/sample2/outs/filtered_feature_bc_matrix/",
    "/path/to/sample3/outs/filtered_feature_bc_matrix/",
    "/path/to/sample4/outs/filtered_feature_bc_matrix/"
  ),
  fragment_file = c(
    "/path/to/sample1/outs/atac_fragments.tsv.gz",
    "/path/to/sample2/outs/atac_fragments.tsv.gz",
    "/path/to/sample3/outs/atac_fragments.tsv.gz",
    "/path/to/sample4/outs/atac_fragments.tsv.gz"
  ),
  stringsAsFactors = FALSE
)

# MACS2 executable
macs2.path <- "/path/to/macs2"

# Output directory
outdir <- "/path/to/output_directory"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. Helper function to create one multiome Seurat object
###############################################################################

create_multiome_object <- function(count_dir, frag_path, sample_name, annotation) {
  counts <- Read10X(count_dir)

  if (!("Gene Expression" %in% names(counts))) {
    stop(sample_name, ": 'Gene Expression' matrix not found in Read10X output.")
  }
  if (!("Peaks" %in% names(counts))) {
    stop(sample_name, ": 'Peaks' matrix not found in Read10X output.")
  }

  obj <- CreateSeuratObject(
    counts = counts[["Gene Expression"]],
    assay = "RNA",
    project = sample_name
  )

  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts[["Peaks"]],
    sep = c(":", "-"),
    fragments = frag_path,
    annotation = annotation
  )

  obj$sampleID <- sample_name
  return(obj)
}

###############################################################################
# 3. Gene annotations
###############################################################################

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
genome(annotation) <- "hg38"

###############################################################################
# 4. Load all samples
###############################################################################

obj_list <- vector("list", length = nrow(sample_info))
names(obj_list) <- sample_info$sample_id

for (i in seq_len(nrow(sample_info))) {
  message("Loading ", sample_info$sample_id[i], " ...")
  obj_list[[i]] <- create_multiome_object(
    count_dir = sample_info$count_dir[i],
    frag_path = sample_info$fragment_file[i],
    sample_name = sample_info$sample_id[i],
    annotation = annotation
  )
}

###############################################################################
# 5. Merge all samples
###############################################################################

combined_obj <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = sample_info$sample_id,
  project = "MultiomeProject"
)

rm(obj_list)
gc()

###############################################################################
# 6. Compute ATAC QC metrics on merged object
###############################################################################

DefaultAssay(combined_obj) <- "ATAC"

combined_obj <- NucleosomeSignal(combined_obj)
combined_obj <- TSSEnrichment(combined_obj)

###############################################################################
# 7. Visual QC plots
###############################################################################

pdf(file.path(outdir, "QC_violin_plots.pdf"), width = 12, height = 8)
print(
  VlnPlot(
    object = combined_obj,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0
  )
)
dev.off()

###############################################################################
# 8. Cell filtering
###############################################################################
# Adjust thresholds as needed for the dataset

combined_obj <- subset(
  x = combined_obj,
  subset =
    nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

###############################################################################
# 9. Call peaks on the filtered merged object using MACS2
###############################################################################

DefaultAssay(combined_obj) <- "ATAC"

called_peaks <- CallPeaks(
  object = combined_obj,
  group.by = "sampleID",
  macs2.path = macs2.path
)

# Keep standard chromosomes only
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
called_peaks <- keepStandardChromosomes(called_peaks, pruning.mode = "coarse")
called_peaks <- called_peaks[as.character(seqnames(called_peaks)) %in% standard_chromosomes]

###############################################################################
# 10. Requantify ATAC counts on the MACS2-derived peak set
###############################################################################

macs2_counts <- FeatureMatrix(
  fragments = Fragments(combined_obj),
  features = called_peaks,
  cells = colnames(combined_obj)
)

stopifnot(identical(colnames(macs2_counts), colnames(combined_obj)))
stopifnot(nrow(macs2_counts) == length(called_peaks))

combined_obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(combined_obj),
  annotation = annotation
)

###############################################################################
# 11. RNA preprocessing
###############################################################################

DefaultAssay(combined_obj) <- "RNA"

combined_obj <- SCTransform(combined_obj, verbose = FALSE)
combined_obj <- RunPCA(combined_obj, assay = "SCT", verbose = FALSE)

###############################################################################
# 12. ATAC preprocessing using MACS2-derived peaks
###############################################################################

DefaultAssay(combined_obj) <- "peaks"

combined_obj <- FindTopFeatures(combined_obj, min.cutoff = 5)
combined_obj <- RunTFIDF(combined_obj)
combined_obj <- RunSVD(combined_obj, reduction.name = "lsi_peaks")

###############################################################################
# 13. Multiome integration with WNN
###############################################################################

combined_obj <- FindMultiModalNeighbors(
  combined_obj,
  reduction.list = list("pca", "lsi_peaks"),
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight"
)

combined_obj <- RunUMAP(
  combined_obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

combined_obj <- FindClusters(
  combined_obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.5
)

###############################################################################
# 14. Save outputs
###############################################################################

saveRDS(combined_obj, file = file.path(outdir, "combined_multiome_processed.rds"))
saveRDS(called_peaks, file = file.path(outdir, "called_peaks_granges.rds"))

pdf(file.path(outdir, "UMAP_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(combined_obj, reduction = "wnn.umap", group.by = "sampleID"))
dev.off()

pdf(file.path(outdir, "UMAP_by_cluster.pdf"), width = 8, height = 6)
print(DimPlot(combined_obj, reduction = "wnn.umap", label = TRUE))
dev.off()