## Generate ER+ cancer cell Atlas


# ============================================================
# Merge ER+ cancer-cell Seurat objects (Seurat v5-safe)
# - Keep RNA raw counts + SCT
# - Merge 100+ samples safely
# - JoinLayers where supported (RNA typically; SCT often not)
# - RNA PCA -> Harmony -> UMAP/Clustering
# - ESR1 low/high labels from RNA
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(Matrix)
})

# ----------------------------
# User settings
# ----------------------------
samples_csv <- "/path/to/sample_list_cancer_cells_102_samples.csv"
out_rds     <- "cancer_cells_merged_harmony_ESR1.rds"
out_failcsv <- "merge_failed_samples.csv"

batch_var   <- "sample_id"

# ESR1 labeling
esr1_gene      <- "ESR1"
low_quantile   <- 0.10
high_quantile  <- 0.90

# Dimensionality / clustering
set.seed(1)
npcs       <- 50
use_dims   <- 1:30
resolution <- 0.4

# Required assays in every object
required_assays <- c("RNA", "SCT")

# Seurat v5 layers to retain (when supported)
keep_layers <- list(
  RNA = c("counts"),
  SCT = c("data")
)

# ----------------------------
# Helpers
# ----------------------------
log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")

read_seurat_any <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("rds", "rda", "rdata")) {
    obj <- readRDS(path)
  } else if (ext == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Need package 'qs' to read: ", path)
    obj <- qs::qread(path)
  } else {
    stop("Unsupported file extension: ", ext, " for file: ", path)
  }
  
  if (!inherits(obj, "Seurat")) stop("Loaded object is not a Seurat object: ", path)
  obj
}

detect_col <- function(df, patterns, required = TRUE, what = "column") {
  hits <- names(df)[Reduce(`|`, lapply(patterns, function(p) grepl(p, names(df), ignore.case = TRUE)))]
  if (length(hits) == 0) {
    if (required) stop("Could not find ", what, ". Tried patterns: ", paste(patterns, collapse = ", "))
    return(NA_character_)
  }
  hits[[1]]
}

assay_supports_layers <- function(assay_obj) inherits(assay_obj, "Assay5")

standardize_one <- function(obj, row, sample_col, path_col, required_assays, keep_layers) {
  sample_id   <- as.character(row[[sample_col]])
  source_path <- as.character(row[[path_col]])
  
  if (!all(required_assays %in% names(obj@assays))) {
    stop("Missing required assays in ", sample_id,
         ". Required: ", paste(required_assays, collapse = ", "),
         " | Has: ", paste(names(obj@assays), collapse = ", "))
  }
  
  obj <- RenameCells(obj, add.cell.id = sample_id)
  
  md <- as.list(row)
  md$sample_id   <- sample_id
  md$source_path <- source_path
  obj <- AddMetaData(obj, metadata = md)
  
  DefaultAssay(obj) <- "RNA"
  
  obj <- DietSeurat(
    object    = obj,
    assays    = required_assays,
    layers    = keep_layers,
    dimreducs = NULL,
    graphs    = NULL,
    misc      = FALSE
  )
  
  obj
}

join_layers_if_supported <- function(seu, assay_name) {
  if (!assay_name %in% names(seu@assays)) return(seu)
  a <- seu[[assay_name]]
  if (!assay_supports_layers(a)) {
    log_msg("JoinLayers skipped for", assay_name, "| class:", paste(class(a), collapse = ", "))
    return(seu)
  }
  log_msg("Joining layers for", assay_name, "| class:", paste(class(a), collapse = ", "))
  JoinLayers(seu, assay = assay_name)
}

# -------- Harmony runner (version-robust, avoids argument partial-matching) --------
run_harmony_seurat <- function(seu, group.by.vars, reduction_name = "pca", dims_use = NULL,
                               assay_use = NULL, harmony_name = "harmony", verbose = FALSE) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Package 'harmony' is required.")
  }
  
  # Get the Seurat S3 method from the harmony package
  fn <- getS3method("RunHarmony", "Seurat", optional = TRUE)
  if (is.null(fn)) {
    # Sometimes it is registered from another namespace
    fn <- getS3method("RunHarmony", "Seurat")
  }
  fml <- names(formals(fn))
  
  args <- list(object = seu, group.by.vars = group.by.vars)
  
  # Different harmony versions use different parameter names; set only exact matches.
  if ("reduction.use" %in% fml) {
    args$reduction.use <- reduction_name
  } else if ("reduction" %in% fml) {
    args$reduction <- reduction_name
  }
  
  if (!is.null(dims_use)) {
    if ("dims.use" %in% fml) {
      args$dims.use <- dims_use
    } else if ("dims" %in% fml) {
      args$dims <- dims_use
    }
  }
  
  if (!is.null(assay_use) && "assay.use" %in% fml) {
    args$assay.use <- assay_use
  }
  
  if ("reduction.save" %in% fml) {
    args$reduction.save <- harmony_name
  } else if ("reduction.name" %in% fml) {
    args$reduction.name <- harmony_name
  }
  
  if ("verbose" %in% fml) args$verbose <- verbose
  
  do.call(fn, args)
}

# ----------------------------
# Preconditions
# ----------------------------
if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("Package 'harmony' is required. Install it in this R environment.")
}

# ----------------------------
# Read sample list
# ----------------------------
log_msg("Reading sample table:", samples_csv)
samples <- read_csv(samples_csv, show_col_types = FALSE)

path_col   <- detect_col(samples, patterns = c("^path$", "seurat", "rds", "qs", "file"), what = "path column")
sample_col <- detect_col(samples, patterns = c("^sample$", "sample_id", "^id$"), what = "sample id column")

log_msg("Using columns:", "sample_col =", sample_col, "| path_col =", path_col)
log_msg("Samples in table:", nrow(samples))

# ----------------------------
# Load + standardize objects
# ----------------------------
objs <- vector("list", nrow(samples))
names(objs) <- as.character(samples[[sample_col]])

failures <- list()

for (i in seq_len(nrow(samples))) {
  sid   <- as.character(samples[[sample_col]][i])
  spath <- as.character(samples[[path_col]][i])
  
  log_msg("Loading", sid, "from", spath)
  
  res <- tryCatch({
    obj <- read_seurat_any(spath)
    standardize_one(
      obj = obj,
      row = samples[i, ],
      sample_col = sample_col,
      path_col = path_col,
      required_assays = required_assays,
      keep_layers = keep_layers
    )
  }, error = function(e) {
    failures[[length(failures) + 1]] <<- tibble(
      sample_id = sid,
      path      = spath,
      error     = conditionMessage(e)
    )
    NULL
  })
  
  objs[[i]] <- res
}

ok_idx  <- which(vapply(objs, function(x) !is.null(x), logical(1)))
objs_ok <- objs[ok_idx]

log_msg("Loaded OK:", length(objs_ok), "/", nrow(samples))

if (length(failures) > 0) {
  fail_tbl <- bind_rows(failures)
  write_csv(fail_tbl, out_failcsv)
  log_msg("Wrote failures table:", out_failcsv)
}

if (length(objs_ok) < 2) stop("Need at least 2 valid Seurat objects to merge. Loaded OK: ", length(objs_ok))

# ----------------------------
# Merge
# ----------------------------
log_msg("Merging objects...")
merged <- merge(x = objs_ok[[1]], y = objs_ok[-1], merge.data = FALSE)

rm(objs, objs_ok)
invisible(gc())

log_msg("Merged cells:", ncol(merged))
log_msg("Assays present:", paste(names(merged@assays), collapse = ", "))
log_msg("Assay classes:", paste(
  names(merged@assays),
  vapply(merged@assays, function(a) paste(class(a), collapse = "/"), character(1)),
  sep = "=", collapse = " | "
))

if (!batch_var %in% colnames(merged@meta.data)) {
  stop("Batch variable '", batch_var, "' not present in merged metadata.")
}

# ----------------------------
# JoinLayers where supported
# ----------------------------
merged <- join_layers_if_supported(merged, "RNA")
merged <- join_layers_if_supported(merged, "SCT")

# ----------------------------
# RNA normalization + HVG + PCA
# ----------------------------
log_msg("Normalizing RNA (LogNormalize)...")
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

log_msg("Finding variable features (RNA)...")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

log_msg("Scaling RNA...")
merged <- ScaleData(merged, features = VariableFeatures(merged), verbose = FALSE)

log_msg("Running PCA (RNA)...")
merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = npcs, verbose = FALSE)

# ----------------------------
# Harmony batch correction (robust call)
# ----------------------------
log_msg("Running Harmony on:", batch_var)

merged <- RunHarmony(
  object        = merged,
  group.by.vars = batch_var,
  reduction.use = "pca",
  dims.use      = use_dims,
  assay.use     = "RNA",
  reduction.save = "harmony",
  verbose       = FALSE
)

# ----------------------------
# UMAP / neighbors / clusters (Harmony space)
# ----------------------------
log_msg("Running UMAP + neighbors + clusters (Harmony)...")
merged <- RunUMAP(merged, reduction = "harmony", dims = use_dims, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = use_dims, verbose = FALSE)
merged <- FindClusters(merged, resolution = resolution, verbose = FALSE)

# ----------------------------
# ESR1 group labeling from RNA
# ----------------------------
log_msg("Assigning ESR1 groups from RNA...")
DefaultAssay(merged) <- "RNA"

if (!esr1_gene %in% rownames(merged[["RNA"]])) {
  stop("Gene '", esr1_gene, "' not found in RNA assay. Check rownames for gene symbols.")
}

get_rna_mat <- function(seu, layer_name, slot_name) {
  a <- seu[["RNA"]]
  if (assay_supports_layers(a)) {
    GetAssayData(seu, assay = "RNA", layer = layer_name)
  } else {
    GetAssayData(seu, assay = "RNA", slot = slot_name)
  }
}

rna_counts <- get_rna_mat(merged, layer_name = "counts", slot_name = "counts")
rna_data   <- get_rna_mat(merged, layer_name = "data",   slot_name = "data")

esr1_counts <- as.numeric(rna_counts[esr1_gene, ])
esr1_logn   <- as.numeric(rna_data[esr1_gene, ])

merged$ESR1_counts <- esr1_counts
merged$ESR1_log1p  <- esr1_logn

nonzero <- which(esr1_counts > 0)
q_low  <- if (length(nonzero) > 50) as.numeric(quantile(esr1_logn[nonzero], probs = low_quantile,  na.rm = TRUE)) else NA_real_
q_high <- if (length(nonzero) > 50) as.numeric(quantile(esr1_logn[nonzero], probs = high_quantile, na.rm = TRUE)) else NA_real_

grp4 <- rep("ESR1_mid", length(esr1_counts))
grp4[esr1_counts == 0] <- "ESR1_zero"
if (!is.na(q_low))  grp4[esr1_counts > 0 & esr1_logn <= q_low]  <- "ESR1_low"
if (!is.na(q_high)) grp4[esr1_counts > 0 & esr1_logn >= q_high] <- "ESR1_high"

merged$ESR1_group_4 <- factor(grp4, levels = c("ESR1_zero", "ESR1_low", "ESR1_mid", "ESR1_high"))
merged$ESR1_low_vs_rest <- factor(ifelse(merged$ESR1_group_4 %in% c("ESR1_zero", "ESR1_low"), "LOW", "REST"),
                                  levels = c("LOW", "REST"))

log_msg("ESR1_group_4 counts:")
print(table(merged$ESR1_group_4, useNA = "ifany"))


##### Recompute SCT


DefaultAssay(merged) <- "RNA"

merged <- SCTransform(
  merged,
  assay = "RNA",
  new.assay.name = "SCT2",
  vst.flavor = "v2",
  return.only.var.genes = FALSE,   # <<< critical
  conserve.memory = TRUE,
  verbose = TRUE
)

# Verify ESR1 is now non-zero in SCT2 data
summary(as.numeric(GetAssayData(merged, assay="SCT2", layer="data")["ESR1", ]))










# ----------------------------
# Save
# ----------------------------
log_msg("Saving merged object:", out_rds)
saveRDS(merged, out_rds)

log_msg("DONE.")