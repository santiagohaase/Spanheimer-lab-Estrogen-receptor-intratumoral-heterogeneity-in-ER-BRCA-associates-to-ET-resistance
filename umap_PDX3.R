##### recaulculate PDX3-control sample

PDX3_Control <- subset(pdx, subset = sampleID == "PDX3_Control")

 

head(PDX3_Control)

# Get the metadata column names
meta_cols <- colnames(PDX3_Control@meta.data)

# Identify columns where "kinase" appears twice
duplicate_cols <- grep(".*kinase.*kinase.*", meta_cols, value = TRUE)

# Remove those columns
PDX3_Control@meta.data <- PDX3_Control@meta.data[, !colnames(PDX3_Control@meta.data) %in% duplicate_cols]

# Verify remaining columns
colnames(PDX3_Control@meta.data)

PDX3_Control@meta.data$ESR1_human_tf_ARCHS4_coexpression_UCell <- NULL
PDX3_Control@meta.data$GATA3_human_tf_ARCHS4_coexpression_UCell <- NULL
PDX3_Control@meta.data$RET_human_kinase_ARCHS4_coexpression_UCell <- NULL
#####################################

library(BiocParallel)

bp <- MulticoreParam(workers = 32)

# Quick test
tmp <- bplapply(1:4, function(x) Sys.getpid(), BPPARAM = bp)
unique(unlist(tmp))
# If you see multiple PIDs, it's actually parallelized



# Now score all missing signatures in parallel, in a single call
PDX3_Control <- AddModuleScore_UCell(
  obj      = PDX3_Control,
  features = signatures,    # named list of all missing signatures
  assay    = "SCT",
  slot     = "counts",
  BPPARAM  = bp,              # pass the parallel backend
  name     = NULL         # each signature gets "<SignatureName>_UCell"
)

bpstop(bp)

###############################

head(PDX3_Control)



saveRDS(PDX3_Control, "/path/to/PDX_analysis_R/PDX3_Control.RDS")

DefaultAssay(PDX3_Control) <- "SCT"

PDX3_Control <- RunPCA(PDX3_Control, features = VariableFeatures(object = PDX3_Control))

ElbowPlot(PDX3_Control)


PDX3_Control <- FindNeighbors(PDX3_Control, dims = 1:15)

PDX3_Control <- FindClusters(PDX3_Control, resolution = 1)

PDX3_Control <- RunUMAP(PDX3_Control, dims = 1:15)


DimPlot(PDX3_Control, reduction = "umap", label = TRUE)

head(PDX3_Control)

Idents(PDX3_Control) <- "seurat_clusters"

DimPlot(PDX3_Control, reduction = "umap", label = TRUE)


# Extract UMAP data and the meta data for the 'mean_expression_F3_TCGA' column
umap_data <- FetchData(PDX3_Control, 
                       vars = c("ESR1", "umap_1", "umap_2"))

# Now plot using ggplot2 directly with specified color limits
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = ESR1)) +
  geom_point(alpha = 0.5, size = 2) +  # Size can be adjusted to be smaller or larger
  scale_color_viridis_c(option = "C", 
                        limits = c(0, 3.5),  # Replace with your chosen min and max values
                        na.value = "lightgrey", 
                        oob = scales::squish) +
  labs(color = "ESR1 Expression") +
  theme_minimal()


############################################################




###########################3 Iteration of configurations to calculate correlations


# Load required libraries
library(Seurat)
library(UCell)
library(dplyr)

head(PDX3_Control)

# =============================================================================
# STEP 2: Modified Metacell Generation Function
# =============================================================================
# This function aggregates both gene expression (from the provided assay) 
# and the pre-computed signature scores (numeric columns from meta.data).
generate_metacells_with_signatures <- function(
    seurat_object,
    assay = "RNA",
    cells_per_metacell,
    reduction = "pca",
    dims_to_use = 1:15,
    verbose = TRUE
) {
  # Ensure the specified reduction exists in the Seurat object
  if (!(reduction %in% names(seurat_object@reductions))) {
    stop(paste("Reduction method", reduction, "not found in the Seurat object."))
  }
  
  # Shuffle cells and split them into non-overlapping groups
  all_cells <- Cells(seurat_object)
  shuffled_cells <- sample(all_cells)
  metacell_groups <- split(shuffled_cells, ceiling(seq_along(shuffled_cells) / cells_per_metacell))
  
  if (verbose) {
    cat("Number of metacells generated:", length(metacell_groups), "\n")
  }
  
  # Extract expression data (assumed to be in the "data" slot; adjust if needed)
  expression_data <- GetAssayData(seurat_object, assay = assay, slot = "data")
  
  # Identify signature score columns from meta.data.
  # Here we assume that signature scores are numeric and contain "_UCell" in their names.
  # You can adjust the pattern if needed.
  metadata <- seurat_object@meta.data
  signature_cols <- grep("_UCell", colnames(metadata), value = TRUE)
  
  if (verbose) {
    cat("Aggregating signature score columns:", paste(signature_cols, collapse = ", "), "\n")
  }
  
  # Initialize lists for aggregated data
  aggregated_expression <- vector("list", length = length(metacell_groups))
  aggregated_signature_scores <- vector("list", length = length(metacell_groups))
  
  for (i in seq_along(metacell_groups)) {
    cells_in_group <- metacell_groups[[i]]
    
    # Aggregate expression data: compute mean for each gene
    aggregated_expression[[i]] <- rowMeans(expression_data[, cells_in_group, drop = FALSE])
    
    # Aggregate signature scores: compute mean for each signature column
    aggregated_signature_scores[[i]] <- colMeans(metadata[cells_in_group, signature_cols, drop = FALSE], na.rm = TRUE)
  }
  
  # Combine the aggregated data into matrices
  expr_mat <- do.call(cbind, aggregated_expression)
  colnames(expr_mat) <- paste0("metacell_", seq_along(metacell_groups))
  
  sig_mat <- do.call(rbind, aggregated_signature_scores)
  rownames(sig_mat) <- paste0("metacell_", seq_along(metacell_groups))
  
  if (verbose) {
    cat("Dimensions of aggregated expression matrix:", dim(expr_mat), "\n")
    cat("Dimensions of aggregated signature matrix:", dim(sig_mat), "\n")
  }
  
  return(list(
    metacells_expression = expr_mat,
    metacells_signature_scores = sig_mat
  ))
}


# =============================================================================
# STEP 3: Grid Search Over Metacell Configurations and Correlation Computation
# =============================================================================

# Define parameter grid
reductions <- c("pca", "umap")
cells_per_metacell_values <- c(2, 5, 10, 20, 50)
dims_values <- c(10, 15, 20, 30)

# Initialize list to store results
results_list <- list()

for (red in reductions) {
  for (cells in cells_per_metacell_values) {
    for (dims in dims_values) {
      
      cat("Processing config: reduction =", red, 
          ", cells_per_metacell =", cells, 
          ", dims =", dims, "\n")
      
      # Generate metacells using the modified function
      metacell_results <- generate_metacells_with_signatures(
        seurat_object = PDX3_Control,
        assay = "SCT",  # Adjust the assay if needed
        cells_per_metacell = cells,
        reduction = red,
        dims_to_use = 1:dims,
        verbose = FALSE
      )
      
      expr_mat <- metacell_results$metacells_expression
      sig_mat <- metacell_results$metacells_signature_scores
      
      # Check that ESR1 is present in the aggregated expression matrix
      if (!"ESR1" %in% rownames(expr_mat)) {
        cat("ESR1 not found in aggregated expression for config:", red, cells, dims, "\n")
        next
      }
      
      # Extract aggregated ESR1 expression (vector of length = number of metacells)
      esr1_metacell <- expr_mat["ESR1", ]
      
      # Compute correlations for each signature
      cor_results <- sapply(colnames(sig_mat), function(sig_col) {
        cor(esr1_metacell, sig_mat[, sig_col], method = "pearson", use = "complete.obs")
      })
      
      # Create a data frame for this configuration
      config_df <- data.frame(
        reduction = red,
        cells_per_metacell = cells,
        dims = dims,
        signature = names(cor_results),
        correlation = cor_results,
        stringsAsFactors = FALSE
      )
      
      results_list[[paste(red, cells, dims, sep = "_")]] <- config_df
    }
  }
}

# Combine all results into a single data frame
final_results <- bind_rows(results_list)

# Save results to CSV
write.csv(final_results, "metacell_correlation_PDX3_Control_ESR1_expression_to_kinases_results_aggregated_scores.csv", row.names = FALSE)
cat("Grid search complete. Results saved to metacell_correlation_results_aggregated_scores.csv\n")



# Extract the metadata from the Seurat object
meta_data <- PDX3_Control@meta.data

# Print all metadata column names
cat("Metadata column names in PDX3_Control:\n")
print(colnames(meta_data))

# Optionally, display the first few rows of the metadata
cat("\nFirst few rows of the metadata:\n")
print(head(meta_data))




# Extract metadata from the Seurat object
meta_data <- PDX3_Control@meta.data

# Get all metadata column names
all_cols <- colnames(meta_data)

# Define columns to exclude
exclude_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA")

# Extract signature columns by excluding the non-signature ones
signature_cols <- setdiff(all_cols, exclude_cols)

# Print the signature column names
cat("Signature columns found:\n")
print(signature_cols)













library(Seurat)
library(dplyr)

# --- Modified metacell generation function that filters to numeric signature columns ---
generate_metacells_with_signatures <- function(
    seurat_object,
    assay = "RNA",
    cells_per_metacell,
    reduction = "pca",
    dims_to_use = 1:15,
    verbose = TRUE
) {
  # Ensure the specified reduction exists in the Seurat object
  if (!(reduction %in% names(seurat_object@reductions))) {
    stop(paste("Reduction method", reduction, "not found in the Seurat object."))
  }
  
  # Shuffle cells and split them into non-overlapping groups
  all_cells <- Cells(seurat_object)
  shuffled_cells <- sample(all_cells)
  metacell_groups <- split(shuffled_cells, ceiling(seq_along(shuffled_cells) / cells_per_metacell))
  
  if (verbose) {
    cat("Number of metacells generated:", length(metacell_groups), "\n")
  }
  
  # Extract expression data (using the specified assay and the "data" slot)
  expression_data <- GetAssayData(seurat_object, assay = assay, slot = "data")
  
  # Extract metadata and then select signature columns by excluding non-signature ones.
  metadata <- seurat_object@meta.data
  exclude_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
  signature_cols <- setdiff(colnames(metadata), exclude_cols)
  
  # Only keep numeric columns
  signature_cols <- signature_cols[sapply(metadata[, signature_cols, drop = FALSE], is.numeric)]
  
  if (verbose) {
    cat("Aggregating signature score columns:\n", paste(signature_cols, collapse = ", "), "\n")
  }
  
  # Initialize lists for aggregated data
  aggregated_expression <- vector("list", length = length(metacell_groups))
  aggregated_signature_scores <- vector("list", length = length(metacell_groups))
  
  for (i in seq_along(metacell_groups)) {
    cells_in_group <- metacell_groups[[i]]
    
    # Aggregate expression: compute mean for each gene
    aggregated_expression[[i]] <- rowMeans(expression_data[, cells_in_group, drop = FALSE])
    
    # Aggregate signature scores: compute mean for each signature column
    aggregated_signature_scores[[i]] <- colMeans(metadata[cells_in_group, signature_cols, drop = FALSE], na.rm = TRUE)
  }
  
  # Combine aggregated data into matrices
  expr_mat <- do.call(cbind, aggregated_expression)
  colnames(expr_mat) <- paste0("metacell_", seq_along(metacell_groups))
  
  sig_mat <- do.call(rbind, aggregated_signature_scores)
  rownames(sig_mat) <- paste0("metacell_", seq_along(metacell_groups))
  
  if (verbose) {
    cat("Dimensions of aggregated expression matrix:", dim(expr_mat), "\n")
    cat("Dimensions of aggregated signature matrix:", dim(sig_mat), "\n")
  }
  
  return(list(
    metacells_expression = expr_mat,
    metacells_signature_scores = sig_mat
  ))
}

# --- Grid search parameters and correlation calculation ---

# Define parameter grid
reductions <- c("pca", "umap")
cells_per_metacell_values <- c(1, 5, 10, 20, 50)
dims_values <- c(15, 20, 30)

# Initialize list to store results
results_list <- list()

for (red in reductions) {
  for (cells in cells_per_metacell_values) {
    for (dims in dims_values) {
      
      cat("Processing config: reduction =", red, 
          ", cells_per_metacell =", cells, 
          ", dims =", dims, "\n")
      
      # Generate metacells using the modified function
      metacell_results <- generate_metacells_with_signatures(
        seurat_object = PDX3_Control,
        assay = "SCT",  # Adjust the assay if needed
        cells_per_metacell = cells,
        reduction = red,
        dims_to_use = 1:dims,
        verbose = FALSE
      )
      
      expr_mat <- metacell_results$metacells_expression
      sig_mat <- metacell_results$metacells_signature_scores
      
      # Check that ESR1 is present in the aggregated expression matrix
      if (!"ESR1" %in% rownames(expr_mat)) {
        cat("ESR1 not found in aggregated expression for config:", red, cells, dims, "\n")
        next
      }
      
      # Extract aggregated ESR1 expression (vector of length = number of metacells)
      esr1_metacell <- expr_mat["ESR1", ]
      
      # Compute Pearson correlations for each signature column
      cor_results <- sapply(colnames(sig_mat), function(sig_col) {
        cor(esr1_metacell, sig_mat[, sig_col], method = "pearson", use = "complete.obs")
      })
      
      # Create a data frame for this configuration
      config_df <- data.frame(
        reduction = red,
        cells_per_metacell = cells,
        dims = dims,
        signature = names(cor_results),
        correlation = cor_results,
        stringsAsFactors = FALSE
      )
      
      results_list[[paste(red, cells, dims, sep = "_")]] <- config_df
    }
  }
}

# Combine all results into a single data frame
final_results <- bind_rows(results_list)

# Save results to CSV
write.csv(final_results, "metacell_correlation_results_aggregated_scores.csv", row.names = FALSE)
cat("Grid search complete. Results saved to metacell_correlation_results_aggregated_scores.csv\n")





