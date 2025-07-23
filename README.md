# CellMentor

CellMentor is a novel supervised cell type aware non-negative matrix factorization (NMF) method designed for enhanced cell type resolution in single-cell RNA sequencing analysis. By integrating cell type annotations into the NMF framework, CellMentor enables improved cell type separation, clustering, and annotation by leveraging latent patterns from reference datasets.
Features

- Improved cell type separation through constrained supervised factorization (CSFNMF)
- Automated parameter optimization for optimal performance
- Efficient projection of query datasets onto learned cell type spaces
- Seamless integration with Seurat for visualization and downstream analysis

## Installation

You can install CellMentor directly from GitHub using devtools:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("petrenkokate/CellMentor", dependencies = TRUE)
```

Load required packages

```R
library(CellMentor)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Seurat)
```

## Basic Usage

Set random seed for reproducibility

```R
set.seed(100)
```

### 1. Load the datasets

```R
# Loading reference dataset (Baron)
baron <- h.baron_dataset()
reference_matrix <- baron$data
reference_celltypes <- baron$celltypes

# Loading query dataset (Muraro)
muraro <- muraro_dataset()
query_matrix <- muraro$data
query_celltypes <- muraro$celltypes # This would be unknown in a real application
                                            # We keep it here for evaluation
```

### (Optional) Create smaller subsets for faster demonstration

```R
# Function to create balanced subsets
create_subset <- function(matrix, celltypes, cells_per_type = 30) {
  # Get unique cell types
  unique_types <- unique(celltypes)
  
  # Select cells for each type
  selected_cells <- c()
  for (cell_type in unique_types) {
    # Get cells of this type
    type_cells <- names(celltypes)[celltypes == cell_type]
    
    # If fewer cells than requested, take all of them
    n_to_select <- min(cells_per_type, length(type_cells))
    
    # Randomly select cells
    selected <- sample(type_cells, n_to_select)
    selected_cells <- c(selected_cells, selected)
  }
  
  # Return subset
  list(
    matrix = matrix[, selected_cells],
    celltypes = celltypes[selected_cells]
  )
}

# Create balanced subsets with 30 cells per type
baron_subset <- create_subset(reference_matrix, reference_celltypes, 30)
muraro_subset <- create_subset(query_matrix, query_celltypes, 30)

# Update variable names for clarity
reference_matrix <- baron_subset$matrix
reference_celltypes <- baron_subset$celltypes
query_matrix <- muraro_subset$matrix
query_celltypes <- muraro_subset$celltypes # This would be unknown in a real application
                                            # We keep it here for evaluation
```

### 2. Create a CellMentor Object

```R
csfnmf_obj <- CreateCSFNMFobject(
  ref_matrix = reference_matrix,
  ref_celltype = reference_celltypes,
  data_matrix = query_matrix,
  norm = TRUE,
  most.variable = TRUE,
  scale = TRUE,
  scale_by = "cells",
  verbose = TRUE,
  num_cores = 1
)
```

***Expected Output::***
```
[11:45:56] Starting CSFNMF object creation
[11:45:59] Validating inputs
[11:45:59] Creating CSFNMF object
[11:45:59] Converting matrices to sparse format
[11:46:00] Setting up annotations
[11:46:00] Cleaning matrices
[11:46:00] Removed 3319 empty rows from reference matrix
[11:46:00] Removed 496 empty rows from query matrix
[11:46:00] Finding common genes
[11:46:00] Found 14904 common genes
[11:46:00] Normalizing data
[11:46:04] Selecting variable genes
[11:46:05] Selected 999 variable genes
[11:46:05] Scaling data by cells
[11:46:06] Encoding cell types
[11:46:06] Reordering data
[11:46:06] Validating final object
[11:46:06] CSFNMF object creation complete
```

### 3. Run CellMentor and choose optimal params

The parameter selection process involves optimizing multiple hyperparameters to achieve the best cell type separation. Note that the function is computationally intensive and because of it we test only limited ranges of parameters.

#### Parameter descriptions:
 - alpha_range: Controls within-class scatter (cell similarity within the same type)
 - beta_range: Controls between-class scatter (cell separation between different types)
 - gamma_range: Controls sparsity of the factorization
 - delta_range: Controls orthogonality between factors

this is parameter settings that work well across most datasets without extensive tuning.
```R
# Find optimal parameters
optimal_params <- CellMentor(
  csfnmf_obj,
  alpha_range = c(1, 5),      # Limited alpha range
  beta_range = c(1, 5),       # Limited beta range
  gamma_range = c(0.1),     # use only one gamma for speed
  delta_range = c(1),       # use only one delta for speed
  num_cores = 5,
  verbose = TRUE
)

# Get best model
best_model <- optimal_params$best_model
K_VALUE <- best_model@parameters$rank
```

***Expected Output:***

```
[11:46:11] Creating training object
[11:46:11] Determining optimal rank
[11:47:02] Optimal rank determined: 80
[11:47:02] Starting parameter grid search
[11:47:02] Testing configuration 1/4
[11:49:09] Testing configuration 2/4
[11:51:20] Testing configuration 3/4
[11:53:25] Testing configuration 4/4
[11:55:40] Training final model with best parameters on full dataset
[11:55:40] Initializing W and H matrices
[11:55:40] Calculating helper matrices
[11:55:40] Calculating alpha
[11:55:40] Calculating H constants
[11:55:47] Updating W and H matrices
```

Check best params

```R
print(optimal_params$best_params)
```
***Expected Output:***

```
$k
[1] 80

$init_method
[1] "regulated"

$alpha
[1] 5

$beta
[1] 5

$gamma
[1] 0.1

$delta
[1] 1
```

### 3. Project Data

```R
# Project query data onto the learned space
h_project <- project_data(
  W = best_model@W, # Learned gene weights
  X = best_model@matrices@data,  # Query data matrix
  num_cores = 5,
  verbose = TRUE
)
```

### 4. Integration with Seurat

```R
rownames(query_matrix) <- make.unique(rownames(query_matrix))
seu_muraro <- CreateSeuratObject(counts = query_matrix)
seu_muraro$celltype <- query_celltypes
# Add CellMentor dimensionality reduction to Seurat object
seu_muraro$CellMentor <- CreateDimReducObject(
    embeddings = t(as.matrix(h_project)),
    key = "CellMentor_",
    assay = DefaultAssay(seu_muraro),
    loadings = as.matrix(best_model@W)
)

# Visualization
seu_muraro <- RunUMAP(seu_muraro, reduction = 'CellMentor', dims= 1:K_VALUE)
DimPlot(seu_muraro, group.by = 'celltype')
```
