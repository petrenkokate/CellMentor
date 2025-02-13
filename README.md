# CellMentor

CellMentor is a novel supervised cell type aware non-negative matrix factorization (NMF) method designed for enhanced cell type resolution in single-cell RNA sequencing analysis. It enables improved cell type separation, clustering, and annotation by leveraging latent patterns from reference datasets.

## Installation

You can install CellMentor directly from GitHub using devtools:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("petrenkokate/CellMentor", dependencies = TRUE)
```

## Basic Usage

### 1. Create a CellMentor Object

```R
# Load required libraries
library(CellMentor)

# Create CellMentor object
object <- CreateCSFNMFobject(
    ref_matrix = reference_counts,    # Reference expression matrix (genes × cells)
    ref_celltype = reference_labels,  # Reference cell type labels
    data_matrix = query_counts        # Query expression matrix (genes × cells)
)
```

### 2. Select Optimal Parameters

The parameter selection process involves optimizing multiple hyperparameters to achieve the best cell type separation. Note that the initial rank (k) selection is computationally intensive and may take some time, but subsequent steps are relatively fast.

```R
# Find optimal parameters
optimal_params <- select_optimal_parameters(
    object,
    num_cores = 10,           # Number of cores for parallel processing
    subset_size = 1,          # Proportion of data to use (1 = all data)
    init_methods = c("uniform", "regulated", "NNDSVD", "skmeanGenes", "skmeanCells"), # initialization methods
    alpha_range = c(1),       # Range of alpha parameters to test
    beta_range = c(1),        # Range of beta parameters to test
    gamma_range = c(1),       # Range of gamma parameters to test
    delta_range = c(1)        # Range of delta parameters to test
)

# Get best model
final_model <- optimal_params$best_model
```

### 3. Project Data

```R
# Project query data onto the learned space
h_test <- project_data(
    W = final_model@W,                 # Learned gene weights
    X = final_model@matrices@data,     # Query data matrix
    seed = 1,
    num_cores = 10,
    chunk_size = NULL,
    verbose = TRUE
)
```


### 4. Integration with Seurat

```R
# Add CellMentor dimensionality reduction to Seurat object
seurat_object$CellMentor <- CreateDimReducObject(
    embeddings = t(as.matrix(h_test)),
    key = "CellMentor_",
    assay = DefaultAssay(seurat_object),
    loadings = as.matrix(final_model@W)
)

# Visualization
seurat_object <- RunUMAP(seurat_object, reduction = 'CellMentor', dims= 1:optimal_params$best_params$k)
DimPlot(seurat_object)
```
