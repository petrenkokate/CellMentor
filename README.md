# CellMentor

CellMentor is a novel supervised cell type aware non-negative matrix factorization (NMF) method designed for enhanced cell type resolution in single-cell RNA sequencing analysis. By integrating cell type annotations into the NMF framework, CellMentor enables improved cell type separation, clustering, and annotation by leveraging latent patterns from reference datasets.
Features

- Improved cell type separation through constrained supervised factorization (CSFNMF)
- Automated parameter optimization for optimal performance
- Efficient projection of query datasets onto learned cell type spaces
- Seamless integration with Seurat for visualization and downstream analysis

# System Requirements
## Hardware requirements
- `CellMentor` package requires only a standard computer with enough RAM to support the in-memory operations.
- R version (e.g., R ≥ 4.3).
- Recommended RAM (e.g., ≥ 16 GB for medium datasets).
- GPU is NOT required.

## Software requirements
### OS Requirements
This package is supported for *macOS*, *Linux* and *Windows*. The package has been tested on the following systems:
+ macOS: Sequoia (15.6.1)
+ Linux: Ubuntu 20.04.6

### R Dependencies
`CellMentor` relies on several R packages, all of which are specified in the package DESCRIPTION file under the Imports field. These dependencies are installed automatically when you install CellMentor.

## Installation

You can install CellMentor from BioConductor:

```R
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("CellMentor")
```

You can install CellMentor directly from GitHub using devtools:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("petrenkokate/CellMentor", dependencies = TRUE)
```


### Typical Installation and Runtime

- **Installation time:** ~5 minutes on a standard desktop computer  
  (may vary depending on the number of dependencies to install).  

- **Demo runtime:** ~8–10 minutes on a standard desktop computer  
  (tested on macOS Sequoia 15.6.1 with R ≥ 4.3, 16 GB RAM).
  

## Basic Usage

Set random seed for reproducibility

```R
set.seed(100)
```

Load required packages

```R
library(CellMentor)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Seurat)
library(scRNAseq)
```

### 1. Load the datasets

```R
# Loading reference dataset (Baron)
baron <- hBaronDataset()
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

***Example Output (illustrative — exact counts depend on the current scRNAseq dataset snapshot):***
```
[hh:mm:ss] Starting CSFNMF object creation
[hh:mm:ss] Validating inputs
[hh:mm:ss] Creating CSFNMF object
[hh:mm:ss] Converting matrices to sparse format
[hh:mm:ss] Setting up annotations
[hh:mm:ss] Cleaning matrices
[hh:mm:ss] Removed ~5000 empty rows from reference matrix
[hh:mm:ss] Removed ~2500 empty rows from query matrix
[hh:mm:ss] Finding common genes
[hh:mm:ss] Found ~13000 common genes
[hh:mm:ss] Normalizing data
[hh:mm:ss] Selecting variable genes
[hh:mm:ss] Selected ~1500 variable genes
[hh:mm:ss] Scaling data by cells
[hh:mm:ss] Encoding cell types
[hh:mm:ss] Reordering data
[hh:mm:ss] Validating final object
[hh:mm:ss] CSFNMF object creation complete
```

### 3. Run CellMentor and choose optimal params

The parameter selection process involves optimizing multiple hyperparameters to achieve the best cell type separation. Note that the function is computationally intensive and because of it we test only limited ranges of parameters.

#### Parameter descriptions:
 - alpha_range: Controls within-class scatter (cell similarity within the same type)
 - beta_range: Controls between-class scatter (cell separation between different types)
 - gamma_range: Controls sparsity of the factorization
 - delta_range: Controls orthogonality between factors

These are parameter settings that work well across most datasets without extensive tuning.

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

***Example Output (illustrative — the selected rank depends on the current data snapshot; expect a value in the tens):***

```
[hh:mm:ss] Creating training object
[hh:mm:ss] Determining optimal rank
[hh:mm:ss] Optimal rank determined: <auto>
[hh:mm:ss] Starting parameter grid search
[hh:mm:ss] Testing configuration 1/4
[hh:mm:ss] Testing configuration 2/4
[hh:mm:ss] Testing configuration 3/4
[hh:mm:ss] Testing configuration 4/4
[hh:mm:ss] Training final model with best parameters on full dataset
[hh:mm:ss] Initializing W and H matrices
[hh:mm:ss] Calculating helper matrices
[hh:mm:ss] Calculating alpha
[hh:mm:ss] Calculating H constants
[hh:mm:ss] Updating W and H matrices
```

Check best params

```R
print(optimal_params$best_params)
```
***Example Output (illustrative — exact `$k` and `$alpha` depend on the data snapshot):***

```
$k
[1] <auto>

$init_method
[1] "regulated"

$alpha
[1] 1   # or 5; whichever scored higher NMI on the held-out split

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

# Run CellMentor on Your Own Data

This quick guide shows how to create a **CSFNMF** object from your data. After creating the object, **follow the same steps as in the demo** (parameter search with `CellMentor()`, projection with `project_data()`, optional Seurat integration).

## Inputs
- **Reference counts matrix**: genes × cells (`ref_counts`)
- **Reference annotations**: vector of length `ncol(ref_counts)` (`ref_celltypes`), names must match `colnames(ref_counts)`
- **Query counts matrix**: genes × cells (`qry_counts`), with overlapping gene IDs

> Tip: rows = genes, columns = cells. Use sparse matrices (`Matrix::dgCMatrix`) for speed/memory.

## Create the object

```r
library(Matrix)
library(CellMentor)

# 1) Build CSFNMF object
csfnmf_obj <- CreateCSFNMFobject(
  ref_matrix    = ref_counts,
  ref_celltype  = ref_celltypes,   # names(ref_celltypes) == colnames(ref_counts)
  data_matrix   = qry_counts,
  norm          = TRUE,
  most.variable = TRUE,
  scale         = TRUE,
  scale_by      = "cells",
  num_cores     = 1,
  verbose       = TRUE
)
```

## Next steps
Proceed exactly as in the **Demo** section:
1. **Hyperparameter search & training:** `optimal <- CellMentor(csfnmf_obj, ...)`  
2. **Best model:** `best_model <- optimal$best_model`  
3. **Projection:** `h_project <- project_data(W = best_model@W, X = best_model@matrices@data)`  
4. *(Optional)* **Seurat integration & UMAP:** same code as in the demo.

# Reproduce the analysis and figures from the paper

All scripts used to generate the analyses and figures reported in the manuscript are openly available at [petrenkokate/CellMentor_paper](https://github.com/petrenkokate/CellMentor_paper).
The repository contains:
- Code for data preprocessing, model training, and evaluation
- Scripts to reproduce each figure in the paper

# Citation

If you use *CellMentor* in your work, please cite:

Hevdeli O†, Petrenko E†, Aran D. CellMentor: cell-type aware dimensionality reduction for single-cell RNA-sequencing data. *Nat Commun* **17**, 396 (2026).  
doi: [https://doi.org/10.1038/s41467-025-67088-7](https://doi.org/10.1038/s41467-025-67088-7)

† These authors contributed equally to this work.

# License

This project is covered under the **Apache 2.0 License**.
