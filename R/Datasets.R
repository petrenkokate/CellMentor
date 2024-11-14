#' Load Baron Human Pancreas Dataset
#' 
#' Loads and processes the Baron et al. human pancreas single-cell RNA-seq dataset
#' 
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom scRNAseq BaronPancreasData
#' @export
h.baron_dataset <- function() {
  # Load raw data
  baron.h <- scRNAseq::BaronPancreasData(which = 'human')
  
  # Extract and process cell type annotations
  celltype_mappings <- c(
    "activated_stellate" = "activated stellate",
    "quiescent_stellate" = "quiescent stellate",
    "t_cell" = "T cell"
  )
  
  baron.h.celltype <- baron.h[["label"]]
  names(baron.h.celltype) = baron.h@colData@rownames
  
  # Update cell type labels using vectorized operation
  for (old_type in names(celltype_mappings)) {
    baron.h.celltype[baron.h.celltype == old_type] <- celltype_mappings[old_type]
  }
  
  # Process expression data
  baron.h.data <- as.matrix(counts(baron.h))
  rownames(baron.h.data) <- gsub(".", "-", rownames(baron.h.data), fixed = TRUE)
  
  list(
    data = baron.h.data,
    celltypes = baron.h.celltype
  )
}

#' Load Muraro Pancreas Dataset
#' 
#' Loads and processes the Muraro et al. pancreas single-cell RNA-seq dataset
#' 
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom scRNAseq MuraroPancreasData
#' @export
muraro_dataset <- function() {
  # Load and filter data
  sceM <- scRNAseq::MuraroPancreasData()
  sceM <- sceM[, !is.na(sceM$label) & sceM$label != "unclear"]
  
  # Process expression data
  sceM.data <- as.matrix(counts(sceM))
  rownames(sceM.data) <- gsub("\\__chr[0-9]*", "", rownames(sceM.data))
  rownames(sceM.data) <- gsub(".", "-", rownames(sceM.data), fixed = TRUE)
  
  # Process cell types
  celltype_mappings <- c(
    "pp" = "gamma",
    "duct" = "ductal"
  )
  
  sceM.celltype <- sceM$label
  for (old_type in names(celltype_mappings)) {
    sceM.celltype[sceM.celltype == old_type] <- celltype_mappings[old_type]
  }
  
  # Remove mesenchymal cells
  keep_cells <- sceM.celltype != "mesenchymal"
  sceM.data <- sceM.data[, keep_cells]
  sceM.celltype <- sceM.celltype[keep_cells]
  
  names(sceM.celltype) <- colnames(sceM.data)
  
  list(
    data = sceM.data,
    celltypes = sceM.celltype
  )
}

#' Load PBMC Dataset from SeuratData
#' 
#' Loads and processes the PBMC single-cell RNA-seq dataset 
#' 
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom SeuratData InstallData
#' @importFrom Seurat UpdateSeuratObject PercentageFeatureSet
#' @export
pbmcsca_data <- function() {
  
  # Define constants
  ACCEPTED_METHODS <- c(
    "10x Chromium (v2) A",
    "10x Chromium (v2) B",
    "10x Chromium (v3)",
    "10x Chromium (v2)"
  )
  
  CELL_TYPES_TO_KEEP <- c(
    "B cell",
    "CD14+ monocyte",
    "CD4+ T cell",
    "Cytotoxic T cell",
    "Natural killer cell"
  )
  
  # Load and process data
  utils::data("pbmcsca", envir = environment())
  if (!exists("pbmcsca")) {
    stop("Could not load zheng_data. Please ensure the package is properly installed.")
  }
  
  # Quality control
  pbmcsca[["percent.mt"]] <- PercentageFeatureSet(pbmcsca, pattern = "^MT-")
  pbmcsca <- subset(pbmcsca,
                    subset = Method %in% ACCEPTED_METHODS & 
                      nFeature_RNA < 4000 & 
                      percent.mt < 5
  )
  
  # Extract data and cell types
  pbmcsca.data <- pbmcsca[["RNA"]]@counts
  pbmcsca.celltype <- pbmcsca$CellType
  names(pbmcsca.celltype) <- colnames(pbmcsca.data)
  rm(pbmcsca)
  # Filter for specific cell types
  keep_cells <- pbmcsca.celltype %in% CELL_TYPES_TO_KEEP
  pbmcsca.data <- pbmcsca.data[, keep_cells]
  pbmcsca.celltype <- pbmcsca.celltype[keep_cells]
  
  list(
    data = pbmcsca.data,
    celltypes = pbmcsca.celltype
  )
}


#' Load Zeisel Brain Dataset
#' 
#' Loads and processes the Zeisel et al. brain single-cell RNA-seq dataset
#' 
#' @param level Integer indicating the annotation level (1 or 2)
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom scRNAseq ZeiselBrainData
#' @export
zeisel_dataset <- function(level = 1) {
  # Load raw data
  zeisel.brain <- scRNAseq::ZeiselBrainData()
  
  if (level == 1) {
    # Level 1 classification with manual corrections
    cell_type_corrections <- c(
      "6" = "Endothelial Cell",
      "9" = "Mural",
      "7" = "Astrocyte",
      "8" = "Ependymal"
    )
    
    zeisel.brain.celltype <- zeisel.brain$level1class
    # Apply corrections based on group numbers
    for (group_num in names(cell_type_corrections)) {
      zeisel.brain.celltype[zeisel.brain$`group #` == as.numeric(group_num)] <- 
        cell_type_corrections[group_num]
    }
    
  } else {
    # Level 2 classification
    zeisel.brain.celltype <- zeisel.brain$level2class
    # Remove cells with undefined types
    valid_cells <- zeisel.brain.celltype != "(none)"
    zeisel.brain.celltype <- zeisel.brain.celltype[valid_cells]
  }
  
  # Extract expression data
  zeisel.brain.data <- counts(zeisel.brain)
  
  # Filter data if using level 2 classification
  if (level != 1) {
    zeisel.brain.data <- zeisel.brain.data[, valid_cells]
  }
  
  # Set cell names
  names(zeisel.brain.celltype) <- zeisel.brain$cell_id
  
  list(
    data = zeisel.brain.data,
    celltypes = zeisel.brain.celltype
  )
}


#' Load Tasic Brain Dataset
#' 
#' Loads and processes the Tasic et al. mouse brain single-cell RNA-seq dataset 
#' 
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom scRNAseq TasicBrainData
#' @importFrom scater addPerCellQC quickPerCellQC calculateTPM
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom GenomicRanges reduce
#' @export
tasic_dataset <- function() {
  # Load initial dataset
  tasic.brain <- scRNAseq::TasicBrainData()
  
  # Quality control
  tasic.brain <- scater::addPerCellQC(
    tasic.brain,
    subsets = list(mito = grep("^mt_", rownames(tasic.brain)))
  )
  
  # Filter cells with missing ERCC data
  valid_ercc <- !is.na(tasic.brain$altexps_ERCC_percent)
  tasic.brain <- tasic.brain[, valid_ercc]
  
  # Apply QC metrics
  qc_metrics <- scater::quickPerCellQC(
    colData(tasic.brain),
    percent_subsets = c("subsets_mito_percent", "altexps_ERCC_percent")
  )
  tasic.brain <- tasic.brain[, !qc_metrics$discard]
  
  # Get transcript length information
  ah <- AnnotationHub()
  mm.db <- ah[["AH73905"]]
  mm.exons <- exonsBy(mm.db, by = "gene")
  mm.exons <- GenomicRanges::reduce(mm.exons)
  mm.len <- sum(width(mm.exons))
  mm.symb <- mapIds(mm.db, 
                    keys = names(mm.len), 
                    keytype = "GENEID", 
                    column = "SYMBOL")
  names(mm.len) <- mm.symb
  
  # Filter and normalize data
  common_genes <- intersect(names(mm.len), rownames(tasic.brain))
  tasic.brain <- tasic.brain[common_genes, ]
  
  # Calculate TPM
  tpm_data <- scater::calculateTPM(tasic.brain, lengths = mm.len[common_genes])
  
  # Process cell types
  tasic.brain.celltype <- tasic.brain$broad_type
  names(tasic.brain.celltype) <- tasic.brain$sample_title
  
  # Remove unclassified cells
  classified_cells <- tasic.brain.celltype != "Unclassified"
  tpm_data <- tpm_data[, classified_cells]
  tasic.brain.celltype <- tasic.brain.celltype[classified_cells]
  
  # Standardize cell type names
  tasic.brain.celltype[tasic.brain.celltype == "Microglia"] <- "microglia"
  
  list(
    data = tpm_data,
    celltypes = tasic.brain.celltype
  )
}


#' Load Zheng Dataset with Downsampling
#' 
#' Loads and processes the Zheng et al. immune cell dataset with optional downsampling
#' 
#' @param min_cell_pre_type Minimum number of cells per cell type (default: 600)
#' @param seed Random seed for reproducibility (default: 1)
#' @param min_cell Minimum cells for gene expression (default: 5)
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom Matrix rowSums
#' @export
zheng_dataset <- function(min_cell_pre_type = 600, seed = 1, min_cell = 5) {
  # Load the data
  # browser()
  utils::data("zheng_data", envir = environment())
  if (!exists("zheng_data")) {
    stop("Could not load zheng_data. Please ensure the package is properly installed.")
  }
  zheng <- zheng_data  # Assign loaded data to working variable
  
  # Update cell type annotations directly in meta.data
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("c34_filtered")] = "HSC"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("cd14_monocytes")] = "Monocytes"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("cd56nk")] = "NK cells"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("jurkat")] = "jurkat cancer cells"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("b_cells")] = "B cells"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("cytotoxic_ cells")] = "CD8+ cytotoxic cells"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("naive_cytotoxic")] = "CD8+ Naive T-cells"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("memory_t")] = "CD4+ Memory T-cells"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("regulatory_t")] = "CD4+ Treg"
  zheng@meta.data[["Annotations"]][zheng@meta.data[["Annotations"]] %in% c("naive_t")] = "CD4+ Naive T-cells"
  
  # Extract data and cell types
  zheng_data_matrix <- zheng@assays[["RNA"]]@counts
  zheng_celltype_list <- zheng@meta.data[["Annotations"]]
  names(zheng_celltype_list) <- colnames(zheng_data_matrix)
  
  # Initialize selected cells
  selected_cells <- c()
  num_cells_per_celltype <- table(zheng_celltype_list)
  
  # Sample cells for each cell type
  set.seed(seed)
  for (type in names(num_cells_per_celltype)) {
    cell_names <- names(zheng_celltype_list[which(zheng_celltype_list == type)])
    selected_names <- sample(cell_names, min(min_cell_pre_type, num_cells_per_celltype[type]))
    selected_cells <- append(selected_cells, selected_names)
  }
  
  # Handle zero-expression genes
  new_data <- zheng_data_matrix[, selected_cells, drop = FALSE]
  zero_ex_gene <- rownames(new_data)[which(Matrix::rowSums(new_data) == 0)]
  
  while (length(zero_ex_gene) > 0) {
    gene <- zero_ex_gene[1]
    data_for_selected_cells <- zheng_data_matrix[, -which(colnames(zheng_data_matrix) %in% selected_cells)]
    nonzero_cell_names <- names(which(data_for_selected_cells[gene, ] > 0))
    
    if (length(nonzero_cell_names) <= 3) {
      break
    }
    
    set.seed(seed)
    cell_ex_gene <- sample(nonzero_cell_names, min(min_cell, length(nonzero_cell_names)))
    selected_cells <- append(selected_cells, cell_ex_gene)
    
    new_data <- zheng_data_matrix[, selected_cells]
    zero_ex_gene <- rownames(new_data)[which(Matrix::rowSums(new_data) == 0)]
  }
  
  # Prepare final dataset
  zheng_downsampled_data <- zheng_data_matrix[, selected_cells]
  zheng_downsampled_celltype <- zheng_celltype_list[selected_cells]
  
  list(
    data = zheng_downsampled_data,
    celltypes = zheng_downsampled_celltype
  )
}

#' Load CITE-seq Dataset
#' 
#' Loads and processes CITE-seq PBMC dataset with protein expression data
#' 
#' @return A list containing:
#'   \item{data}{RNA expression matrix}
#'   \item{ADT}{Antibody-derived tag expression matrix}
#'   \item{ADT scale}{Scaled antibody-derived tag data}
#' @importFrom Seurat CreateSeuratObject CreateAssayObject NormalizeData ScaleData Read10X_h5
#' @export
CITEseqDataset <- function() {
  # Create temporary directory for downloads
  temp_dir <- tempdir()
  h5_file <- file.path(temp_dir, "pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")
  
  # Download data
  tryCatch({
    download.file(
      "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5",
      destfile = h5_file
    )
    
    # Load and process data
    expression_data <- Seurat::Read10X_h5(filename = h5_file)
    pbmc <- CreateSeuratObject(counts = expression_data$`Gene Expression`, 
                               project = "pbmc")
    
    # Process protein expression data
    adt_counts <- expression_data$`Antibody Capture`
    common_cells <- intersect(colnames(adt_counts), colnames(pbmc))
    adt_counts <- adt_counts[, common_cells]
    
    # Create protein assay
    pbmc[["ADT"]] <- CreateAssayObject(counts = adt_counts)
    pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
    pbmc <- ScaleData(pbmc, assay = "ADT")
    
    # Extract relevant data
    result <- list(
      data = pbmc[["RNA"]]$counts,
      ADT = pbmc[["ADT"]]$data,
      "ADT scale" = pbmc[["ADT"]]$scale.data
    )
    
    return(result)
    
  }, finally = {
    # Clean up
    if (file.exists(h5_file)) {
      file.remove(h5_file)
    }
  })
}