# CellMentor 1.1.1 (2026-05-19)

## Bug fixes

* `cluster_nmf()` now falls back to k-means clustering when the NMF
  embedding has too few features or cells for SNN (variable-feature
  selection on low-dim embeddings was failing under Seurat 5.x). The
  number of k-means / hclust centers is now clamped to a feasible range,
  so grid search no longer fails silently on small reference splits.
* `CellMentor()` now stops with a clear message if every parameter
  configuration in the grid search fails, instead of cascading to a
  cryptic `switch()` error during final-model training.
* `evaluate_params()` warnings now include the underlying error message,
  so build logs surface the real cause of failures.
* `initialize_wh()` validates that `method` is a single character string
  and reports a helpful error if not.

## Documentation

* Updated citation to the published Nature Communications paper
  (Hevdeli et al., 2026, doi:10.1038/s41467-025-67088-7).
* Added `inst/CITATION` so `citation("CellMentor")` returns a structured
  bibliographic entry.

# CellMentor 0.99.2 (2024-12-10)

## CHANGES FOR BIOCONDUCTOR REVIEW

* Replaced `parallel` package with `BiocParallel` for Bioconductor compliance
* Renamed `h.baron_dataset()` to `hBaronDataset()` (old name deprecated)
* Updated R version requirement to >= 4.5.0
* Removed `GeneExpression` from biocViews (redundant with Transcriptomics)

# CellMentor 0.99.1 (2024-10-27)

* Initial Bioconductor submission of CellMentor.
