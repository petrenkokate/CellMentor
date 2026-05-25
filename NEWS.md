# CellMentor 1.0.1 (2026-05-25)

## Bug fixes

* `calculate_performance()` now passes `as.factor(true_labels)` and
  `as.factor(clusters)` to `aricode::NMI()`. As of `aricode` 1.1.0,
  `NMI()` routes mixed-type inputs through `sort_pairs()`, which calls
  `as.integer()` on character labels — turning cell-type names such as
  `"alpha"` into `NA` and aborting with `"NA are not supported."`.
  Coercing both sides to factor keeps the computation well-defined.
  Backport of the fix shipped in 1.1.2 on the devel branch.

# CellMentor 0.99.2 (2024-12-10)

## CHANGES FOR BIOCONDUCTOR REVIEW

* Replaced `parallel` package with `BiocParallel` for Bioconductor compliance
* Renamed `h.baron_dataset()` to `hBaronDataset()` (old name deprecated)
* Updated R version requirement to >= 4.5.0
* Removed `GeneExpression` from biocViews (redundant with Transcriptomics)

# CellMentor 0.99.1 (2024-10-27)

* Initial Bioconductor submission of CellMentor.
