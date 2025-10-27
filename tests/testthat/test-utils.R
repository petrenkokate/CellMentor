# Tests for utility functions

test_that("utility functions handle basic inputs", {
  test_data <- create_minimal_test_data()
  
  # Test basic matrix operations if you have utility functions
  # Adjust based on what's actually in your utils.R
  
  # Example: if you have a function to check matrix validity
  # expect_true(is_valid_matrix(test_data$ref_matrix))
  
  # Example: if you have normalization utilities
  # normalized <- normalize_matrix(test_data$ref_matrix)
  # expect_true(all(Matrix::colSums(normalized) > 0))
})

test_that("scaling functions work correctly", {
  test_data <- create_minimal_test_data()
  
  # Test that scaling produces different output
  obj_scaled <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = FALSE,
    most.variable = FALSE,
    scale = TRUE,
    scale_by = "cells",
    verbose = FALSE
  )
  
  obj_unscaled <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = FALSE,
    most.variable = FALSE,
    scale = FALSE,
    verbose = FALSE
  )
  
  # Matrices should differ
  expect_false(identical(
    as.matrix(obj_scaled@matrices@ref),
    as.matrix(obj_unscaled@matrices@ref)
  ))
})