# Tests for CreateCSFNMFobject function

test_that("CreateCSFNMFobject creates valid object", {
  test_data <- create_minimal_test_data()
  
  obj <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = FALSE,
    most.variable = FALSE,
    scale = FALSE,
    verbose = FALSE,
    num_cores = 1
  )
  
  # Check it's the right class
  expect_s4_class(obj, "csfnmf")
  
  # Check matrices exist
  expect_true(!is.null(obj@matrices))
  expect_true(!is.null(obj@matrices@ref))
  expect_true(!is.null(obj@matrices@data))
})

test_that("CreateCSFNMFobject validates input dimensions", {
  test_data <- create_minimal_test_data()
  
  # Wrong celltype length should error
  wrong_celltypes <- test_data$ref_celltypes[1:10]
  
  expect_error(
    CreateCSFNMFobject(
      ref_matrix = test_data$ref_matrix,
      ref_celltype = wrong_celltypes,
      data_matrix = test_data$query_matrix,
      verbose = FALSE
    )
  )
})

test_that("CreateCSFNMFobject handles normalization", {
  test_data <- create_minimal_test_data()
  
  obj <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = TRUE,
    most.variable = FALSE,
    scale = FALSE,
    verbose = FALSE,
    num_cores = 1
  )
  
  # Should create valid object
  expect_s4_class(obj, "csfnmf")
})