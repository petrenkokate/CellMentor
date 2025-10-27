# Tests for project_data function

test_that("project_data works with valid inputs", {
  skip_if_slow()
  
  test_data <- create_minimal_test_data()
  
  obj <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = FALSE,
    most.variable = FALSE,
    scale = FALSE,
    verbose = FALSE
  )
  
  result <- CellMentor(
    obj,
    k = 5,
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  # Project data
  h_project <- project_data(
    W = result$best_model@W,
    X = result$best_model@matrices@data,
    num_cores = 1,
    verbose = FALSE
  )
  
  # Check output is matrix-like
  expect_true(is.matrix(h_project) || inherits(h_project, "Matrix"))
  
  # Check dimensions: k x cells
  expect_equal(nrow(h_project), result$best_params$k)
  expect_equal(ncol(h_project), ncol(result$best_model@matrices@data))
})

test_that("project_data produces non-negative values", {
  skip_if_slow()
  
  test_data <- create_minimal_test_data()
  
  obj <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = FALSE,
    most.variable = FALSE,
    scale = FALSE,
    verbose = FALSE
  )
  
  result <- CellMentor(
    obj,
    k = 5,
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  h_project <- project_data(
    W = result$best_model@W,
    X = result$best_model@matrices@data,
    num_cores = 1,
    verbose = FALSE
  )
  
  # NMF should produce non-negative results
  expect_true(all(h_project >= 0))
})

test_that("project_data validates matrix dimensions", {
  test_data <- create_minimal_test_data()
  
  # Create mismatched matrices
  W <- Matrix::Matrix(runif(50 * 10), nrow = 50, ncol = 10)
  X_wrong <- Matrix::Matrix(runif(100 * 30), nrow = 100, ncol = 30)
  
  # Should error due to dimension mismatch
  expect_error(
    project_data(W = W, X = X_wrong, num_cores = 1, verbose = FALSE)
  )
})

test_that("project_data handles sparse matrices", {
  skip_if_slow()
  
  test_data <- create_minimal_test_data()
  
  obj <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = FALSE,
    most.variable = FALSE,
    scale = FALSE,
    verbose = FALSE
  )
  
  result <- CellMentor(
    obj,
    k = 5,
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  # Ensure X is sparse
  X_sparse <- Matrix::Matrix(result$best_model@matrices@data, sparse = TRUE)
  
  h_project <- project_data(
    W = result$best_model@W,
    X = X_sparse,
    num_cores = 1,
    verbose = FALSE
  )
  
  expect_true(!is.null(h_project))
  expect_equal(ncol(h_project), ncol(X_sparse))
})

test_that("project_data results are finite", {
  skip_if_slow()
  
  test_data <- create_minimal_test_data()
  
  obj <- CreateCSFNMFobject(
    ref_matrix = test_data$ref_matrix,
    ref_celltype = test_data$ref_celltypes,
    data_matrix = test_data$query_matrix,
    norm = FALSE,
    most.variable = FALSE,
    scale = FALSE,
    verbose = FALSE
  )
  
  result <- CellMentor(
    obj,
    k = 5,
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  h_project <- project_data(
    W = result$best_model@W,
    X = result$best_model@matrices@data,
    num_cores = 1,
    verbose = FALSE
  )
  
  # No NaN, Inf, or -Inf values
  expect_true(all(is.finite(h_project)))
})