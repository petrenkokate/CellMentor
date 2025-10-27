# Tests for CellMentor parameter optimization

test_that("CellMentor runs with minimal parameters", {
  skip_if_slow()
  
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
  
  result <- CellMentor(
    obj,
    k = 5,
    init_methods = c("regulated"),
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  # Check return structure
  expect_type(result, "list")
  expect_true("best_params" %in% names(result))
  expect_true("best_model" %in% names(result))
  
  # Check best_params has required fields
  expect_true("k" %in% names(result$best_params))
  expect_true("alpha" %in% names(result$best_params))
  expect_true("beta" %in% names(result$best_params))
  
  # Check best_model is valid (it's a traincsfnmf object)
  expect_s4_class(result$best_model, "traincsfnmf")
})

test_that("CellMentor explores single parameter combination", {
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
    init_methods = c("regulated"),
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  # Should complete successfully
  expect_type(result, "list")
  expect_true("best_params" %in% names(result))
})

test_that("CellMentor best_model has W and H matrices", {
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
  
  # Check W and H exist
  expect_true(!is.null(result$best_model@W))
  expect_true(!is.null(result$best_model@H))
  
  # Check W dimensions: genes x k
  expect_equal(nrow(result$best_model@W), nrow(obj@matrices@data))
  expect_equal(ncol(result$best_model@W), result$best_params$k)
  
  # Check H dimensions: k x cells
  expect_equal(nrow(result$best_model@H), result$best_params$k)
})

test_that("CellMentor validates k parameter", {
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
  
  # k = 0 should fail
  expect_error(
    CellMentor(
      obj,
      k = 0,
      alpha_range = c(1),
      beta_range = c(1),
      gamma_range = c(0.1),
      delta_range = c(1),
      verbose = FALSE
    )
  )
})