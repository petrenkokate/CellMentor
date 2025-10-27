# Tests for initialization methods in InitializeWH.R

test_that("initialize_regulated works", {
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
  
  # Run CellMentor with regulated initialization
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
  
  expect_equal(result$best_params$init_method, "regulated")
  expect_s4_class(result$best_model, "traincsfnmf")
})

test_that("initialize_uniform works", {
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
  
  # Run CellMentor with uniform initialization
  result <- CellMentor(
    obj,
    k = 5,
    init_methods = c("uniform"),
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  expect_equal(result$best_params$init_method, "uniform")
  expect_s4_class(result$best_model, "traincsfnmf")
})

test_that("initialize_nndsvd works", {
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
  
  # Run CellMentor with NNDSVD initialization
  result <- CellMentor(
    obj,
    k = 5,
    init_methods = c("NNDSVD"),
    alpha_range = c(1),
    beta_range = c(1),
    gamma_range = c(0.1),
    delta_range = c(1),
    verbose = FALSE,
    num_cores = 1
  )
  
  expect_equal(result$best_params$init_method, "NNDSVD")
  expect_s4_class(result$best_model, "traincsfnmf")
})