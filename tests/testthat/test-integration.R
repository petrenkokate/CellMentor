# End-to-end integration tests

test_that("Complete workflow executes successfully", {
  skip_if_slow()
  
  test_data <- create_minimal_test_data()
  
  # Step 1: Create object
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
  
  expect_s4_class(obj, "csfnmf")
  
  # Step 2: Optimize parameters (minimal grid for speed)
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
  
  expect_type(result, "list")
  expect_true("best_model" %in% names(result))
  
  # Step 3: Project data
  h_project <- project_data(
    W = result$best_model@W,
    X = result$best_model@matrices@data,
    num_cores = 1,
    verbose = FALSE
  )
  
  expect_true(!is.null(h_project))
  expect_true(all(h_project >= 0))  # NMF non-negativity
})