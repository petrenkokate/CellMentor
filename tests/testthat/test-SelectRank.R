test_that("SelectRank works with minimal data", {
  skip_if_slow()
  
  test_data <- create_minimal_test_data()
  train_matrix <- test_data$query_matrix  # use plain matrix
  
  # Run SelectRank
  rank <- SelectRank(
    train_matrix = train_matrix,
    max_p_value = 0.05,
    jump_step = 0.2,
    verbose = FALSE
  )
  
  # Expect a numeric rank
  expect_type(rank, "integer")
  expect_true(rank > 0)
})

test_that("SelectRank validates input parameters", {
  test_data <- create_minimal_test_data()
  train_matrix <- test_data$query_matrix
  
  # invalid range should error
  expect_error(
    SelectRank(
      train_matrix = train_matrix,
      range = c(1, 0),
      verbose = FALSE
    )
  )
  
  # invalid matrix input should error
  expect_error(
    SelectRank(
      train_matrix = NULL,
      verbose = FALSE
    )
  )
})

test_that("SelectRank returns reasonable rank values", {
  skip_if_slow()
  
  test_data <- create_minimal_test_data()
  train_matrix <- test_data$query_matrix
  
  rank <- SelectRank(
    train_matrix = train_matrix,
    max_p_value = 0.01,
    jump_step = 0.1,
    verbose = FALSE
  )
  
  n_cells <- ncol(train_matrix)
  expect_lte(rank, n_cells)
  expect_gte(rank, 1)
})