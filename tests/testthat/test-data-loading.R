# Tests for data loading functions - force execution for coverage

test_that("hBaronDataset exists and is callable", {
  expect_true(exists("hBaronDataset"))
  expect_type(hBaronDataset, "closure")
})

test_that("muraro_dataset exists and is callable", {
  expect_true(exists("muraro_dataset"))
  expect_type(muraro_dataset, "closure")
})

test_that("hBaronDataset returns correct structure", {
  # DON'T skip - run it for coverage
  baron <- hBaronDataset()
  
  # Check structure
  expect_type(baron, "list")
  expect_true("data" %in% names(baron))
  expect_true("celltypes" %in% names(baron))
  
  # Check types
  expect_true(is.matrix(baron$data) || inherits(baron$data, "Matrix"))
  expect_true(is.vector(baron$celltypes) || is.factor(baron$celltypes))
  
  # Check dimensions match
  expect_equal(ncol(baron$data), length(baron$celltypes))
})

test_that("muraro_dataset returns correct structure", {
  # DON'T skip - run it for coverage
  muraro <- muraro_dataset()
  
  # Check structure
  expect_type(muraro, "list")
  expect_true("data" %in% names(muraro))
  expect_true("celltypes" %in% names(muraro))
  
  # Check types
  expect_true(is.matrix(muraro$data) || inherits(muraro$data, "Matrix"))
  expect_true(is.vector(muraro$celltypes) || is.factor(muraro$celltypes))
  
  # Check dimensions match
  expect_equal(ncol(muraro$data), length(muraro$celltypes))
})