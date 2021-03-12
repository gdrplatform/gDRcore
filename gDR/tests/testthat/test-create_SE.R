library(testthat); library(gDR);
source("setUp.R")

test_that("create_SE2 works as expected", {
  create_SE2(test_df)
})


testthat::test_that("'create_SE' works as expected", {
  testDataDir <- system.file(package = "gDR", "testdata", "data7")
  lRef <- read_ref_data(testDataDir)
  
  # create SE containing 'matrix' assay
  se_matrix <-
    create_SE(df_data = lRef$df_raw_data, assay_type = "matrix")
  expect_true(inherits(SummarizedExperiment::assay(se_matrix, 1), "matrix"))
  
  # create SE containing 'BumpyMatrix' assay
  se_bumpy_matrix <-
    create_SE(df_data = lRef$df_raw_data, assay_type = "BumpyMatrix")
  expect_true(inherits(
    SummarizedExperiment::assay(se_bumpy_matrix, 1),
    "BumpyMatrix"
  ))
  
  # same dimensions of the assays in both SEs
  expect_identical(dim(se_matrix), dim(se_bumpy_matrix))
})

testthat::test_that("create_SE throwing expected errors", {
  # Test assertion:
  expect_error(
    create_SE(df_data = list(a = 1),
      "any(inherits(df_data, \"data.frame\"), inherits(df_data, \"DataFrame\")) is not TRUE",
      fixed = TRUE)
  )
  expect_error(
    create_SE(df_data = data.frame(a = 1), data_type = 1),
    "Assertion on 'data_type' failed: Must be of type 'character', not 'double'.",
    fixed = TRUE
  )
  expect_error(
    create_SE(
      df_data = data.frame(a = 1),
      data_type = "all",
      readout = c("a", "b")
    ),
    "Assertion on 'readout' failed: Must have length 1.",
    fixed = TRUE
  )
  expect_error(
    create_SE(
      df_data = data.frame(a = 1),
      data_type = "all",
      discard_keys = 1
    ),
    "Assertion on 'discard_keys' failed: Must be of type 'character' (or 'NULL'), not 'double'.",
    fixed = TRUE
  )
  expect_error(
    create_SE(
      df_data = data.frame(a = 1),
      data_type = "all",
      assay_type = "newMatrixObject"
    ),
    "'arg' should be one of \"matrix\", \"BumpyMatrix\"",
    fixed = TRUE
  )
})
