testthat::context("Test assay-related functions")

testthat::test_that("'createSE' works as expected", {
  testDataDir <- system.file(package = "gDR", "testdata", "data7")
  lRef <- read_ref_data(testDataDir)
  
  # create SE containing 'matrix' assay
  se_matrix <-
    createSE(df_data = lRef$df_raw_data, assay_type = "matrix")
  expect_true(inherits(SummarizedExperiment::assay(se_matrix, 1), "matrix"))
  
  # create SE containing 'BumpyMatrix' assay
  se_bumpy_matrix <-
    createSE(df_data = lRef$df_raw_data, assay_type = "BumpyMatrix")
  expect_true(inherits(
    SummarizedExperiment::assay(se_bumpy_matrix, 1),
    "BumpyMatrix"
  ))
  
  # same dimensions of the assays in both SEs
  expect_identical(dim(se_matrix), dim(se_bumpy_matrix))
})
