skip("obsolete; refactored runDrugResponsePipeline2 as per GDR-621")

library(testthat)
context("Test update_experiment_metadata")

# get_cell_lines tests
test_that("update_experiment_metadata updates metadata properly", {
  testDataDir <- system.file(package = "gDR", "testdata", "data7")
  lRef <- read_ref_data(testDataDir)
  normSE <- gDR::normalize_SE(lRef$df_raw_data, discard_keys = 'Replicate')
  avgSE <- gDR::average_SE(normSE)
  metricsSE <- gDR::metrics_SE(avgSE)
  expect_error(S4Vectors::metadata(metricsSE) <- update_experiment_metadata(S4Vectors::metadata(metricsSE), list(description = "Testing description")), NA)
  expect_error(S4Vectors::metadata(metricsSE) <- update_experiment_metadata(S4Vectors::metadata(metricsSE), list(description = "Testing description")), NA)
})
