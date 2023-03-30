test_that("wLigand: test_synthetic_data", {
  data <- "finalMAE_wLigand.RDS"
  original <- gDRutils::get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(gDRtestData::generateLigandData)(
    cell_lines, drugs, FALSE
  )
  expect_length(mae$warnings, 2)
  
  test_synthetic_data(original, mae$result, data)
})
