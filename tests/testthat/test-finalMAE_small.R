test_that("small: test_synthetic_data", {
  data <- "finalMAE_small.RDS"
  original <- get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(gDRtestData::generateNoiseRawData)(
    cell_lines, drugs, FALSE
  )
  expect_length(mae$warnings, 2)
  
  test_synthetic_data(original, mae$result, data)
})
