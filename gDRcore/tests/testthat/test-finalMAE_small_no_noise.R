test_that("small_no_noise: test_synthetic_data", {
  data <- "finalMAE_small_no_noise.RDS"
  original <- get_synthetic_data(data)
  
  mae <- purrr::quietly(gDRtestData::generateNoNoiseRawData)(
    cell_lines, drugs, FALSE
  )
  expect_length(mae$warnings, 2)

  test_synthetic_data(original, mae$result, data)
})
