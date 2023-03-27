test_that("small: test_synthetic_data", {
  data <- "finalMAE_small.RDS"
  original <- get_synthetic_data(data)
  
  mae <- purrr::quietly(gDRtestData::generateNoiseRawData)(
    cell_lines, drugs, FALSE
  )
  expect_length(mae$warnings, 2)
  
  mae$result
  
  
  test_synthetic_data(original, mae$result, data)
})


original_se <- original[[1]]
new_se <- mae$result[[1]]
