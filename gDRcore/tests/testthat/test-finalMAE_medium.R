test_that("medium: test_synthetic_data", {
  data <- "finalMAE_medium.RDS"
  original <- get_synthetic_data(data)
  
  mae <- purrr::quietly(gDRtestData::generateMediumData)(
    cell_lines, drugs, e_inf, ec50, hill_coef
  )
  expect_length(mae$warnings, 3)
  
  test_synthetic_data(original, mae$result, data)
})
