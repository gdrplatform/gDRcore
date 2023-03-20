test_that("many_lines: test_synthetic_data", {
  data <- "finalMAE_many_lines.RDS"
  original <- get_synthetic_data(data)
  
  mae <-
    purrr::quietly(gDRtestData::generateManyLinesData)(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
  expect_length(mae$warnings, 2)
  
  test_synthetic_data(original, mae$result, data)
})
