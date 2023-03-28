test_that("combo_matrix: test_synthetic_data", {
  data <- "finalMAE_combo_matrix.RDS"
  original <- get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(gDRtestData::generateComboMatrix)(cell_lines, drugs, FALSE)
  expect_length(mae$warnings, 4)
  test_synthetic_data(original, mae$result, data)
})
