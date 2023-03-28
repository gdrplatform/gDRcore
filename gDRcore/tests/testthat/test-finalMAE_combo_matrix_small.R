test_that("combo_matrix_small: test_synthetic_data", {
  data <- "finalMAE_combo_matrix_small.RDS"
  original <- get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(gDRtestData::generateComboMatrixSmall)(
    cell_lines, drugs, FALSE
  )
  expect_length(mae$warnings, 4)

  test_synthetic_data(original, mae$result, data)
})
