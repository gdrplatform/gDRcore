test_that("combo_2dose_nonoise2: test_synthetic_data", {
  data <- "finalMAE_combo_2dose_nonoise2.RDS"
  original <- get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(gDRtestData::generateComboNoNoiseData2)(
    cell_lines, drugs, FALSE
  )
  expect_length(mae$warnings, 4)

  test_synthetic_data(original, mae$result, data)
})
