test_that("combo_2dose_nonoise3: test_synthetic_data", {
  data <- "finalMAE_combo_2dose_nonoise3.RDS"
  original <- get_synthetic_data(data)
  
  mae <- suppressWarnings(# expected warnings
    gDRtestData::generateComboNoNoiseData3(cell_lines, drugs, e_inf, ec50, hill_coef)
  )

  test_synthetic_data(original, mae, data)
})
