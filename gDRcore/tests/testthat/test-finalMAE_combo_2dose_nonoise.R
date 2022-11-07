test_that("combo_2dose_nonoise: test_synthetic_data", {
  data <- "finalMAE_combo_2dose_nonoise.RDS"
  original <- get_synthetic_data(data)
  
  mae <-
    gDRtestData::generateComboNoNoiseData(cell_lines, drugs, e_inf, ec50, hill_coef)
  
  test_synthetic_data(original, mae, data)
})