test_that("combo_matrix: test_synthetic_data", {
  data <- "finalMAE_combo_matrix.RDS"
  original <- get_synthetic_data(data)
  
  mae <-
    gDRtestData::generateComboMatrix(cell_lines, drugs, e_inf, ec50, hill_coef)
  
  test_synthetic_data(original, mae, data)
})