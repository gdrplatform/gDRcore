test_that("wLigand: test_synthetic_data", {
  data <- "finalMAE_wLigand.RDS"
  original <- get_synthetic_data(data)
  
  mae <-
    gDRtestData::generateLigandData(cell_lines, drugs, e_inf, ec50, hill_coef)
  
  test_synthetic_data(original, mae, data)
})