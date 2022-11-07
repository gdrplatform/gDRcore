test_that("many_drugs: test_synthetic_data", {
  data <- "finalMAE_many_drugs.RDS"
  original <- get_synthetic_data(data)
  
  mae <-
    gDRtestData::generateManyDrugsData(cell_lines, drugs, e_inf, ec50, hill_coef)
  
  test_synthetic_data(original, mae, data)
})