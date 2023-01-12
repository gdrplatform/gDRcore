test_that("medium: test_synthetic_data", {
  data <- "finalMAE_medium.RDS"
  original <- get_synthetic_data(data)
  
  mae <- suppressWarnings({ # expected warnings
    gDRtestData::generateMediumData(cell_lines, drugs, e_inf, ec50, hill_coef)
  })
  
  test_synthetic_data(original, mae, data)
})
