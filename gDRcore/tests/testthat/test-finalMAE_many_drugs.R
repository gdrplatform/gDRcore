test_that("many_drugs: test_synthetic_data", {
  data <- "finalMAE_many_drugs.RDS"
  original <- get_synthetic_data(data)
  
  mae <-
    purrr::quietly(gDRtestData::generateManyDrugsData)(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
  expect_length(mae$warnings, 2)
  
  test_synthetic_data(original, mae$result, data)

})
