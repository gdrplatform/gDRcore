test_that("calculate_excess works as expected", {
  metric <- data.frame(Concentration = c(1, 2, 3, 1, 2, 3),
    Concentration_2 = c(1, 1, 1, 2, 2, 2),
    GRvalue = c(100, 200, 300, 400, 500, 600))
  measured <- data.frame(Concentration = c(3, 1, 2, 2, 1, 3),
    Concentration_2 = c(1, 1, 1, 2, 2, 2),
    testvalue = c(200, 0, 100, 400, 300, 500))
  series_identifiers <- c("Concentration", "Concentration_2")
  metric_col <- "GRvalue"
  measured_col <- "testvalue"
  obs <- calculate_excess(metric, measured, series_identifiers, metric_col, measured_col)
  
  expect_equal(obs[, series_identifiers], measured[, series_identifiers])
  expect_equal(obs$excess, rep(100, 6))
})
