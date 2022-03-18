test_that(".calculate_matrix_metric works as expected", {
  n <- 10
  sa1 <- data.frame(conc = seq(n), conc2 = rep(0, n), metric = seq(n))
  sa2 <- data.frame(conc = rep(0, n), conc2 = seq(n), metric = seq(n))
  obs <- gDRcore:::.calculate_matrix_metric(sa1, series_id1 = "conc", sa2, series_id2 = "conc2", "metric", sum)
  expect_equal(dim(obs), c(100, 5))

  # Validates data.
  temp1 <- sa2
  temp2 <- sa1
  expect_error(gDRcore:::.calculate_matrix_metric(temp1, series_id1 = "conc", temp2, series_id2 = "conc2", "metric", sum))
})
