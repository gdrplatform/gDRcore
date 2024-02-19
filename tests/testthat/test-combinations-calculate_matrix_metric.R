test_that("calculate_HSA works as expected", {
  n <- 10
  sa1 <- data.table::data.table(conc = seq(n), conc2 = rep(0, n), smooth = seq(n))
  sa2 <- data.table::data.table(conc = rep(0, n), conc2 = seq(n), smooth = seq(n))
  hsa <- calculate_HSA(sa1, "conc", sa2, "conc2", "smooth")
  expect_equal(dim(hsa), c(100, 5))
})

test_that("calculate_Bliss works as expected", {
  n <- 10
  sa1 <- data.table::data.table(conc = seq(n), conc2 = rep(0, n), smooth = seq(n))
  sa2 <- data.table::data.table(conc = rep(0, n), conc2 = seq(n), smooth = seq(n))
  bliss <- calculate_Bliss(sa1, "conc", sa2, "conc2", "smooth")
  expect_equal(dim(bliss), c(100, 5))
})

test_that(".calculate_matrix_metric works as expected", {
  n <- 10
  sa1 <- data.table::data.table(conc = seq(n), conc2 = rep(0, n), smooth = seq(n))
  sa2 <- data.table::data.table(conc = rep(0, n), conc2 = seq(n), smooth = seq(n))
  obs <- gDRcore:::.calculate_matrix_metric(sa1, series_id1 = "conc", sa2, series_id2 = "conc2",
                                            "smooth", sum, measured_col = "smooth")
  expect_equal(dim(obs), c(n ^ 2, 5))

  # Validates data.
  temp1 <- sa2
  temp2 <- sa1
  expect_error(
    gDRcore:::.calculate_matrix_metric(temp1, series_id1 = "conc", temp2,
                                       series_id2 = "conc2", "smooth", sum, measured_col = "smooth")
  )
})
