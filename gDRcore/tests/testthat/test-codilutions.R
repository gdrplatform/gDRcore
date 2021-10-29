library(gDRcore)
library(testthat)

test_that("fit_combo_codilutions works as expected", {
  n <- 8
  concs <- 10 * 3 ^ seq(0, -n, -1)
  start <- gDRutils::logistic_4parameters(concs, Vinf = 0.1, V0 = 1, EC50 = 0.5, h = 2)
  vals  <- NULL
  for (i in seq(n)) {
    vals <- c(vals, start * start[i])
  }
  nested_identifiers <- c("Concentration", "Concentration_2")
  measured <- DataFrame(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    GRvalue = vals)

  obs <- gDRcore:::fit_combo_codilutions(measured, nested_identifiers, "GRvalue")
  expect_equal(dim(obs), c(7, 17))
  expect_true("ratio" %in% colnames(obs))
  expect_equal(obs$ratio, c(.04, .1, .3, 1, 3, 10, 30), tolerance = 10e-3)
})


test_that("fit_codilution_series works as expected", {
  n <- 10
  concs <- seq(1, n, 1)
  start <- gDRutils::logistic_4parameters(concs, Vinf = 0.1, V0 = 1, EC50 = 0.5, h = 2)
  vals  <- NULL
  ratio <- 0.5
  for (i in seq(n)) {
    vals <- c(vals, start * start[i])
  }
  nested_identifiers <- c("Concentration", "Concentration_2")
  measured <- DataFrame(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    GRvalue = vals)

  ratios <- measured$Concentration_2 / measured$Concentration
  keep <- !is.na(ratios) & ratios == ratio
  codilution <- measured[keep, ]

  obs <- gDRcore:::fit_codilution_series(codilution, 
    series_1 = "Concentration",
    series_2 = "Concentration_2",
    e_0 = 1,
    GR_0 = 1,
    normalization_type = "GRvalue")

  expect_equal(dim(obs), c(1, 17))
  expect_true("ratio" %in% colnames(obs))
  expect_equal(obs$ratio, ratio, tolerance = 10e-3)
})
