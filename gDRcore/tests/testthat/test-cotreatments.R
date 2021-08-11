test_that("fit_combo_cotreatments works as expected", {
  n <- 10
  concs <- seq(0, n, 1)
  start <- gDRutils::logistic_4parameters(concs, Vinf = 0.1, V0 = 1, EC50 = 0.5, h = 2)
  vals  <- NULL
  for (i in seq(n)) {
    vals <- c(vals, start + i*concs)
  }
  nested_identifiers <- c("Concentration", "Concentration_2")
  measured <- DataFrame(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    GRvalue = vals)

  obs <- gDRcore:::fit_combo_cotreatments(measured, nested_identifiers, "GR")
  expect_equal(dim(obs), c(20, 18))
  expect_true(nested_identifiers %in% colnames(obs))
  expect_equal(unname(vapply(obs[, nested_identifiers], class, "")), c("SimpleNumericList", "SimpleNumericList"))
})


test_that("fit_cotreatment_series works as expected", {
  n <- 10
  concs <- seq(1, n, 1)
  start <- gDRutils::logistic_4parameters(concs, Vinf = 0.1, V0 = 1, EC50 = 0.5, h = 2)
  vals  <- NULL
  for (i in seq(n)) {
    vals <- c(vals, start + i*concs)
  }
  nested_identifiers <- c("Concentration", "Concentration_2")
  measured <- DataFrame(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    GRvalue = vals)

  obs <- gDRcore:::fit_cotreatment_series(measured, 
    series_id = "Concentration",
    cotrt_id = "Concentration_2",
    e_0 = 1,
    GR_0 = 1,
    cotrt_value = 1,
    normalization_type = "GR")

  expect_equal(dim(obs), c(1, 18))
  expect_true(all(nested_identifiers %in% colnames(obs)))
  expect_equal(unname(vapply(obs[, nested_identifiers], class, "")), c("SimpleNumericList", "SimpleNumericList"))
})
