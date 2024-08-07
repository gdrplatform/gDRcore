test_that("fit_combo_cotreatments works as expected", {
  n <- 10
  concs <- seq(0, n, 1)
  start <- gDRutils::predict_efficacy_from_conc(concs, x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
  vals  <- NULL
  for (i in seq(n)) {
    vals <- c(vals, start + i * concs)
  }
  measured <- data.table::data.table(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    GRvalue = vals)

  obs <- suppressWarnings(# warnings are expected
    fit_combo_cotreatments(
      measured, 
      series_id = "Concentration",
      cotrt_id = "Concentration_2",
      normalization_type = "GR"
  ))
  expect_equal(dim(obs), c(11, 19))
  expect_true("cotrt_value" %in% colnames(obs))
  expect_equal(class(obs$cotrt_value), "numeric")
  # TODO: test that all the fits are the same.
})


test_that("fit_cotreatment_series works as expected", {
  n <- 10
  concs <- seq(1, n, 1)
  start <- gDRutils::predict_efficacy_from_conc(concs, x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
  vals  <- NULL
  for (i in seq(n)) {
    vals <- c(vals, start + i * concs)
  }
  nested_identifiers <- c("Concentration", "Concentration_2")
  measured <- data.table::data.table(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    x = vals,
    normalization_type = "GR")

  obs <- suppressWarnings(# warnings are expected
    fit_cotreatment_series(
      measured, 
      series_id = "Concentration",
      cotrt_id = "Concentration_2",
      e_0 = 1,
      GR_0 = 1,
      cotrt_value = 1,
      normalization_type = "GR"
  ))

  expect_equal(dim(obs), c(1, 19))
  expect_true("cotrt_value" %in% colnames(obs))
  expect_equal(class(obs$cotrt_value), "numeric")
})
