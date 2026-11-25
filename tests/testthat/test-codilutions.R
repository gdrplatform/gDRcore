test_that("fit_combo_codilutions works as expected", {
  n <- 8
  concs <- 10 * 3 ^ seq(0, -n, -1)
  start <- gDRutils::predict_efficacy_from_conc(concs, x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
  vals <- NULL
  for (i in seq(n)) {
    vals <- c(vals, start * start[i])
  }
  nested_identifiers <- c("Concentration", "Concentration_2")
  measured <- data.table::data.table(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    x = vals,
    normalization_type = "GR")

  obs <- fit_combo_codilutions(measured, nested_identifiers, "GR")
  expect_equal(dim(obs), c(7, 20))
  expect_true(all(c("ratio", "fit_message") %in% colnames(obs)))
  expect_equal(obs$ratio, c(.04, .1, .3, 1, 3, 10, 30), tolerance = 10e-3)
})


test_that("fit_codilution_series works as expected", {
  n <- 10
  concs <- seq(1, n, 1)
  start <- gDRutils::predict_efficacy_from_conc(concs, x_inf = 0.1, x_0 = 1, ec50 = 0.5, h = 2)
  vals <- NULL
  ratio <- 0.5
  for (i in seq(n)) {
    vals <- c(vals, start * start[i])
  }
  nested_identifiers <- c("Concentration", "Concentration_2")
  measured <- data.table::data.table(Concentration = rep(concs, n),
    Concentration_2 = rep(concs, each = n),
    x = vals,
    normalization_type = "GR")

  ratios <- measured$Concentration_2 / measured$Concentration
  keep <- !is.na(ratios) & ratios == ratio
  codilution <- measured[keep, ]

  res <- purrr::quietly(fit_codilution_series)(
    codilution,
    series_1 = "Concentration",
    series_2 = "Concentration_2",
    e_0 = 1,
    GR_0 = 1,
    normalization_type = "GR"
  )

  # the single gDR fit labels its outcome in fit_type instead of warning per curve
  expect_false(any(grepl("overriding original x_0", res$warnings)))

  obs <- res$result
  expect_equal(dim(obs), c(1, 20))
  expect_true(all(c("ratio", "fit_message") %in% colnames(obs)))
  expect_equal(obs$ratio, ratio, tolerance = 10e-3)
})
