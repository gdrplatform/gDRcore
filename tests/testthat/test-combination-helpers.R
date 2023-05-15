test_that("map_ids_to_fits works as expected", {
  pred <- c(1, 5, 5)
  match_col <- c(1, 1, 2)
  fitting_id_col <- "match_on_me"
  
  fit1 <- data.table::data.table(h = 2.09, x_inf = 0.68, x_0 = 1, ec50 = 0.003)
  fit2 <- data.table::data.table(h = 0.906, x_inf = 0.46, x_0 = 1, ec50 = 0.001)
  fittings <- do.call(rbind, list(fit1, fit2))
  fittings[[fitting_id_col]] <- c(1, 2)
  
  obs <- map_ids_to_fits(pred, match_col, fittings, fitting_id_col)
  exp1 <- gDRutils::predict_efficacy_from_conc(pred[1], fit1$x_inf, fit1$x_0, fit1$ec50, fit1$h)
  exp2 <- gDRutils::predict_efficacy_from_conc(pred[2], fit1$x_inf, fit1$x_0, fit1$ec50, fit1$h)
  exp3 <- gDRutils::predict_efficacy_from_conc(pred[3], fit2$x_inf, fit2$x_0, fit2$ec50, fit2$h)
  
  expect_equal(obs, c(exp1, exp2, exp3))
})

test_that("calculate_excess works as expected", {
  metric <- data.table::data.table(Concentration = c(1, 2, 3, 1, 2, 3),
                       Concentration_2 = c(1, 1, 1, 2, 2, 2),
                       GRvalue = c(100, 200, 300, 400, 500, 600))
  measured <- data.table::data.table(Concentration = c(3, 1, 2, 2, 1, 3),
                         Concentration_2 = c(1, 1, 1, 2, 2, 2),
                         testvalue = c(200, 0, 100, 400, 300, 500))
  series_identifiers <- c("Concentration", "Concentration_2")
  metric_col <- "GRvalue"
  measured_col <- "testvalue"
  obs <- calculate_excess(metric, measured, series_identifiers, metric_col, measured_col)
  
  expect_equal(obs[, ..series_identifiers], measured[, ..series_identifiers])
  expect_equal(obs$x, rep(100, 6))
})



test_that(".calculate_dilution_ratio works as expected", {
  ratio <- 0.5
  concs <- 10 ^ (seq(-3, 1, ratio))
  obs <- gDRcore:::.calculate_dilution_ratio(concs)
  expect_equal(obs, ratio)
})


test_that("map_conc_to_standardized_conc works as expected", {
  ratio <- 0.5
  conc1 <- c(0, 10 ^ (seq(-3, 1, ratio)))

  shorter_range <- conc1[-1]
  noise <- runif(length(shorter_range), 1e-12, 1e-11)
  conc2 <- shorter_range + noise

  obs <- map_conc_to_standardized_conc(conc1, conc2)
  expect_true(methods::is(obs, "data.table"))
})


test_that("replace_conc_with_standardized_conc works as expected", {
  conc_map <- data.table::data.table(orig = c(0.99, 0.6, 0.456, 0.4), std = c(1, 0.6, 0.46, 0.4))
  original_concs <- c(0.456, 0.456, 0.4, 0.99)
  exp <- c(0.46, 0.46, 0.4, 1)
  obs <- replace_conc_with_standardized_conc(original_concs, conc_map,
                                             original_conc_col = "orig", standardized_conc_col = "std")
  expect_equal(unname(obs), exp)
})

test_that(".standardize_conc works as expected", {
  concs <- 10 ^ (seq(-1, 1, 0.9))
  obs <- .standardize_conc(concs)
  expect_equal(obs, c(0.1, 0.794, 6.31))
})
