.round_params <- function(df) {
  df[] <- lapply(df, round, 4)
  df
}


test_that("fit_curves fails with expected errors", {
  expect_error(fit_curves(list()))

  # Log10 concentrations.
  df_resp2 <- df_resp
  df_resp2$Concentration <- df_resp2$Concentration * -1
  expect_error(fit_curves(df_resp2, series_identifiers = "Concentration"),
    reg = "logisticFit accepts only unlogged concentrations")

  # Invalid normalization_type.
  expect_error(fit_curves(df_resp, series_identifiers = "Concentration", normalization_type = "BOGUS"),
    reg = "unknown curve type")
  expect_error(fit_curves(df_resp, series_identifiers = "Concentration", normalization_type = c("GRvalue", "BOGUS")),
    reg = "unknown curve type")

  expect_error(fit_curves(
    df_resp,
    series_identifiers = c("Concentration", "add_conc"),
    normalization_type = c("GRvalue")
  ),
  reg = "gDR does not yet support multiple series_identifiers, feature coming soon")
})

test_that("NA values are handled correctly", {
  df_resp_NA <- df_resp
  df_resp_NA[, "x"] <- NA
  expect_warning(fit_curves(df_resp_NA, series_identifiers = "Concentration"))
  
  df_result_NA <- purrr::quietly(fit_curves)(df_resp_NA, series_identifiers = "Concentration")
  expect_length(df_result_NA$warnings, 2)
  expect_true(all(is.na(df_result_NA$result[, "xc50"])))
})


test_that("appropriate fit type is assigned for various use cases", {
  set.seed(1112020) # For reproducibility.

  # Test a 3P fit.
  ## Note that this should correspond to a cytotoxic response.
  df_result <- fit_curves(df_resp, series_identifiers = "Concentration")
  expect_equal(round(df_result[, names(params), with = FALSE], 4), expected, tolerance = 1e-5)

  obs_fit <- unname(unlist(unique(df_result[, "fit_type"])))
  expect_equal(obs_fit, "DRC3pHillFitModelFixS0")
  expect_equal(dim(df_result), expected_dims)

  # Test a 4P fit (without the x_0 value).
  df_result <- fit_curves(df_resp, series_identifiers = "Concentration", e_0 = NA, GR_0 = NA)
  expect_equal(round(df_result[, names(params), with = FALSE], 4), expected, tolerance = 1e-5)
  obs_fit <- unname(unlist(unique(df_result[, "fit_type"])))
  expect_equal(obs_fit, "DRC4pHillFitModel")
  expect_equal(dim(df_result), expected_dims)

  # Test for constant fit.
  df_resp4 <- df_resp
  df_resp4[df_resp4$normalization_types == "RV", "x"] <-
    df_resp4[df_resp4$normalization_types == "GR", "x"] <- 0.5

  expect_warning(df_result <- fit_curves(df_resp4, series_identifiers = "Concentration"),
    reg = "overriding original x_0 argument") # Override.
  obs_fit <- unname(unlist(unique(df_result[, "fit_type"])))
  expect_equal(obs_fit, "DRCConstantFitResult")
  expect_equal(unname(unlist(df_result[1, c("x_0", "x_inf", "x_mean", "x_AOC", "x_AOC_range")])),
    rep(unname(unlist(unique(df_resp4[df_resp4$normalization_types == "RV", "x"]))), 5))
  expect_equal(dim(df_result), expected_dims)

  ## Test for all values below 0.5.
  df_resp5 <- df_resp

  # Scale all readouts.
  max_rv <- max(df_resp[df_resp$normalization_types == "RV", "x"])
  min_rv <- min(df_resp[df_resp$normalization_types == "RV", "x"])

  max_gr <- max(df_resp[df_resp$normalization_types == "GR", "x"])
  min_gr <- min(df_resp[df_resp$normalization_types == "GR", "x"])

  df_resp5[df_resp5$normalization_types == "RV", "x"] <-
    (df_resp[df_resp$normalization_types == "RV", "x"] / (max_rv - min_rv)) *
    ((max_rv - min_rv) / 2)
  df_resp5[df_resp5$normalization_types == "GR", "x"] <-
    (df_resp[df_resp$normalization_types == "GR", "x"] / (max_gr - min_gr)) *
    ((max_gr - min_gr) / 2)
  df_result5 <- fit_curves(df_resp5, series_identifiers = "Concentration")
  expect_equal(unname(unlist(df_result5[, c("x_inf")])), c(0, -1))
  obs_fit <- unname(unlist(unique(df_result5[, "fit_type"])))
  expect_equal(obs_fit, "DRC3pHillFitModelFixS0")

  # Test for all values above 0.5.
  ## Note that this corresponds to partial growth inhibition.
  df_result6 <- fit_curves(df_resp_above, series_identifiers = "Concentration")
  expect_equal(unname(unlist(df_result6[, c("x_inf")])), c(0.5, 0.5), tolerance = 1e-5)
  expect_true(all(df_result6[, c("xc50")] > 500))
  obs_fit <- unname(unlist(unique(df_result6[, "fit_type"])))
  expect_equal(obs_fit, "DRC3pHillFitModelFixS0")

  # Test for a pushed constant fit by adding noise.
  df_resp7 <- df_resp_above
  noise <- sample(seq(-1, 1, 0.1), nrow(df_resp7) / 2)
  emax <- 0.8
  df_resp7[df_resp7$normalization_types == "RV", "x"] <-
    pmin(unname(unlist(df_resp7[df_resp7$normalization_types == "RV", "x"])) + noise, emax)
  df_resp7[df_resp7$normalization_types == "GR", "x"] <-
    pmin(unname(unlist(df_resp7[df_resp7$normalization_types == "GR", "x"])) + noise, emax)

  expect_warning(
    df_result7 <- fit_curves(df_resp7, series_identifiers = "Concentration", force_fit = FALSE),
    ".*overriding original x_0 argument.*"
  )
  expect_equal(unname(unlist(unique(df_result7[, "fit_type"]))), "DRCConstantFitResult")
  expect_equal(unique(unname(unlist(df_result7[, c("x_mean", "x_inf", "x_0")]))),
    0.70781, tolerance = 1e-5)

  # Test that force argument overrides as expected.
  df_result8 <- fit_curves(df_resp7, series_identifiers = "Concentration", force_fit = TRUE)
  expect_equal(unname(unlist(unique(df_result8[, "fit_type"]))), "DRC3pHillFitModelFixS0")

  # Test that pcutoff argument works as expected.
  expect_warning(
    df_result9 <- fit_curves(df_resp7, series_identifiers = "Concentration", force_fit = FALSE, pcutoff = 1),
    ".*overriding original x_0 argument.*"
  )

  expect_equal(unname(unlist(df_result9[, "fit_type"])), c("DRC3pHillFitModelFixS0", "DRCConstantFitResult"))
  df_result10 <- fit_curves(df_resp7, series_identifiers = "Concentration",
    force_fit = FALSE, pcutoff = 1.01) # Essentially equivalent to a 'force = TRUE'.
  expect_equal(unname(unlist(unique(df_result10[, "fit_type"]))), "DRC3pHillFitModelFixS0")

  # Test for GR values from 0-1.
  ## Note that this correspond to a fully cytostatic response (no growth).
  df_resp11 <- df_resp
  df_resp11[df_resp11$normalization_types == "GR", "x"] <-
    df_resp11[df_resp11$normalization_types == "RV", "x"]

  df_result11 <- fit_curves(df_resp11, series_identifiers = "Concentration")
  expect_equal(unname(unlist(unique(df_result11[, "x_inf"]))), c(0.1, 0.1), tolerance = 1e-5)

  # Test for too few points.
  df_result <- purrr::quietly(fit_curves)(df_resp[c(3:5, 12:14), ],
                          series_identifiers = "Concentration", n_point_cutoff = 4)
  expect_length(df_result$warnings, 2)
  obs_fit <- unique(unlist(df_result$result[, "fit_type"]))
  expect_equal(obs_fit, "DRCTooFewPointsToFit")
  expect_equal(dim(df_result$result), expected_dims)

  #nolint start
    # TODO: Test for invalid fit. Maybe try a bunch of noise.
    #  expect_warning(df_result <- fit_curves(df_resp[3:5, ], series_identifiers = "Concentration", n_point_cutoff = 1),
      #reg = "fitting failed with error")
    #
    #  obs_fit <- unique(df_result[, "fit_type"])
    #  expect_equal(obs_fit, "DRCInvalidFitResult")
    #  expect_equal(dim(df_result), expected_dims)
  #nolint end
})


test_that("normalization_type can be specified", {
  GR_df_result <- fit_curves(df_resp, series_identifiers = "Concentration", normalization_type = "GR")
  expect_equal(rownames(GR_df_result), "GR_gDR")
  expect_equal(round(GR_df_result[, names(params), with = FALSE], 4), expected[2, ], tolerance = 1e-5)

  RV_df_result <- fit_curves(df_resp, series_identifiers = "Concentration", normalization_type = "RV")
  expect_equal(rownames(RV_df_result), "RV_gDR")
  expect_equal(round(RV_df_result[, names(params), with = FALSE], 4), expected[1, ], tolerance = 1e-5)
})

test_that("predict_efficacy_from_conc works as expected", {
  c <- 1
  Vinf <- 0.1
  V0 <- 1
  h <- 2
  EC50 <- 0.5

  # Non - numeric values cause an error.
  expect_error(
    predict_efficacy_from_conc(
      c = "non-numeric_entry",
      x_inf = Vinf,
      x_0 = V0,
      ec50 = EC50,
      h = h
    )
  )

  # Normal fit.
  v <- predict_efficacy_from_conc(
    c = c,
    x_inf = Vinf,
    x_0 = V0,
    ec50 = EC50,
    h = h
  )
  expect_equal(v, 0.28)

  # Flat fit.
  EC50 <- 0
  v <- predict_efficacy_from_conc(
    c = c,
    x_inf = Vinf,
    x_0 = V0,
    ec50 = EC50,
    h = h
  )
  expect_equal(v, Vinf)

  # Multiple concentrations.
  conc <- c(1, 1.5)
  obs <- predict_conc_from_efficacy(conc, Vinf, V0, EC50, h)
  expect_equal(length(obs), length(conc))

})


test_that("predict_conc_from_efficacy works as expected", {
  c <- 1
  Vinf <- 0.1
  V0 <- 1
  h <- 2
  EC50 <- 0.5

  efficacy <- predict_efficacy_from_conc(c, Vinf, V0, EC50, h)

  # Normal in-range.
  expect_equal(predict_conc_from_efficacy(efficacy, Vinf, V0, EC50, h), c)

  # Edge case: efficacy > x_0.
  V0 <- 0.8
  efficacy <- 0.9
  expect_equal(predict_conc_from_efficacy(efficacy, Vinf, V0, EC50, h), 0)

  # Edge case: efficacy < x_inf.
  V0 <- 1
  efficacy <- 0.05
  expect_equal(predict_conc_from_efficacy(efficacy, Vinf, V0, EC50, h), Inf)

  # Multiple efficacies.
  efficacy <- c(0.05, 0.15)
  obs <- predict_conc_from_efficacy(efficacy, Vinf, V0, EC50, h)
  expect_equal(length(obs), length(efficacy))
  expect_equal(obs[1], Inf)
})


###################
# Helper functions
###################

test_that(".setup_metric_output works as expected", {
  obs <- .setup_metric_output()
  expect_true(is.list(obs))
  expect_equal(length(obs), 16)
})


test_that(".estimate_xc50 works as expected", {
  expect_equal(.estimate_xc50(c(NA, NA, NA, NA)), NA)
  expect_equal(.estimate_xc50(c(0.6, 0.7, 0.8, 0.9)), Inf)
  expect_equal(.estimate_xc50(c(0.1, 0.2, 0.3, 0.4)), -Inf)
  expect_equal(.estimate_xc50(c(0.1, 0.2, 0.6, 0.7)), NA)
})


test_that("average_dups works as expected", {
  df <- data.table::data.table(concs = rep(seq(5), each = 2),
                               norm_value = seq(10))
  expect_equal(average_dups(df, "concs"), 
               data.table::data.table(concs = seq(5), norm_value = seq(1.5, 9.5, 2)))
})


###########
# Setters
###########

test_that(".set_mean_params works as expected", {
  out <- list()

  # Above 0.5.
  v <- 0.6
  above <- .set_mean_params(out, mean_norm_value = v)
  expect_true(all(c(above$x_0, above$x_inf, above$x_mean) == v))
  expect_true(all(c(above$x_AOC, above$x_AOC_range) == (1 - v)))
  expect_equal(above$xc50, Inf)

  # Below 0.5.
  v <- 0.4
  below <- .set_mean_params(out, mean_norm_value = v)
  expect_true(all(c(below$x_0, below$x_inf, below$x_mean) == v))
  expect_equal(below$xc50, -Inf)

  v <- seq(3)
  expect_error(.set_mean_params(out, mean_norm_value = v))
})


test_that("set_constant_fit_params works as expected", {
  x_0 <- NA
  na <- list(x_0 = x_0)
  na <- set_constant_fit_params(na, mean_norm_value = 0.6)

  expect_equal(na$ec50, 0)
  expect_equal(na$r2, 0)
  expect_equal(na$h, 0.0001)
  expect_equal(na$fit_type, "DRCConstantFitResult")
  expect_true(all(c(na$x_0, na$x_inf, na$x_mean) == 0.6))
  expect_true(all(c(na$x_AOC, na$x_AOC_range) == (1 - 0.6)))
  expect_equal(na$xc50, Inf)

  x_0 <- 1
  out <- list(x_0 = x_0)
  one <- set_constant_fit_params(out, mean_norm_value = 0.6)

  expect_equal(one$ec50, 0)
  expect_equal(na$r2, 0)
  expect_equal(one$h, 0.0001)
  expect_equal(one$fit_type, "DRCConstantFitResult")
  expect_true(all(c(one$x_0, one$x_inf, one$x_mean) == 0.6))
  expect_true(all(c(one$x_AOC, one$x_AOC_range) == (1 - 0.6)))
  expect_equal(one$xc50, Inf)
})

test_that(".set_invalid_fit_params works as expected", {
  out <- list()

  # All NA.
  obs_na <- .set_invalid_fit_params(out, norm_values = rep(NA, 6))
  expect_equal(obs_na$xc50, NA)

  # > 0.5.
  obs_high <- .set_invalid_fit_params(out, norm_values = rep(0.7, 6))
  expect_equal(obs_high$xc50, Inf)

  # < 0.5.
  obs_low <- .set_invalid_fit_params(out, norm_values = rep(0.3, 6))
  expect_equal(obs_low$xc50, -Inf)

  # Mixed.
  obs_mixed <- .set_invalid_fit_params(out, norm_values = rep(c(0.3, 0.7), 3))
  expect_equal(obs_mixed$xc50, NA)
  expect_equal(obs_mixed$r2, NA)
})

###############
# Calculations
###############

test_that(".calculate_x_max works as expected", {
  n <- 4
  df <- data.frame(concs = c(0.01, 1, 0.03, 5), norm_values = seq(n))
  expect_equal(.calculate_x_max(df), 2)

  df$norm_values <- rep(NA, n)
  expect_equal(.calculate_x_max(df), NA)
})


test_that(".calculate_xc50 works as expected", {
  Vinf <- 0.1
  V0 <- 1
  h <- 2
  EC50 <- 0.5
  obs <- .calculate_xc50(ec50 = EC50, x0 = V0, xInf = Vinf, h = h)
  expect_equal(obs, 0.559017)
})


test_that("cap_xc50 works as expected", {
  expect_error(cap_xc50(xc50 = c(1, 2), max_conc = c(10, 10), capping_fold = 5))
  
  expect_equal(cap_xc50(xc50 = 26, max_conc = 5, capping_fold = 5), Inf)
  expect_equal(cap_xc50(xc50 = 1e-6, max_conc = 5, capping_fold = 5), -Inf)
  expect_equal(cap_xc50(xc50 = 1, max_conc = 5, capping_fold = 5), 1)
})

test_that("predict_efficacy_from_conc works as expected", {
  # params
  h <- 2
  x_inf <- 0.1
  x_0 <- 1
  ec50 <- 0.5
  conc <- c(0, 10 ^ (seq(-3, 1, 0.5)))
  
  out <- predict_efficacy_from_conc(conc, x_inf, x_0, ec50, h)
  
  res <- c(x_0, 
           vapply(conc[2:NROW(conc)], 
                  function(c) x_inf + (x_0 - x_inf) / (1 + (c / ec50) ^ h), numeric(1)))
  
  expect_equal(out, res)
})
