# Python with the synergy package is not available in CI, so the tests that fit
# a surface are skipped there; everything before the Python call is covered.

.synergy_available <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_module_available("synergy")
}

# Combination grid in the shape apply_fit() hands to a fit function
.musyc_grid_dt <- function(norm_type = "GR", conc = c(0, 0.01, 0.1, 1, 10), x = NULL) {
  dt <- data.table::CJ(Concentration = conc, Concentration_2 = conc)
  dt <- dt[!(Concentration == 0 & Concentration_2 == 0)]
  dt[, normalization_type := norm_type]
  dt[, x := if (is.null(x)) seq(1, 0.2, length.out = .N) else x]
  dt
}

# Grid sampled from a MuSyC surface with known parameters
.musyc_truth_dt <- function(params, norm_type = "GR") {
  synergy <- reticulate::import("synergy.combination")
  truth <- synergy$MuSyC(fit_gamma = FALSE)
  for (name in names(params)) {
    truth[[name]] <- params[[name]]
  }
  truth$gamma12 <- 1
  truth$gamma21 <- 1
  conc <- c(0, 10^seq(-3, 1, length.out = 8))
  dt <- .musyc_grid_dt(norm_type, conc)
  dt[, x := as.numeric(truth$E(
    reticulate::np_array(Concentration), reticulate::np_array(Concentration_2)
  ))]
  dt
}

.fake_musyc_rows <- function(beta = c(0.4, 0.2), converged = c(TRUE, TRUE)) {
  data.table::data.table(
    normalization_type = c("GR", "RV")[seq_along(beta)],
    alpha12 = 2, alpha21 = 1.5, beta = beta, r2 = 0.9, converged = converged
  )
}


test_that("musyc_fit_fn validates its arguments", {
  dt <- .musyc_grid_dt()
  expect_error(musyc_fit_fn(as.list(dt)), "avg_dt")
  expect_error(musyc_fit_fn(dt, fit_gamma = "yes"), "fit_gamma")
  expect_error(musyc_fit_fn(dt, n_restarts = -1), "n_restarts")
  expect_error(musyc_fit_fn(dt, seed = 1.5), "seed")
})

test_that("musyc_fit_fn returns an NA result for a grid it cannot fit", {
  dt <- data.table::data.table(
    normalization_type = "GR",
    Concentration = c(0.1, 1),
    Concentration_2 = c(0, 0),
    x = c(0.8, 0.5)
  )
  result <- expect_warning(musyc_fit_fn(dt), "insufficient data")
  expect_false(result$converged)
  expect_true(is.na(result$alpha12))
  expect_identical(names(result), names(.musyc_na_result("GR")))
})

test_that("musyc_fit_fn stops with an install hint when synergy is missing", {
  testthat::skip_if(.synergy_available())
  expect_error(musyc_fit_fn(.musyc_grid_dt()), "synergy")
})

test_that(".musyc_grid adds the untreated corner and drops missing responses", {
  obs <- .musyc_grid(c(0, 0.1, 1), c(0.1, 0, 1), c(0.5, NA, 0.2))
  expect_equal(NROW(obs), 3L)
  expect_equal(obs$E[obs$d1 == 0 & obs$d2 == 0], 1)
  expect_false(anyNA(obs$E))
})

test_that(".musyc_bounds lets E go negative for GR only", {
  expect_equal(.musyc_bounds(-1)$E3_bounds, c(-1, 1.5))
  expect_equal(.musyc_bounds(0)$E3_bounds, c(0, 1.5))
  expect_equal(.musyc_bounds(0)$E0_bounds, c(0.99, 1.01))
})

test_that(".with_local_seed leaves the global random number generator untouched", {
  set.seed(7)
  before <- .Random.seed
  first <- .with_local_seed(1L, stats::runif(3))
  expect_identical(.Random.seed, before)
  expect_identical(.with_local_seed(1L, stats::runif(3)), first)
})

test_that("musyc_fit_fn recovers the parameters of a known surface", {
  testthat::skip_if_not(.synergy_available(), "Python synergy package not available")
  params <- list(
    E0 = 1, E1 = 0.4, E2 = 0.3, E3 = -0.2, h1 = 1.2, h2 = 0.9,
    C1 = 0.5, C2 = 0.2, alpha12 = 3, alpha21 = 0.5
  )
  dt <- .musyc_truth_dt(params)
  result <- musyc_fit_fn(dt)

  expect_true(result$converged)
  expect_gt(result$r2, 0.999)
  expect_equal(result$E3, params$E3, tolerance = 0.05)
  expect_equal(log10(result$alpha12), log10(params$alpha12), tolerance = 0.1)
  expect_identical(musyc_fit_fn(dt), result)
})

test_that("musyc_fit_fn fits through apply_fit and records the Python runtime", {
  testthat::skip_if_not(.synergy_available(), "Python synergy package not available")
  mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
  se_out <- gDRutils::apply_fit(
    combo_se, musyc_fit_fn, "combination",
    output_assay = "musyc_params",
    summary_fn = musyc_summary_fn, summary_assay = "musyc_summary",
    fit_source = "musyc"
  )

  params <- gDRutils::convert_se_assay_to_dt(se_out, "musyc_params")
  expect_true(all(c("alpha12", "beta", "r2", "converged") %in% names(params)))
  expect_true(any(params$converged))
  expect_true("musyc_summary" %in% SummarizedExperiment::assayNames(se_out))
  runtime <- unlist(S4Vectors::metadata(se_out)$provenance, use.names = TRUE)
  expect_true(any(grepl("synergy", names(runtime))))
})

test_that("musyc_summary_fn calls synergy from the mean beta", {
  expect_equal(musyc_summary_fn(.fake_musyc_rows(beta = c(0.5, 0.3)))$synergy_call, "synergistic")
  expect_equal(musyc_summary_fn(.fake_musyc_rows(beta = c(-0.5, -0.3)))$synergy_call, "antagonistic")
  expect_equal(musyc_summary_fn(.fake_musyc_rows(beta = c(0.05, 0)))$synergy_call, "additive")
  expect_equal(musyc_summary_fn(.fake_musyc_rows(converged = c(FALSE, FALSE)))$synergy_call, "unknown")
})

test_that("musyc_summary_fn averages converged fits only", {
  result <- musyc_summary_fn(.fake_musyc_rows(beta = c(0.6, 0.2), converged = c(TRUE, FALSE)))
  expect_equal(result$mean_beta, 0.6)
  expect_true(result$any_converged)
  expect_gt(result$mean_log2_alpha, 0)
  expect_error(musyc_summary_fn(list(converged = TRUE)), "fit_dt")
})
