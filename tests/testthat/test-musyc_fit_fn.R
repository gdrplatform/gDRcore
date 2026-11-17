# Combination grid in the shape apply_fit() hands to a fit function
.musyc_grid_dt <- function(norm_type = "GR", conc = c(0, 0.01, 0.1, 1, 10), x = NULL) {
  dt <- data.table::CJ(Concentration = conc, Concentration_2 = conc)
  dt <- dt[!(Concentration == 0 & Concentration_2 == 0)]
  dt[, normalization_type := norm_type]
  dt[, x := if (is.null(x)) seq(1, 0.2, length.out = .N) else x]
  dt
}

.musyc_theta <- function(params) {
  p <- unlist(params[.musyc_param_names])
  c(p[1:4], log(p[5:10]))
}

# Grid sampled from a MuSyC surface with known parameters
.musyc_truth_dt <- function(params, norm_type = "GR", conc = c(0, 10^seq(-3, 1, length.out = 8))) {
  dt <- .musyc_grid_dt(norm_type, conc)
  dt[, x := .musyc_surface(.musyc_theta(params), Concentration, Concentration_2)]
  dt
}

.musyc_known_params <- list(
  E0 = 1, E1 = 0.4, E2 = 0.3, E3 = -0.2, h1 = 1.2, h2 = 0.9,
  C1 = 0.5, C2 = 0.2, alpha12 = 3, alpha21 = 0.5
)

.fake_musyc_rows <- function(beta = c(0.4, 0.2), converged = c(TRUE, TRUE)) {
  data.table::data.table(
    normalization_type = c("GR", "RV")[seq_along(beta)],
    alpha12 = 2, alpha21 = 1.5, beta = beta, r2 = 0.9, converged = converged
  )
}


test_that("musyc_fit_fn validates its arguments", {
  dt <- .musyc_grid_dt()
  expect_error(musyc_fit_fn(as.list(dt)), "avg_dt")
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

test_that("musyc_available is deprecated and always TRUE", {
  expect_true(expect_warning(musyc_available(), class = "deprecatedWarning"))
})

test_that(".musyc_surface matches the MuSyC model of the synergy package", {
  # synergy 1.0.0, MuSyC(fit_gamma = FALSE)$E() with .musyc_known_params
  d1 <- c(0, 0.1, 0, 1, 10, 0.05)
  d2 <- c(0, 0, 0.1, 1, 0.01, 5)
  expected <- c(1, 0.9240376357249616, 0.7557627758278932,
                0.10644444828648601, 0.3272282638114358, 0.321605942435626)
  expect_equal(.musyc_surface(.musyc_theta(.musyc_known_params), d1, d2), expected, tolerance = 1e-12)
})

test_that("predict_musyc_surface evaluates a fit at the given concentrations", {
  d1 <- c(0, 0.1, 0, 1, 10, 0.05)
  d2 <- c(0, 0, 0.1, 1, 0.01, 5)
  expected <- .musyc_surface(.musyc_theta(.musyc_known_params), d1, d2)
  expect_equal(predict_musyc_surface(d1, d2, .musyc_known_params), expected)
  expect_equal(predict_musyc_surface(d1, d2, data.table::as.data.table(.musyc_known_params)), expected)
  expect_error(predict_musyc_surface(d1, d2[-1], .musyc_known_params), "conc_2")
  expect_error(predict_musyc_surface(d1, d2, .musyc_known_params[-10]), "alpha21")
})

test_that(".musyc_surface returns the Jacobian of the surface", {
  d1 <- rep(c(0, 10^seq(-3, 1, length.out = 6)), each = 7)
  d2 <- rep(c(0, 10^seq(-3, 1, length.out = 6)), times = 7)
  theta <- .musyc_theta(.musyc_known_params)
  numeric_jacobian <- vapply(seq_along(theta), function(k) {
    step <- replace(numeric(length(theta)), k, 1e-6)
    (.musyc_surface(theta + step, d1, d2) - .musyc_surface(theta - step, d1, d2)) / 2e-6
  }, numeric(length(d1)))
  expect_equal(unname(.musyc_surface(theta, d1, d2, jacobian = TRUE)), numeric_jacobian,
               tolerance = 1e-7)
})

test_that(".musyc_surface returns a finite Jacobian for a very steep curve", {
  d1 <- rep(c(0, 10^seq(-3, 1, length.out = 6)), each = 7)
  d2 <- rep(c(0, 10^seq(-3, 1, length.out = 6)), times = 7)
  # with h1 = 400 the dose ratio powers reach exp(200), so the squared denominator overflows
  theta <- .musyc_theta(replace(.musyc_known_params, "h1", 400))
  numeric_jacobian <- vapply(seq_along(theta), function(k) {
    step <- replace(numeric(length(theta)), k, 1e-6)
    (.musyc_surface(theta + step, d1, d2) - .musyc_surface(theta - step, d1, d2)) / 2e-6
  }, numeric(length(d1)))
  expect_equal(unname(.musyc_surface(theta, d1, d2, jacobian = TRUE)), numeric_jacobian,
               tolerance = 1e-7)
})

test_that(".musyc_least_squares returns the best point it evaluated", {
  # a flat, noisy surface on which nlminb ends on singular convergence with a
  # par worse than the objective it reports
  conc <- c(0, 10^seq(-3, 1, length.out = 8))
  dt <- data.table::CJ(Concentration = conc, Concentration_2 = conc)[!(Concentration == 0 & Concentration_2 == 0)]
  set.seed(36)
  obs <- .musyc_grid(dt$Concentration, dt$Concentration_2, stats::runif(NROW(dt), 0.88, 1.07))
  limits <- .musyc_fit_limits(.musyc_bounds(-1, conc[-1], conc[-1]))
  set.seed(36003)
  start <- c(1, stats::runif(3, -1, 1.5), stats::rnorm(2, 1, 1), stats::rnorm(4, -1, 2))
  start <- pmin(pmax(start, limits$lower), limits$upper)
  residuals <- function(theta) .musyc_surface(theta, obs$d1, obs$d2) - obs$E
  fit <- .musyc_least_squares(
    residuals, function(theta) .musyc_surface(theta, obs$d1, obs$d2, jacobian = TRUE),
    start, limits$lower, limits$upper
  )
  expect_equal(fit$objective, sum(residuals(fit$par)^2))
  expect_lt(fit$objective, sum(residuals(start)^2))
})

test_that(".musyc_hill_fit recovers the parameters of a Hill curve", {
  d <- c(0, 10^seq(-3, 1, length.out = 9))
  E <- 1 + (0.2 - 1) * d^1.5 / (0.05^1.5 + d^1.5)
  fit <- .musyc_hill_fit(d, E)
  expect_equal(unlist(fit), c(E0 = 1, Emax = 0.2, h = 1.5, C = 0.05), tolerance = 1e-4)
})

test_that(".musyc_best_fit prefers a converged fit to a better one that ran out of iterations", {
  fit <- function(r2, converged) list(params = list(), r2 = r2, converged = converged)
  expect_equal(.musyc_best_fit(list(fit(0.9, TRUE), fit(0.95, FALSE), fit(0.8, TRUE)))$r2, 0.9)
  expect_equal(.musyc_best_fit(list(fit(0.9, FALSE), fit(0.95, FALSE)))$r2, 0.95)
  expect_equal(.musyc_best_fit(list(fit(NA_real_, FALSE), fit(0.7, TRUE)))$r2, 0.7)
  expect_true(is.na(.musyc_best_fit(list(fit(NA_real_, FALSE), fit(NA_real_, FALSE)))$r2))
})

test_that(".musyc_grid adds the untreated corner and drops missing responses", {
  obs <- .musyc_grid(c(0, 0.1, 1), c(0.1, 0, 1), c(0.5, NA, 0.2))
  expect_equal(NROW(obs), 3L)
  expect_equal(obs$E[obs$d1 == 0 & obs$d2 == 0], 1)
  expect_false(anyNA(obs$E))
})

test_that(".musyc_bounds lets E go negative for GR only and ties C to the doses", {
  doses <- c(0.01, 0.1, 1)
  expect_equal(.musyc_bounds(-1, doses, doses)$E3_bounds, c(-1, 1.5))
  expect_equal(.musyc_bounds(0, doses, doses)$E3_bounds, c(0, 1.5))
  expect_equal(.musyc_bounds(0, doses, doses)$E0_bounds, c(0.99, 1.01))
  expect_equal(.musyc_bounds(0, doses, c(1, 10))$C1_bounds, c(1e-5, 1e3))
  expect_equal(.musyc_bounds(0, doses, c(1, 10))$C2_bounds, c(1e-3, 1e4))
  expect_null(.musyc_bounds(0, doses, doses)$alpha12_bounds)
})

test_that(".musyc_fit_limits puts C on a log scale and leaves h and alpha free", {
  limits <- .musyc_fit_limits(.musyc_bounds(0, c(0.01, 1), c(0.1, 10)))
  expect_equal(unname(limits$lower[c("E0", "E3", "C1", "C2")]), c(0.99, 0, log(1e-5), log(1e-4)))
  expect_equal(unname(limits$upper[c("E1", "C1", "C2")]), c(1.5, log(1e3), log(1e4)))
  expect_true(all(is.infinite(limits$lower[c("h1", "h2", "alpha12", "alpha21")])))
  expect_true(all(is.infinite(limits$upper[c("h1", "h2", "alpha12", "alpha21")])))
})

test_that(".with_local_seed leaves the global random number generator untouched", {
  set.seed(7)
  before <- .Random.seed
  first <- .with_local_seed(1L, stats::runif(3))
  expect_identical(.Random.seed, before)
  expect_identical(.with_local_seed(1L, stats::runif(3)), first)
})

test_that("musyc_fit_fn stays close to the synergy fits of public Goetz pairs", {
  # averaged grids of Goetz et al., Cancers 2024, and the fits of synergy 1.0.0 on them
  grids <- data.table::fread(test_path("testdata", "musyc_goetz_grids.csv"))
  ref <- data.table::fread(test_path("testdata", "musyc_goetz_synergy.csv"))
  keys <- c("CellLineName", "DrugName", "DrugName_2", "normalization_type")
  fits <- grids[, suppressWarnings(musyc_fit_fn(.SD))[c("r2", "beta")], by = keys,
                .SDcols = c("Concentration", "Concentration_2", "x", "normalization_type")]
  fits <- merge(fits, ref, by = keys)
  expect_identical(NROW(fits), 4L)
  expect_true(all(fits$r2 > fits$r2_synergy - 0.02))
  # A2058 GR is flat: both fits reach a similar r2 with a different beta
  regular <- fits[fits$CellLineName != "A2058"]
  expect_true(all(abs(regular$beta - regular$beta_synergy) < 0.02))
})

test_that("musyc_fit_fn recovers the parameters of a known surface", {
  dt <- .musyc_truth_dt(.musyc_known_params)
  result <- musyc_fit_fn(dt)

  expect_true(result$converged)
  expect_gt(result$r2, 0.999)
  expect_equal(result$E3, .musyc_known_params$E3, tolerance = 0.05)
  expect_equal(log10(result$alpha12), log10(.musyc_known_params$alpha12), tolerance = 0.1)
  expect_identical(musyc_fit_fn(dt), result)
})

test_that("musyc_fit_fn recovers alpha for a drug more potent than any fixed bound", {
  params <- list(
    E0 = 1, E1 = 0.2, E2 = 0.3, E3 = 0, h1 = 1, h2 = 1,
    C1 = 5e-5, C2 = 5e-5, alpha12 = 3, alpha21 = 3
  )
  dt <- .musyc_truth_dt(params, "RV", conc = c(0, 10^seq(-7, -3, length.out = 9)))
  result <- musyc_fit_fn(dt)

  expect_gt(result$r2, 0.999)
  expect_equal(result$alpha12, 3, tolerance = 0.01)
  expect_equal(result$C1, 5e-5, tolerance = 0.01)
})

test_that("musyc_fit_fn fits through apply_fit", {
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

test_that("musyc_summary_fn clips extreme alpha values before averaging", {
  rows <- .fake_musyc_rows()
  rows[, alpha21 := c(2e-13, 1)]
  result <- musyc_summary_fn(rows)
  expect_equal(result$mean_log2_alpha, mean(c(1, 1, -10, 0)))
})
