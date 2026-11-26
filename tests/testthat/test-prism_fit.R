# The fit recipe of a fit source: gDRutils keys the recipes by fit_source, or by fit_id from 1.13.
.fit_recipe_of <- function(se, fit_source, assay = "Metrics") {
  fits <- gDRutils::get_SE_provenance(se)$fits[[assay]]
  if (!is.null(fits[[fit_source]])) {
    return(fits[[fit_source]])
  }
  Filter(function(f) identical(f$fit_source, fit_source), fits)[1][[1]]
}

# replicate points of one decreasing curve, upper limit 1, lower limit 0.2, inflection 0.1
prism_points_curve <- function(n_doses = 8L, inflection = 0.1, seed = 1L) {
  set.seed(seed)
  dose <- rep(10^seq(-3, 1, length.out = n_doses), each = 3)
  viability <- 0.2 + 0.8 / (1 + (dose / inflection)^1.5) + stats::rnorm(length(dose), sd = 0.03)
  list(dose = dose, viability = viability)
}

test_that(".prism_auc integrates the clamped curve over log2 dose", {
  fit <- list(upper_limit = 1, lower_limit = 0.2, inflection = 0.1, slope = -1.5)
  grid <- seq(log2(1e-3), log2(10), length.out = 1e5)
  expect_equal(.prism_auc(fit, 1e-3, 10), mean(.prism_curve(2^grid, fit)), tolerance = 1e-4)
  # the curve is clamped to [0, 1]
  expect_equal(.prism_auc(list(upper_limit = 1.3, lower_limit = 1.2, inflection = 1, slope = -1), 0.01, 1), 1)
  expect_equal(.prism_auc(list(upper_limit = -0.1, lower_limit = -0.2, inflection = 1, slope = -1), 0.01, 1), 0)
})

test_that(".prism_log2_ic50 gives the crossing of 0.5 within the tested doses only", {
  fit <- list(upper_limit = 1, lower_limit = 0, inflection = 0.1, slope = -1.5)
  expect_equal(.prism_log2_ic50(fit, 1e-3, 10), log2(0.1), tolerance = 1e-4)
  # outside the tested doses
  expect_true(is.na(.prism_log2_ic50(fit, 1, 10)))
  # a lower limit above 0.5 or an upper limit below it
  expect_true(is.na(.prism_log2_ic50(modifyList(fit, list(lower_limit = 0.6)), 1e-3, 10)))
  expect_true(is.na(.prism_log2_ic50(modifyList(fit, list(upper_limit = 0.4)), 1e-3, 10)))
  expect_true(is.na(.prism_log2_ic50(modifyList(fit, list(lower_limit = NA)), 1e-3, 10)))
})

test_that("fit_prism_curve recovers a curve and keeps the fit that explains the most variance", {
  skip_if_not(prism_fit_available(), "dr4pl is not installed")
  p <- prism_points_curve()
  fit <- fit_prism_curve(p$dose, p$viability)

  expect_s3_class(fit, "data.table")
  expect_identical(names(fit), c("fit_name", "lower_limit", "upper_limit", "slope", "inflection", "mse",
                                 "mad", "frac_var_explained", "successful_fit", "auc_riemann",
                                 "minimum_dose", "maximum_dose", "auc", "log2_auc", "log2_ic50", "ic50"))
  expect_true(fit$successful_fit)
  expect_equal(fit$lower_limit, 0.2, tolerance = 0.1)
  expect_equal(fit$inflection, 0.1, tolerance = 0.3)
  expect_lt(fit$slope, 0)
  expect_gt(fit$frac_var_explained, 0.9)
  expect_equal(fit$auc, .prism_auc(as.list(fit), min(p$dose), max(p$dose)))
  expect_equal(fit$auc_riemann, mean(pmin(p$viability, 1)))
  expect_equal(fit$ic50, 2^fit$log2_ic50)

  # the kept fit is the accepted candidate with the highest explained variance
  candidates <- .prism_candidates(p$dose, pmin(p$viability, 1.5), c(0.8, 1.01))
  explained <- vapply(candidates, function(f) {
    1 - mean((pmin(p$viability, 1.5) - .prism_curve(p$dose, f))^2) / stats::var(pmin(p$viability, 1.5))
  }, numeric(1))
  expect_identical(fit$fit_name, names(candidates)[which.max(explained)])
  expect_equal(fit$frac_var_explained, max(explained))
})

test_that("fit_prism_curve caps viability and needs enough doses", {
  skip_if_not(prism_fit_available(), "dr4pl is not installed")
  p <- prism_points_curve()
  high <- p$viability
  high[1:3] <- 3
  expect_equal(fit_prism_curve(p$dose, high)$auc_riemann, mean(pmin(high, 1)))

  # four doses: no fit, the data values only
  p4 <- prism_points_curve(n_doses = 4L)
  fit4 <- fit_prism_curve(p4$dose, p4$viability)
  expect_false(fit4$successful_fit)
  expect_true(is.na(fit4$fit_name) && is.na(fit4$auc) && is.na(fit4$ic50))
  expect_equal(fit4$minimum_dose, min(p4$dose))

  # five doses: a fit, but no IC50 below six doses
  p5 <- prism_points_curve(n_doses = 5L)
  fit5 <- fit_prism_curve(p5$dose, p5$viability)
  expect_true(fit5$successful_fit)
  expect_true(is.na(fit5$log2_ic50))
  expect_false(is.na(fit_prism_curve(p5$dose, p5$viability, ic50_min_doses = 5L)$log2_ic50))

  # doses of 0 and missing values are left out
  fit0 <- fit_prism_curve(c(0, NA, p$dose), c(1, 1, p$viability))
  expect_equal(fit0$auc, fit_prism_curve(p$dose, p$viability)$auc)
})

test_that("refit_PRISM_curves fits one curve per barcode and leaves out curves with too few doses", {
  skip_if_not(prism_fit_available(), "dr4pl is not installed")
  curve <- function(lua, n_doses, inflection) {
    p <- prism_points_curve(n_doses, inflection)
    data.table::data.table(depmap_id = "ACH-000004", lua = lua, pert_name = "G00106", day = 5,
                           pert_dose = p$dose, l2fc = log2(p$viability), clid = "CL1")
  }
  points <- rbind(curve("LUA-1", 8L, 0.1), curve("LUA-2", 8L, 1), curve("LUA-3", 4L, 0.1))
  data.table::setattr(points, "gDR_import", list(files = data.frame(role = "prism_l2fc")))

  curves <- refit_PRISM_curves(points, keep_cols = "clid")
  expect_identical(curves$lua, c("LUA-1", "LUA-2"))
  expect_identical(curves$clid, c("CL1", "CL1"))
  expect_true(all(curves$successful_fit))
  expect_lt(curves$inflection[1], curves$inflection[2])
  expect_identical(attr(curves, "gDR_import")$files$role, "prism_l2fc")
  expect_equal(curves$auc[1], fit_prism_curve(points$pert_dose[1:24], 2^points$l2fc[1:24])$auc)
  # without further columns
  plain <- refit_PRISM_curves(points)
  expect_identical(names(plain), c("depmap_id", "lua", "pert_name", "day", names(fit_prism_curve(1, 1))))
  expect_equal(plain$auc, curves$auc)
})

# replicate points for the three cell lines of finalMAE_prism used in test-prism_metrics.R
prism_refit_points <- function() {
  clids <- c("HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
             "MONOMAC6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "NOT_IN_THE_DATA")
  data.table::rbindlist(lapply(seq_along(clids), function(i) {
    p <- prism_points_curve(inflection = 0.05 * i, seed = i)
    data.table::data.table(pool_id = "P1", depmap_id = sprintf("ACH-00000%d", i + 3), lua = sprintf("LUA-%d", i),
                           pert_name = "G00106", day = 5, pert_dose = p$dose, l2fc = log2(p$viability),
                           screen = "TEST_SCREEN", clid = clids[i], Duration = 120)
  }))
}

test_that("add_PRISM_refit_metrics adds the refitted curves next to the gDR and PRISM fits", {
  skip_if_not(prism_fit_available(), "dr4pl is not installed")
  se <- MultiAssayExperiment::experiments(gDRutils::get_synthetic_data("finalMAE_prism"))[["single-agent"]]
  gdr_before <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  points <- prism_refit_points()

  out <- expect_warning(add_PRISM_refit_metrics(se, points), "1 of 4 PRISM curves have no matching cell")
  metrics <- gDRutils::convert_se_assay_to_dt(out, "Metrics")
  refit <- metrics[metrics$fit_source == "PRISM_refit", ]
  expect_identical(NROW(refit), 3L)
  expect_identical(unique(refit$normalization_type), "RV")
  expect_true(all(refit$h > 0))
  expect_true(all(!is.na(refit$prism_auc)))
  expect_identical(NROW(metrics[metrics$fit_source == "gDR", ]), NROW(gdr_before))

  steps <- gDRutils::get_SE_provenance(out)$steps
  expect_identical(steps$prism_refit$settings$n_added, 3L)
  expect_identical(steps$prism_refit$settings$n_left_out, 1L)
  expect_identical(steps$prism_refit$settings$cap, 1.5)
  expect_true(all(c("drc", "dr4pl") %in% names(steps$prism_refit$settings$packages)))
  expect_false(is.null(.fit_recipe_of(out, "PRISM_refit")))

  # running it again replaces the refitted rows
  again <- suppressWarnings(add_PRISM_refit_metrics(out, points))
  again_metrics <- gDRutils::convert_se_assay_to_dt(again, "Metrics")
  expect_identical(NROW(again_metrics), NROW(metrics))
})

test_that("prism_fit_fn fits gDR-processed data through apply_fit", {
  skip_if_not(prism_fit_available(), "dr4pl is not installed")
  se <- MultiAssayExperiment::experiments(gDRutils::get_synthetic_data("finalMAE_prism"))[["single-agent"]]
  out <- gDRutils::apply_fit(se, prism_fit_fn, "single-agent", input_assay = "Normalized",
                             output_assay = "Metrics", slicing_values = "RV", fit_source = "PRISM_refit")
  metrics <- gDRutils::convert_se_assay_to_dt(out, "Metrics")
  refit <- metrics[metrics$fit_source == "PRISM_refit", ]
  gdr <- metrics[metrics$fit_source == "gDR" & metrics$normalization_type == "RV", ]

  expect_identical(NROW(refit), NROW(gdr))
  expect_identical(unique(refit$normalization_type), "RV")
  expect_true(all(refit$N_conc > 0))
  expect_true(any(!is.na(refit$prism_auc)))
  recipe <- .fit_recipe_of(out, "PRISM_refit")
  expect_false(is.null(recipe))
  expect_true(any(grepl("dr4pl", c(names(unlist(recipe)), unlist(recipe)))))
})
