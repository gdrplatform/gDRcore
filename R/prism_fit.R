.ref_prism_fit <- paste(
  "PRISM dose-response fitting as in the PRISM sequencing pipeline (theprismlab/sushi,",
  "get_best_fit) and the DepMap PRISM Oncology Reference; this is an independent",
  "implementation of the published algorithm."
)

# The columns that identify a curve in PRISM's replicate-level data: the cell line barcode
# (pool, LUA, cell set) and the compound plate, as PRISM fits them.
.prism_curve_cols <- c("pool_id", "depmap_id", "lua", "cell_set", "growth_condition", "pert_type",
                       "pert_name", "pert_id", "day", "x_project_id", "pert_plate", "pert_vehicle")


#' Fit a dose-response curve with PRISM's algorithm
#'
#' Fits the 4-parameter log-logistic curve
#' \code{upper_limit + (lower_limit - upper_limit) / (1 + (dose / inflection)^slope)} the way
#' PRISM fits its screens: every replicate point is used, viability is capped at \code{cap}, and
#' five fits are tried in this order, keeping the accepted one that explains the most variance
#' (the earlier one on a tie):
#' \enumerate{
#'   \item \code{drc_drm_constrained}: \code{drc::drm()} with \code{LL.4()}, the upper limit in
#'     \code{ul_bounds}, the lower limit in [0, 1.01] and a decreasing curve;
#'   \item \code{drc_drm_unconstrained}: the same with the slope unbounded, kept only when the
#'     curve decreases;
#'   \item \code{dr4pl_initMan_constrained_optNM}: \code{dr4pl::dr4pl()} from
#'     (1, 8 x the lowest dose, -3, 0.01) with Nelder-Mead and the bounds of fit 1;
#'   \item \code{dr4pl_initL_unconstrained}: \code{dr4pl} with the logistic start, upper limit
#'     1 (in [0.99, 1.01]) and lower limit 0.3, kept when the curve decreases;
#'   \item \code{dr4pl_initL_constrained}: \code{dr4pl} with the logistic start, upper limit 1
#'     and lower limit 0.005, and the bounds of fit 1.
#' }
#' A \code{drc} fit is accepted when it converged with finite parameters; a \code{dr4pl} fit that
#' did not converge is replaced by its robust refit, and fits 3 and 5 are then accepted when they
#' converged with finite parameters.
#'
#' \code{auc} is the integral over log2 dose of the curve clamped to [0, 1], divided by the
#' log2 dose range; \code{auc_riemann} is the mean of the points capped at 1. \code{log2_ic50} is
#' where the curve crosses 0.5 within the tested doses, and \code{NA} when it does not, when the
#' lower limit is at least 0.5 or the upper limit at most 0.5, or with fewer than
#' \code{ic50_min_doses} doses.
#'
#' Needs the \code{dr4pl} package (see \code{\link{prism_fit_available}}).
#'
#' @param dose numeric vector of doses (above 0), one per point.
#' @param viability numeric vector of viabilities (e.g. \code{2^l2fc}), one per point.
#' @param cap number, the viabilities above it are set to it before the fit.
#' @param ul_bounds numeric vector of two, the bounds of the upper limit.
#' @param min_doses integer, the smallest number of distinct doses fitted.
#' @param ic50_min_doses integer, the smallest number of distinct doses that gives an IC50.
#'
#' @return One-row \code{data.table} with the columns of PRISM's DRC table: \code{fit_name},
#'   \code{lower_limit}, \code{upper_limit}, \code{slope}, \code{inflection}, \code{mse},
#'   \code{mad}, \code{frac_var_explained}, \code{successful_fit}, \code{auc_riemann},
#'   \code{minimum_dose}, \code{maximum_dose}, \code{auc}, \code{log2_auc}, \code{log2_ic50} and
#'   \code{ic50}. Without an accepted fit, \code{successful_fit} is \code{FALSE} and the curve
#'   values are \code{NA}.
#'
#' @examples
#' if (prism_fit_available()) {
#'   dose <- rep(10^seq(-3, 1, length.out = 8), each = 3)
#'   viability <- 0.2 + 0.8 / (1 + (dose / 0.1)^1.5) + stats::rnorm(24, sd = 0.03)
#'   fit_prism_curve(dose, viability)
#' }
#'
#' @keywords prism_fit
#' @export
fit_prism_curve <- function(dose,
                            viability,
                            cap = 1.5,
                            ul_bounds = c(0.8, 1.01),
                            min_doses = 5L,
                            ic50_min_doses = 6L) {
  checkmate::assert_numeric(dose)
  checkmate::assert_numeric(viability, len = length(dose))
  checkmate::assert_number(cap)
  checkmate::assert_numeric(ul_bounds, len = 2L, any.missing = FALSE, sorted = TRUE)
  checkmate::assert_count(min_doses, positive = TRUE)
  checkmate::assert_count(ic50_min_doses, positive = TRUE)
  .assert_prism_fit_available()

  keep <- is.finite(dose) & dose > 0 & is.finite(viability)
  dose <- dose[keep]
  y <- pmin(viability[keep], cap)
  out <- data.table::data.table(
    fit_name = NA_character_, lower_limit = NA_real_, upper_limit = NA_real_, slope = NA_real_,
    inflection = NA_real_, mse = NA_real_, mad = NA_real_, frac_var_explained = NA_real_,
    successful_fit = FALSE,
    auc_riemann = if (length(y)) mean(pmin(y, 1)) else NA_real_,
    minimum_dose = if (length(dose)) min(dose) else NA_real_,
    maximum_dose = if (length(dose)) max(dose) else NA_real_,
    auc = NA_real_, log2_auc = NA_real_, log2_ic50 = NA_real_, ic50 = NA_real_
  )
  n_doses <- data.table::uniqueN(dose)
  if (n_doses < min_doses) {
    return(out)
  }
  fits <- .prism_candidates(dose, y, ul_bounds)
  if (!length(fits)) {
    return(out)
  }

  var_y <- stats::var(y)
  fit_stats <- lapply(fits, function(fit) {
    res <- y - .prism_curve(dose, fit)
    mse <- mean(res^2)
    c(mse = mse, mad = stats::median(abs(res)), frac_var_explained = 1 - mse / var_y)
  })
  explained <- vapply(fit_stats, function(s) s[["frac_var_explained"]], numeric(1))
  best <- if (all(is.na(explained))) 1L else which.max(explained)
  fit <- fits[[best]]

  out$fit_name <- names(fits)[best]
  out$lower_limit <- fit$lower_limit
  out$upper_limit <- fit$upper_limit
  out$slope <- fit$slope
  out$inflection <- fit$inflection
  out$mse <- fit_stats[[best]][["mse"]]
  out$mad <- fit_stats[[best]][["mad"]]
  out$frac_var_explained <- fit_stats[[best]][["frac_var_explained"]]
  out$successful_fit <- TRUE
  out$auc <- .prism_auc(fit, out$minimum_dose, out$maximum_dose)
  out$log2_auc <- log2(out$auc)
  if (n_doses >= ic50_min_doses) {
    out$log2_ic50 <- .prism_log2_ic50(fit, out$minimum_dose, out$maximum_dose)
    out$ic50 <- 2^out$log2_ic50
  }
  out
}


#' Whether the PRISM fit can run
#'
#' \code{\link{fit_prism_curve}} and the functions built on it fit with the \code{drc} and
#' \code{dr4pl} packages; \code{dr4pl} is only suggested by gDRcore.
#'
#' @return Logical scalar.
#'
#' @examples
#' prism_fit_available()
#'
#' @keywords prism_fit
#' @export
prism_fit_available <- function() {
  requireNamespace("dr4pl", quietly = TRUE) && requireNamespace("drc", quietly = TRUE)
}


# Stops with an install hint when dr4pl is missing.
#' @keywords internal
.assert_prism_fit_available <- function() {
  if (!prism_fit_available()) {
    stop("the PRISM fit needs the 'dr4pl' package: install.packages(\"dr4pl\")", call. = FALSE)
  }
}


# The package versions of the PRISM fit, for the provenance.
#' @keywords internal
.prism_fit_runtime <- function() {
  vapply(c("drc", "dr4pl"), function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))
}


# PRISM's curve at the given doses.
#' @keywords internal
.prism_curve <- function(dose, fit) {
  fit$upper_limit + (fit$lower_limit - fit$upper_limit) / (1 + (dose / fit$inflection)^fit$slope)
}


# The accepted candidate fits, named by fit_name, in PRISM's order.
#' @keywords internal
.prism_candidates <- function(dose, y, ul_bounds) {
  # the slope of a decreasing curve is negative; the bound leaves room for flat fits
  slope_bound <- 1e-5
  constrained_lower <- c(ul_bounds[1], -Inf, -Inf, 0)
  constrained_upper <- c(ul_bounds[2], Inf, slope_bound, 1.01)
  fits <- list(
    drc_drm_constrained = .prism_drc_fit(dose, y, ul_bounds, min_b = -slope_bound),
    drc_drm_unconstrained = .prism_drc_fit(dose, y, ul_bounds, min_b = -Inf, decreasing = TRUE),
    dr4pl_initMan_constrained_optNM = .prism_dr4pl_fit(
      dose, y, start = c(1, 8 * min(dose), -3, 0.01), init_method = "Mead",
      lower = constrained_lower, upper = constrained_upper, accept = "converged"
    ),
    dr4pl_initL_unconstrained = .prism_dr4pl_fit(
      dose, y, start = c(1, NA, NA, 0.3), init_method = "logistic",
      lower = c(0.99, -Inf, -Inf, 0), upper = c(1.01, Inf, Inf, 1.01),
      accept = "decreasing", slope_bound = slope_bound
    ),
    dr4pl_initL_constrained = .prism_dr4pl_fit(
      dose, y, start = c(1, NA, NA, 0.005), init_method = "logistic",
      lower = constrained_lower, upper = constrained_upper, accept = "converged"
    )
  )
  Filter(Negate(is.null), fits)
}


# A drc LL.4 fit. drc's curve is c + (d - c) / (1 + (dose / e)^b): b is minus PRISM's slope,
# c the lower and d the upper limit. NULL when the fit fails or is not accepted.
#' @keywords internal
.prism_drc_fit <- function(dose, y, ul_bounds, min_b, decreasing = FALSE) {
  fit <- .quietly(drc::drm(
    y ~ dose, data = data.frame(dose = dose, y = y), fct = drc::LL.4(),
    lowerl = c(min_b, 0, ul_bounds[1], -Inf), upperl = c(Inf, 1.01, ul_bounds[2], Inf)
  ))
  if (is.null(fit)) {
    return(NULL)
  }
  p <- unname(stats::coef(fit))
  if (!isTRUE(fit$fit$convergence) || !all(is.finite(p)) || (decreasing && !(p[1] > 0))) {
    return(NULL)
  }
  list(lower_limit = p[2], upper_limit = p[3], slope = -p[1], inflection = p[4])
}


# A dr4pl fit; dr4pl's parameters are the upper limit, IC50, slope and lower limit of PRISM's
# curve. A fit that did not converge is replaced by its robust refit. NULL when the fit fails or
# is not accepted.
#' @keywords internal
.prism_dr4pl_fit <- function(dose, y, start, init_method, lower, upper,
                             accept = c("converged", "decreasing"), slope_bound = 1e-5) {
  accept <- match.arg(accept)
  fit <- .quietly(dr4pl::dr4pl(
    dose, y, init.parm = dr4pl::dr4pl_theta(start[1], start[2], start[3], start[4]),
    method.init = init_method, method.optim = "Nelder-Mead", lowerl = lower, upperl = upper
  ))
  if (is.null(fit)) {
    return(NULL)
  }
  if (!isTRUE(fit$convergence) && !is.null(fit$dr4pl.robust)) {
    fit <- fit$dr4pl.robust
  }
  p <- unname(as.numeric(fit$parameters))
  ok <- switch(accept,
    converged = isTRUE(fit$convergence) && length(p) == 4L && all(is.finite(p)),
    decreasing = length(p) == 4L && isTRUE(p[3] < slope_bound)
  )
  if (!ok) {
    return(NULL)
  }
  list(lower_limit = p[4], upper_limit = p[1], slope = p[3], inflection = p[2])
}


# Evaluates a fit call without its output: NULL on an error; the errors drc and dr4pl print
# from their own try() calls are dropped.
#' @keywords internal
.quietly <- function(expr) {
  out <- NULL
  utils::capture.output(
    out <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) NULL),
    type = "message"
  )
  out
}


# PRISM's AUC: the curve clamped to [0, 1], integrated over log2 dose and divided by the range.
#' @keywords internal
.prism_auc <- function(fit, min_dose, max_dose) {
  span <- log2(max_dose) - log2(min_dose)
  clamped <- function(x) pmax(pmin(.prism_curve(2^x, fit), 1), 0)
  tryCatch(
    stats::integrate(clamped, log2(min_dose), log2(max_dose))$value / span,
    error = function(e) NA_real_
  )
}


# log2 of the dose where the curve crosses 0.5 within the tested doses, or NA.
#' @keywords internal
.prism_log2_ic50 <- function(fit, min_dose, max_dose) {
  if (!isTRUE(fit$lower_limit < 0.5) || !isTRUE(fit$upper_limit > 0.5)) {
    return(NA_real_)
  }
  tryCatch(
    stats::uniroot(function(x) .prism_curve(2^x, fit) - 0.5, c(log2(min_dose), log2(max_dose)))$root,
    error = function(e) NA_real_
  )
}


#' Refit PRISM's replicate-level data with PRISM's algorithm
#'
#' Fits every curve of PRISM's replicate-level log2 fold changes (read with
#' \code{gDRimport::read_PRISM_l2fc()}) with \code{\link{fit_prism_curve}}, one curve per cell
#' line barcode, compound and plate as PRISM fits them. Curves with fewer than \code{min_doses}
#' doses are left out, as in PRISM's pipeline. The result has the shape of PRISM's DRC table, so
#' it can be compared with a delivered table or added to the \code{Metrics} assay with
#' \code{\link{add_PRISM_refit_metrics}}.
#'
#' @param prism_points \code{data.table} with one row per replicate point: the curve columns,
#'   \code{dose_col} and \code{l2fc_col}.
#' @param curve_cols character vector, the columns that identify a curve (those present are used);
#'   by default PRISM's: \code{pool_id}, \code{depmap_id}, \code{lua}, \code{cell_set},
#'   \code{growth_condition}, \code{pert_type}, \code{pert_name}, \code{pert_id}, \code{day},
#'   \code{x_project_id}, \code{pert_plate} and \code{pert_vehicle}.
#' @param dose_col,l2fc_col strings, the columns with the dose and the log2 fold change.
#' @param keep_cols character vector, further columns kept with each curve (one value per curve),
#'   e.g. the cell line identifier and Duration added by \code{gDRimport::read_PRISM_l2fc()}.
#' @inheritParams fit_prism_curve
#'
#' @return \code{data.table} with one row per curve: \code{curve_cols}, \code{keep_cols} and the
#'   columns of \code{\link{fit_prism_curve}}. The \code{gDR_import} attribute of
#'   \code{prism_points} (the input files) is kept.
#'
#' @examples
#' if (prism_fit_available()) {
#'   points <- data.table::data.table(
#'     depmap_id = "ACH-000001", pert_name = "drug A", day = 5,
#'     pert_dose = rep(10^seq(-3, 1, length.out = 8), each = 3)
#'   )
#'   points$l2fc <- log2(0.2 + 0.8 / (1 + (points$pert_dose / 0.1)^1.5))
#'   refit_PRISM_curves(points)
#' }
#'
#' @keywords prism_fit
#' @export
refit_PRISM_curves <- function(prism_points,
                               curve_cols = .prism_curve_cols,
                               dose_col = "pert_dose",
                               l2fc_col = "l2fc",
                               keep_cols = NULL,
                               cap = 1.5,
                               ul_bounds = c(0.8, 1.01),
                               min_doses = 5L,
                               ic50_min_doses = 6L) {
  checkmate::assert_multi_class(prism_points, c("data.table", "data.frame"))
  checkmate::assert_character(curve_cols, min.len = 1, any.missing = FALSE)
  checkmate::assert_string(dose_col)
  checkmate::assert_string(l2fc_col)
  checkmate::assert_character(keep_cols, any.missing = FALSE, null.ok = TRUE)
  .assert_prism_fit_available()

  points <- data.table::as.data.table(prism_points)
  checkmate::assert_names(names(points), must.include = c(dose_col, l2fc_col))
  keys <- intersect(curve_cols, names(points))
  if (!length(keys)) {
    stop("none of 'curve_cols' is a column of 'prism_points'")
  }
  extra <- setdiff(intersect(keep_cols, names(points)), keys)
  keep <- is.finite(points[[dose_col]]) & is.finite(points[[l2fc_col]])
  points <- points[keep]
  points$.dose <- points[[dose_col]]
  points$.viability <- 2^points[[l2fc_col]]
  n_doses <- points[, list(n = data.table::uniqueN(.dose)), by = keys]
  fitted <- n_doses[n_doses$n >= min_doses, keys, with = FALSE]
  points <- points[fitted, on = keys]

  fit_one <- function(dose, viability) {
    fit_prism_curve(dose, viability, cap = cap, ul_bounds = ul_bounds, min_doses = min_doses,
                    ic50_min_doses = ic50_min_doses)
  }
  curves <- if (length(extra)) {
    points[, c(lapply(.SD, `[`, 1L), fit_one(.dose, .viability)), by = keys, .SDcols = extra]
  } else {
    points[, fit_one(.dose, .viability), by = keys]
  }
  data.table::setattr(curves, "gDR_import", attr(prism_points, "gDR_import", exact = TRUE))
  curves
}


#' Add PRISM-algorithm fits of PRISM's replicate-level data to the Metrics assay
#'
#' Refits PRISM's replicate-level data with \code{\link{refit_PRISM_curves}} and adds the curves
#' to the \code{Metrics} assay of processed PRISM data as \code{fit_source}, next to the gDR fits
#' and the curves PRISM delivered (\code{\link{add_PRISM_metrics}}, \code{fit_source = "PRISM"}).
#' The curves are matched to the cells as in \code{\link{add_PRISM_metrics}}, and the columns of a
#' curve are mapped onto the gDR model the same way.
#'
#' Run it after \code{\link{fit_SE}}, which replaces the whole \code{Metrics} assay.
#'
#' @param se \code{SummarizedExperiment} of single-agent PRISM data, processed by gDR.
#' @param prism_points \code{data.table} of replicate points, as returned by
#'   \code{gDRimport::read_PRISM_l2fc()}: PRISM's curve columns, \code{pert_dose}, \code{l2fc},
#'   the cell line identifier and Duration.
#' @param fit_source string, the fit source of the refitted curves.
#' @inheritParams add_PRISM_metrics
#' @inheritParams fit_prism_curve
#'
#' @return \code{se} with the refitted curves in \code{metrics_assay}. The fit is recorded in the
#'   provenance under \code{fit_source}, and the step \code{"prism_refit"} holds the input files,
#'   the settings, the \code{drc} and \code{dr4pl} versions and the number of curves fitted,
#'   added and left out.
#'
#' @examples
#' if (prism_fit_available() && requireNamespace("gDRimport", quietly = TRUE) &&
#'     exists("read_PRISM_l2fc", asNamespace("gDRimport"))) {
#'   mae <- gDRutils::get_synthetic_data("finalMAE_prism")
#'   se <- MultiAssayExperiment::experiments(mae)[["single-agent"]]
#'   points <- gDRimport::read_PRISM_l2fc(
#'     l2fc_path = system.file("testdata/prism_l2fc.csv", package = "gDRimport"),
#'     meta_data_path = system.file("testdata/prism_drc_model.csv", package = "gDRimport")
#'   )
#'   se <- add_PRISM_refit_metrics(se, points)
#'   metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#'   metrics[metrics$fit_source == "PRISM_refit", ]
#' }
#'
#' @keywords prism_fit
#' @export
add_PRISM_refit_metrics <- function(se,
                                    prism_points,
                                    drug_map = NULL,
                                    cap = 1.5,
                                    normalization_type = "RV",
                                    fit_source = "PRISM_refit",
                                    metrics_assay = "Metrics",
                                    averaged_assay = "Averaged") {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_multi_class(prism_points, c("data.table", "data.frame"))
  checkmate::assert_character(drug_map, names = "unique", any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_number(cap)
  checkmate::assert_string(normalization_type)
  checkmate::assert_string(fit_source)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_string(averaged_assay)

  idfs <- gDRutils::get_SE_identifiers(se, simplify = TRUE)
  curves <- refit_PRISM_curves(prism_points, keep_cols = c(idfs[["cellline"]], idfs[["duration"]], "screen"),
                               cap = cap)
  added <- .add_prism_curves(se, curves, drug_map, normalization_type, metrics_assay,
                             averaged_assay, fit_source = fit_source)
  se <- added$se

  se <- gDRutils::record_fit_provenance(
    se, assay = metrics_assay, fit_source = fit_source, fn = add_PRISM_refit_metrics,
    fit_args = list(drug_map = drug_map, cap = cap, normalization_type = normalization_type,
                    fit_source = fit_source, metrics_assay = metrics_assay,
                    averaged_assay = averaged_assay)
  )
  files <- attr(prism_points, "gDR_import", exact = TRUE)$files
  gDRutils::record_processing_step(
    se, "prism_refit",
    settings = list(
      source = "PRISM replicate-level log2 fold changes",
      files = if (!is.null(files)) as.list(files[, intersect(c("role", "name", "md5"), names(files))]),
      screen = if ("screen" %in% names(curves)) unique(curves$screen),
      cap = cap,
      packages = as.list(.prism_fit_runtime()),
      n_curves = added$n_curves,
      n_fitted = sum(curves$successful_fit),
      n_added = sum(added$matched),
      n_left_out = sum(!added$matched)
    )
  )
}


#' Fit function: PRISM's algorithm for gDR-processed data
#'
#' A fit function for \code{gDRutils::apply_fit()} that fits one cell of gDR-processed
#' single-agent data with \code{\link{fit_prism_curve}}, using every replicate point. Use it on
#' the \code{Normalized} assay, the relative viability only:
#'
#' \preformatted{gDRutils::apply_fit(se, prism_fit_fn, "single-agent", input_assay = "Normalized",
#'                     output_assay = "Metrics", slicing_values = "RV",
#'                     fit_source = "PRISM_refit")}
#'
#' This applies PRISM's method to the data as gDR processed them (after masking, with gDR's curve
#' identity). It does not reproduce PRISM's numbers: gDR's PRISM import caps viability at 1.03
#' (PRISM fits up to 1.5), and gDR merges the curves of a cell line's barcodes and compound plates
#' into one cell. For PRISM's own curves, refit PRISM's file with
#' \code{\link{add_PRISM_refit_metrics}}.
#'
#' @param dt \code{data.table} with the points of one cell: the concentration column and \code{x}
#'   (and \code{normalization_type}).
#' @param ... further arguments of \code{\link{fit_prism_curve}} (\code{cap}, \code{ul_bounds},
#'   \code{min_doses}, \code{ic50_min_doses}).
#'
#' @return One-row \code{data.table} with the \code{Metrics} columns of the curve (\code{x_0},
#'   \code{x_inf}, \code{ec50}, \code{h}, \code{r2}, \code{fit_type}, \code{N_conc}, ...) and
#'   PRISM's values (\code{prism_auc}, \code{prism_ic50}, ...), mapped as in
#'   \code{\link{add_PRISM_metrics}}.
#'
#' @examples
#' if (prism_fit_available()) {
#'   mae <- gDRutils::get_synthetic_data("finalMAE_prism")
#'   se <- MultiAssayExperiment::experiments(mae)[["single-agent"]]
#'   se <- gDRutils::apply_fit(se, prism_fit_fn, "single-agent", input_assay = "Normalized",
#'                             output_assay = "Metrics", slicing_values = "RV",
#'                             fit_source = "PRISM_refit")
#' }
#'
#' @keywords prism_fit
#' @export
prism_fit_fn <- function(dt, ...) {
  conc <- gDRutils::get_env_identifiers("concentration")
  dt <- data.table::as.data.table(dt)
  checkmate::assert_names(names(dt), must.include = c(conc, "x"))
  norm_type <- if ("normalization_type" %in% names(dt)) dt$normalization_type[1] else "RV"
  points <- dt[dt[[conc]] > 0 & is.finite(dt$x)]
  curve <- fit_prism_curve(points[[conc]], points$x, ...)
  out <- .prism_curves_to_metrics(curve, normalization_type = norm_type)
  out$fit_source <- NULL
  out$N_conc <- data.table::uniqueN(points[[conc]])
  out
}

prism_fit_fn <- gDRutils::set_fit_fn_info(
  prism_fit_fn,
  model = "PRISM 4-parameter log-logistic (best of five drc and dr4pl fits)",
  description = paste(
    "Fits every replicate point with viability capped at 1.5, tries two drc and three dr4pl",
    "fits of the 4-parameter log-logistic curve, keeps the one that explains the most variance,",
    "and gives PRISM's AUC (curve clamped to [0, 1] over log2 dose) and IC50."
  ),
  reference = .ref_prism_fit,
  env = function() as.list(.prism_fit_runtime())
)
