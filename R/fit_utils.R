####
# Public API
####

#' apply_fit
#'
#' Generic layer for applying a user-supplied fit function across all
#' (row x column x slicing_col) triplets in a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} and persisting
#' the results into a named output assay.
#'
#' @details
#' The experiment type is declared via \code{data_type}, which resolves a
#' built-in slicing profile (default \code{slicing_cols}, \code{slicing_values},
#' and \code{input_assay}).  Any profile field can be overridden with the
#' corresponding argument.
#'
#' An optional \code{summary_fn} may be provided.  It is called once per
#' (row x column) cell and receives all rows that \code{fit_fn} produced for
#' that cell, returning a single summary row written to \code{summary_assay}.
#' This is the right place for aggregated metrics (e.g. mean synergy across
#' normalization types) that span multiple slices.
#'
#' Use the pipe to apply several custom fits to the same SE:
#' \preformatted{
#'   se |>
#'     apply_fit(bliss_fn, "combination",
#'                      output_assay = "custom_bliss", ...) |>
#'     apply_fit(musyc_fn,  "combination",
#'                      output_assay = "musyc_params",
#'                      summary_fn   = musyc_summary_fn,
#'                      summary_assay = "musyc_summary", ...)
#' }
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @param fit_fn function(\code{data.table}) returning a named list or
#'   single-row data.frame.  Receives data already filtered to one
#'   \code{slicing_cols} combination.
#' @param data_type one of \code{"single-agent"}, \code{"combination"},
#'   \code{"time-course"}.  Resolves profile defaults.
#' @param slicing_cols character vector of column name(s) used to split each
#'   BumpyMatrix cell into sub-experiments.  \code{NULL} uses the profile
#'   default.
#' @param slicing_values character vector of values to iterate over in
#'   \code{slicing_cols}.  \code{NULL} iterates all unique values found in
#'   each cell.
#' @param input_assay name of the source BumpyMatrix assay.  \code{NULL} uses
#'   the profile default.
#' @param output_assay name of the assay to write fit results into.  Required.
#'   Any assay name is accepted, including the native \code{"Metrics"} assay.
#'   With \code{merge = "merge"} (default), existing rows keyed by
#'   \code{fit_source} are replaced while all other rows are preserved —
#'   safe for co-existing alongside native gDR metrics.
#'   With \code{merge = "replace"}, the entire assay is overwritten.
#' @param summary_fn optional function(\code{data.table}) → named list called
#'   once per (row x column) cell on all rows returned by \code{fit_fn} for
#'   that cell.  Requires \code{summary_assay}.
#' @param summary_assay name of the assay to write summary results into.
#'   Required when \code{summary_fn} is provided.
#' @param merge \code{"merge"} (idempotent upsert keyed by \code{fit_source}
#'   + \code{slicing_cols}) or \code{"replace"} (overwrite the whole assay).
#' @param on_error \code{"warn"} (skip failed cells with a warning) or
#'   \code{"stop"} (propagate the error).
#' @param fit_source character string recorded in the \code{fit_source} column
#'   of every output row.  Forms part of the upsert key.
#'
#' @return updated \code{SummarizedExperiment} with \code{output_assay} (and
#'   optionally \code{summary_assay}) added or updated.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' mean_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
#' se_out <- apply_fit(
#'   se, mean_fn, "single-agent",
#'   output_assay = "custom_mean", fit_source = "demo"
#' )
#'
#' @keywords metrics
#' @export
#'
apply_fit <- function(se,
                      fit_fn,
                      data_type,
                      slicing_cols = NULL,
                      slicing_values = NULL,
                      input_assay = NULL,
                      output_assay,
                      summary_fn = NULL,
                      summary_assay = NULL,
                      merge = "merge",
                      on_error = "warn",
                      fit_source) {

  checkmate::assert_string(data_type)
  profile <- get_fit_profile(data_type)  # errors with list of valid profiles if unknown

  if (is.null(slicing_cols)) slicing_cols <- profile$slicing_cols
  if (is.null(slicing_values)) slicing_values <- profile$slicing_values
  if (is.null(input_assay)) input_assay <- profile$input_assay

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_function(fit_fn)
  checkmate::assert_character(slicing_cols, min.len = 1L)
  checkmate::assert_character(slicing_values, null.ok = TRUE, min.len = 1L)
  checkmate::assert_string(input_assay)
  checkmate::assert_string(output_assay)
  checkmate::assert_choice(merge, c("merge", "replace"))
  checkmate::assert_choice(on_error, c("warn", "stop"))
  checkmate::assert_string(fit_source)
  if (!is.null(summary_fn)) {
    checkmate::assert_function(summary_fn)
    checkmate::assert_string(summary_assay)
  }

  if (!input_assay %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("%s assay is required", input_assay))
  }

  tmp_assay <- "__custom_fit_tmp__"
  wrapper_fn <- .make_fit_wrapper(
    fit_fn, slicing_cols, slicing_values, on_error, fit_source
  )

  se_out <- tryCatch(
    gDRutils::apply_bumpy_function(
      se             = se,
      FUN            = wrapper_fn,
      req_assay_name = input_assay,
      out_assay_name = tmp_assay
    ),
    error = function(e) {
      # apply_bumpy_function errors on SE dimnames mismatch when all cells
      # return empty data.tables (nothing to split back into BumpyMatrix).
      # Distinguish from user fit_fn errors propagated with on_error="stop".
      msg <- conditionMessage(e)
      if (grepl("rownames|colnames|dimnames|withDimnames", msg,
                ignore.case = TRUE)) {
        se   # treat as all-empty — return original SE
      } else {
        stop(e)
      }
    }
  )

  if (!tmp_assay %in% SummarizedExperiment::assayNames(se_out)) {
    return(se)
  }

  new_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, tmp_assay),
    row.field = "row", column.field = "column"
  )
  new_metrics <- data.table::as.data.table(new_df)

  if (NROW(new_metrics) == 0L) {
    return(se)
  }

  upsert_key <- c("fit_source", slicing_cols)
  se <- .persist_assay(se, new_metrics, merge, output_assay, "row", "column", upsert_key)

  if (!is.null(summary_fn)) {
    se <- .apply_summary_fn(
      se, new_metrics, summary_fn, summary_assay, merge, fit_source, on_error
    )
  }

  se
}


#' apply_fit_to_se
#'
#' Apply a user-supplied fit function per (drug x cell line x normalization_type)
#' triplet and persist results into the Metrics assay.
#'
#' This is a convenience wrapper around \code{\link{apply_fit}} for the
#' single-agent use case with the standard Metrics assay.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   an Averaged assay
#' @param fit_fn function taking a \code{data.table} and returning a named list
#' @param normalization_types character vector of types to iterate
#' @param averaged_assay name of the input assay
#' @param metrics_assay name of the output assay
#' @param merge \code{"merge"} (idempotent upsert) or \code{"replace"}
#' @param on_error \code{"warn"} (skip + warning) or \code{"stop"}
#' @param fit_source string recorded in the \code{fit_source} column
#'
#' @return updated \code{SummarizedExperiment}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' simple_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
#' se_out <- apply_fit_to_se(se, simple_fn, fit_source = "demo")
#'
#' @keywords metrics
#' @export
#'
apply_fit_to_se <- function(se,
                            fit_fn,
                            normalization_types = c("GR", "RV"),
                            averaged_assay = "Averaged",
                            metrics_assay = "Metrics",
                            merge = "merge",
                            on_error = "warn",
                            fit_source = "custom") {

  # Warn for missing standard metric columns before delegating
  # (apply_fit does not impose column expectations on custom assays)
  result <- apply_fit(
    se             = se,
    fit_fn         = fit_fn,
    data_type      = "single-agent",
    slicing_values = normalization_types,
    input_assay    = averaged_assay,
    output_assay   = metrics_assay,
    merge          = merge,
    on_error       = on_error,
    fit_source     = fit_source
  )

  if (metrics_assay %in% SummarizedExperiment::assayNames(result)) {
    new_df <- BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(result, metrics_assay),
      row.field = "row", column.field = "column"
    )
    missing <- setdiff(
      gDRutils::get_header("response_metrics"),
      names(new_df)
    )
    if (length(missing) > 0L) {
      warning(sprintf(
        "fit_fn output missing columns: %s",
        paste(missing, collapse = ", ")
      ))
    }
  }

  result
}


####
# Multi-fit (single-pass)
####

#' apply_fits
#'
#' Apply multiple fit functions in a \strong{single pass} over each BumpyMatrix
#' cell, writing each function's results into its own named output assay.
#'
#' @details
#' \code{apply_fits()} is the performance-efficient alternative to
#' chaining multiple \code{\link{apply_fit}} calls when two or more fit
#' functions operate on the \strong{same input assay}.  Instead of unsplitting
#' the BumpyMatrix K times (once per function), it traverses each cell once and
#' applies all functions in that single pass.
#'
#' Use this when:
#' \itemize{
#'   \item You have two or more independent fit functions on the same data
#'     (e.g. Bliss + HSS on combination data).
#'   \item A single fit function produces results for multiple output assays
#'     (shared pre-computation pattern — see below).
#' }
#'
#' \subsection{Independent fit functions (most common)}{
#' Names of \code{fit_fns} become the output assay names:
#' \preformatted{
#'   apply_fits(
#'     combo_se,
#'     fit_fns = list(
#'       custom_bliss = bliss_fit_fn,
#'       custom_hss   = hss_fit_fn
#'     ),
#'     data_type  = "combination",
#'     fit_source = "synergy"
#'   )
#' }
#' }
#'
#' \subsection{Shared pre-computation (advanced)}{
#' A single function can return a \strong{named list of named lists} to write
#' multiple assays while computing expensive intermediates only once:
#' \preformatted{
#'   bliss_and_hss <- function(dt) {
#'     sa_curves <- fit_hill_curves(dt)  # expensive — done once
#'     list(
#'       custom_bliss = list(bliss_score = compute_bliss(sa_curves, dt)),
#'       custom_hss   = list(hss_score   = compute_hss(sa_curves, dt))
#'     )
#'   }
#'   apply_fits(
#'     combo_se,
#'     fit_fns    = list(bliss_and_hss = bliss_and_hss),
#'     output_assay_map = c(bliss_and_hss = NA),  # ignored; keys from return value
#'     ...
#'   )
#' }
#' The multi-output pattern is detected automatically when a fit function
#' returns a named list whose values are themselves named lists.  Each top-level
#' name maps to an assay; the inner named list provides the row columns.
#' }
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @param fit_fns named list of functions.  Each name becomes an output assay
#'   name; each value is a function(\code{data.table}) → named list (or a named
#'   list of named lists for the shared pre-computation pattern).
#' @param data_type one of \code{"single-agent"}, \code{"combination"},
#'   \code{"time-course"}.
#' @param slicing_cols character vector; \code{NULL} uses profile default.
#' @param slicing_values character vector; \code{NULL} uses profile default.
#' @param input_assay string; \code{NULL} uses profile default.
#' @param merge \code{"merge"} or \code{"replace"}.
#' @param on_error \code{"warn"} or \code{"stop"}.
#' @param fit_source character string stamped as \code{fit_source} in every
#'   output row.
#'
#' @return updated \code{SummarizedExperiment} with one new (or updated) assay
#'   per entry in \code{fit_fns}.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
#' se <- mae[["single-agent"]]
#' fn_a <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
#' fn_b <- function(dt) list(x_sd = sd(dt$x, na.rm = TRUE))
#' se_out <- apply_fits(
#'   se,
#'   fit_fns = list(mean_metrics = fn_a, sd_metrics = fn_b),
#'   data_type = "single-agent",
#'   fit_source = "demo"
#' )
#'
#' @keywords metrics
#' @export
#'
apply_fits <- function(se,
                       fit_fns,
                       data_type,
                       slicing_cols = NULL,
                       slicing_values = NULL,
                       input_assay = NULL,
                       merge = "merge",
                       on_error = "warn",
                       fit_source) {

  checkmate::assert_string(data_type)
  profile <- get_fit_profile(data_type)

  if (is.null(slicing_cols)) slicing_cols <- profile$slicing_cols
  if (is.null(slicing_values)) slicing_values <- profile$slicing_values
  if (is.null(input_assay)) input_assay <- profile$input_assay

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_list(fit_fns, min.len = 1L, types = "function", names = "unique")
  checkmate::assert_character(slicing_cols, min.len = 1L)
  checkmate::assert_character(slicing_values, null.ok = TRUE, min.len = 1L)
  checkmate::assert_string(input_assay)
  checkmate::assert_choice(merge, c("merge", "replace"))
  checkmate::assert_choice(on_error, c("warn", "stop"))
  checkmate::assert_string(fit_source)

  if (!input_assay %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("%s assay is required", input_assay))
  }

  out_assay_names <- names(fit_fns)
  slice_col <- slicing_cols[1L]

  # Unsplit the BumpyMatrix ONCE into a flat data.table, then split by
  # (row × column) — avoids repeated S4 dispatch per cell that dominates
  # wall-clock time when fit_fn is fast (e.g. Bliss, HSS, summary stats).
  asy <- SummarizedExperiment::assay(se, input_assay)
  dt_all <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(asy, row.field = "row", column.field = "column")
  )
  cell_list <- split(dt_all, by = c("row", "column"), sorted = FALSE)

  # Accumulate per-assay results: list(assay_name -> list of row data.tables)
  assay_results <- stats::setNames(
    lapply(out_assay_names, function(nm) list()),
    out_assay_names
  )

  for (cell_dt in cell_list) {
    r <- cell_dt[["row"]][1L]
    cc <- cell_dt[["column"]][1L]

    vals <- if (!is.null(slicing_values)) {
      slicing_values
    } else if (slice_col %in% names(cell_dt)) {
      unique(cell_dt[[slice_col]])
    } else {
      NA_character_
    }

    for (val in vals) {
      sub_dt <- if (slice_col %in% names(cell_dt) && !is.na(val)) {
        cell_dt[cell_dt[[slice_col]] == val, ]
      } else {
        cell_dt
      }
      if (NROW(sub_dt) == 0L) next
      label <- if (!is.na(val)) sprintf("%s=%s", slice_col, val) else "all"

      for (fn_nm in out_assay_names) {
        fn <- fit_fns[[fn_nm]]
        result_dt <- tryCatch({
          res <- fn(sub_dt)
          # Detect multi-output pattern: named list of named lists
          if (is.list(res) && length(res) > 0L &&
              is.list(res[[1L]]) && !is.null(names(res))) {
            for (assay_nm in names(res)) {
              row_dt <- data.table::as.data.table(as.list(res[[assay_nm]]))
              row_dt[["fit_source"]] <- fit_source
              row_dt[["row"]] <- r
              row_dt[["column"]] <- cc
              if (!is.na(val) && !slice_col %in% names(row_dt)) {
                row_dt[[slice_col]] <- val
              }
              if (assay_nm %in% out_assay_names) {
                assay_results[[assay_nm]] <- c(
                  assay_results[[assay_nm]], list(row_dt)
                )
              }
            }
            next
          }
          data.table::as.data.table(as.list(res))
        }, error = function(e) {
          if (on_error == "stop") stop(e)
          warning(sprintf("fit_fn '%s' failed for %s row=%s col=%s: %s",
                          fn_nm, label, r, cc, conditionMessage(e)))
          NULL
        })

        if (is.null(result_dt)) next

        result_dt[["fit_source"]] <- fit_source
        result_dt[["row"]] <- r
        result_dt[["column"]] <- cc
        if (!is.na(val) && !slice_col %in% names(result_dt)) {
          result_dt[[slice_col]] <- val
        }
        assay_results[[fn_nm]] <- c(assay_results[[fn_nm]], list(result_dt))
      }
    }
  }

  # Persist each assay's accumulated rows
  upsert_key <- c("fit_source", slicing_cols)
  for (assay_nm in out_assay_names) {
    rows <- assay_results[[assay_nm]]
    if (length(rows) == 0L) next
    merged_dt <- data.table::rbindlist(rows, fill = TRUE)
    se <- .persist_assay(se, merged_dt, merge, assay_nm, "row", "column", upsert_key)
  }

  se
}


####
# Reference fit functions
####

#' fit_drug_response_metrics
#'
#' Reference fit function replicating standard gDR Hill curve fitting
#' (3-parameter model, \code{x_0} fixed at 1).  Produces results comparable
#' to \code{\link[gDRcore]{fit_SE}} and \code{\link[gDRutils]{logisticFit}}.
#' For use with \code{\link{apply_fit_to_se}} or \code{\link{apply_fit}}
#' on single-agent data.
#'
#' @param avg_dt \code{data.table} of averaged data for one
#'   (drug x cell line x normalization_type) triplet
#' @param capping_fold numeric capping fold passed to
#'   \code{\link[gDRutils]{cap_xc50}}
#'
#' @return Named list of fit metrics compatible with the standard gDR
#'   \code{Metrics} assay schema.
#'
#' @seealso \code{\link{fit_drug_response_metrics_4p}} for a 4-parameter
#'   variant that fits \code{x_0} freely.
#'
#' @examples
#' dt <- data.table::data.table(
#'   Concentration = c(0.001, 0.01, 0.1, 1, 10),
#'   x = c(0.95, 0.8, 0.5, 0.2, 0.1),
#'   normalization_type = "RV"
#' )
#' fit_drug_response_metrics(dt)
#'
#' @export
fit_drug_response_metrics <- function(avg_dt, capping_fold = 5) {
  .fit_drug_response_metrics_impl(avg_dt, x_0 = 1, capping_fold = capping_fold)
}


#' fit_drug_response_metrics_4p
#'
#' Variant of \code{\link{fit_drug_response_metrics}} using a 4-parameter
#' Hill model — \code{x_0} is fitted freely rather than fixed at 1.
#' Useful when the upper asymptote is expected to deviate from 1 (e.g.
#' partial agonists, non-standard assay readouts).
#'
#' @inheritParams fit_drug_response_metrics
#'
#' @return Named list of fit metrics.  \code{fit_type} will be
#'   \code{"DRC4pHillFitModel"} on success.
#'
#' @seealso \code{\link{fit_drug_response_metrics}} for the default
#'   3-parameter variant that matches \code{fit_SE()} output.
#'
#' @examples
#' dt <- data.table::data.table(
#'   Concentration = c(0.001, 0.01, 0.1, 1, 10),
#'   x = c(0.95, 0.8, 0.5, 0.2, 0.1),
#'   normalization_type = "RV"
#' )
#' fit_drug_response_metrics_4p(dt)
#'
#' @export
fit_drug_response_metrics_4p <- function(avg_dt, capping_fold = 5) {
  .fit_drug_response_metrics_impl(avg_dt, x_0 = NA_real_, capping_fold = capping_fold)
}


#' @keywords internal
.fit_drug_response_metrics_impl <- function(avg_dt, x_0 = 1, capping_fold = 5) {
  norm_type <- avg_dt$normalization_type[1]
  conc_col <- gDRutils::get_env_identifiers("concentration")
  conc <- avg_dt[[conc_col]]
  x <- avg_dt$x

  keep <- !is.na(x) & !is.na(conc)
  x <- x[keep]
  conc <- conc[keep]

  if (length(x) == 0L) {
    return(.empty_fit_result(norm_type))
  }

  x_mean_obs <- mean(x, na.rm = TRUE)
  N_conc <- length(unique(conc))
  maxlog10Conc <- log10(max(conc, na.rm = TRUE))

  # Parameters matching logisticFit() in gDRutils (fit_curves.R:113-131):
  #   RV: priors c(2, 0.4, 1, med), lower c(0.1,  0,  0, min/10)
  #   GR: priors c(2, 0.1, 1, med), lower c(0.1, -1, -1, min/10)
  x_inf_prior <- if (norm_type == "GR") 0.1 else 0.4
  lower_x_inf <- if (norm_type == "GR") -1 else 0
  lower_x_0   <- if (norm_type == "GR") -1 else 0
  med_conc <- stats::median(conc)
  min_conc <- min(conc)
  max_conc <- max(conc)
  controls <- drc::drmc(relTol = 1e-04, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)
  cap <- 0.1  # same default as fit_SE / logisticFit

  # Cap x values as logisticFit does (pmin to x_0 + cap)
  x_cap_limit <- if (is.na(x_0)) 1 + cap else x_0 + cap
  x <- pmin(x, x_cap_limit)

  four_param <- is.na(x_0)
  fit_param <- c("h", "x_inf", "x_0", "ec50")

  if (four_param) {
    fct <- drc::LL.4(names = fit_param)
    start  <- c(2, x_inf_prior, 1, med_conc)
    lowerl <- c(0.1, lower_x_inf, lower_x_0, min_conc / 10)
    upperl <- c(5, 1, 1 + cap, max_conc * 10)
  } else {
    fit_param_3p <- fit_param[-3]   # drop x_0
    fct <- drc::LL.3u(upper = x_0, names = fit_param_3p)
    start  <- c(2, x_inf_prior, med_conc)
    lowerl <- c(0.1, lower_x_inf, min_conc / 10)
    upperl <- c(5, min(x_0 + cap, 1), max_conc * 10)
  }

  fit <- tryCatch(
    drc::drm(x ~ conc,
             data = data.table::data.table(x = x, conc = conc),
             logDose = NULL,
             fct = fct,
             start = start,
             lowerl = lowerl,
             upperl = upperl,
             control = controls,
             na.action = stats::na.omit),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    coefs <- stats::coef(fit)
    rss <- sum(stats::residuals(fit)^2, na.rm = TRUE)
    r2 <- 1 - rss / sum((x - x_mean_obs)^2)
    # x_mean from model prediction (matches fit_SE / logisticFit behaviour)
    x_mean_model <- mean(
      stats::predict(fit, data.table::data.table(
        conc = 10^seq(log10(min(conc)), log10(max(conc)), length.out = 101L)
      )),
      na.rm = TRUE
    )
    if (four_param) {
      ec50_val <- coefs[["ec50:(Intercept)"]]
      x_0_val <- coefs[["x_0:(Intercept)"]]
      x_inf_val <- coefs[["x_inf:(Intercept)"]]
      fit_type <- "DRC4pHillFitModel"
    } else {
      ec50_val <- coefs[["ec50:(Intercept)"]]
      x_0_val <- x_0
      x_inf_val <- coefs[["x_inf:(Intercept)"]]
      fit_type <- "DRC3pHillFitModelFixS0"
    }
    list(
      normalization_type = norm_type,
      x_mean = x_mean_model,
      x_AOC = 1 - x_mean_model,
      N_conc = N_conc,
      maxlog10Concentration = maxlog10Conc,
      ec50 = ec50_val,
      xc50 = gDRutils::cap_xc50(ec50_val, max(conc), capping_fold = capping_fold),
      h = coefs[1],
      r2 = r2,
      x_0 = x_0_val,
      x_inf = x_inf_val,
      fit_type = fit_type
    )
  } else {
    list(
      normalization_type = norm_type,
      x_mean = x_mean_obs,
      x_AOC = 1 - x_mean_obs,
      N_conc = N_conc,
      maxlog10Concentration = maxlog10Conc,
      ec50 = NA_real_,
      xc50 = .estimate_xc50_fallback(x),
      h = NA_real_,
      r2 = NA_real_,
      x_0 = if (four_param) NA_real_ else x_0,
      x_inf = NA_real_,
      fit_type = "DRCInvalidFitResult"
    )
  }
}


####
# apply_combo_scores — high-level wrapper replicating fit_SE.combinations scoring
####

#' apply_combo_scores
#'
#' Compute Bliss and HSA synergy scores for a combination SE, replicating the
#' exact scoring logic of \code{\link{fit_SE.combinations}}.
#'
#' Unlike the low-level \code{\link{bliss_fit_fn}} and \code{\link{hss_fit_fn}}
#' (which work on raw Averaged data per triplet), this function uses the SA fit
#' parameters from the \code{Metrics} assay to generate curve-smoothed single-agent
#' responses via \code{\link[gDRutils]{predict_efficacy_from_conc}} before computing
#' excess.  This produces results numerically identical to \code{fit_SE.combinations}.
#'
#' @details
#' The function requires a \code{Metrics} assay containing columns
#' \code{dilution_drug}, \code{cotrt_value}, \code{ec50}, \code{h},
#' \code{x_inf}, \code{x_0}, and \code{normalization_type} — as produced by
#' \code{fit_SE.combinations} or \code{\link{apply_fit_to_se}} with
#' \code{fit_drug_response_metrics}.
#'
#' Scoring steps (per drug-combo × cell-line × normalization_type):
#' \enumerate{
#'   \item Predict smooth SA responses at every combo concentration using
#'     \code{predict_efficacy_from_conc} on \code{drug_1} and \code{drug_2}
#'     parameter sets from \code{Metrics}.
#'   \item Average \code{col_values} (drug-1-along-conc1) and
#'     \code{row_values} (drug-2-along-conc2) → \code{smooth}.
#'   \item Compute Bliss-expected using \code{\link{calculate_Bliss}} and
#'     HSA-expected using \code{\link{calculate_HSA}} on the smooth SA edges.
#'   \item Compute excess via \code{\link{calculate_excess}}.
#'   \item Score = mean of top 10-percentile excess via \code{\link{calculate_score}}.
#' }
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   \code{"Averaged"} and \code{"Metrics"} assays (combination experiment).
#' @param scores_assay string; name of the output assay to write scores into.
#'   Default \code{"scores"}.
#' @param averaged_assay string; name of the input assay. Default \code{"Averaged"}.
#' @param metrics_assay string; name of the assay containing SA fit parameters.
#'   Default \code{"Metrics"}.
#' @param normalization_types character vector of normalization types to process.
#'   Default \code{c("GR", "RV")}.
#' @param fit_source string recorded in the \code{fit_source} column. Default
#'   \code{"gDR"}.
#'
#' @return Updated \code{SummarizedExperiment} with a \code{scores_assay}
#'   assay containing \code{bliss_score} and \code{hsa_score} per triplet.
#'
#' @seealso \code{\link{bliss_fit_fn}}, \code{\link{hss_fit_fn}} for the
#'   simplified raw-data variants.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_out <- apply_combo_scores(combo_se)
#' "scores" %in% SummarizedExperiment::assayNames(combo_se_out)
#'
#' @importFrom checkmate assert_class assert_string assert_character
#' @export
apply_combo_scores <- function(se,
                               scores_assay = "scores",
                               averaged_assay = "Averaged",
                               metrics_assay = "Metrics",
                               normalization_types = c("GR", "RV"),
                               fit_source = "gDR") {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(scores_assay)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_character(normalization_types, min.len = 1L)
  gDRutils::validate_se_assay_name(se, averaged_assay)
  gDRutils::validate_se_assay_name(se, metrics_assay)

  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")
  series_identifiers <- c(conc1_col, conc2_col)

  avg_bm <- SummarizedExperiment::assay(se, averaged_assay)
  met_bm <- SummarizedExperiment::assay(se, metrics_assay)

  # Iterate over each (drug-combo × cell-line) BumpyMatrix cell
  score_rows <- vector("list", prod(dim(se)) * length(normalization_types))
  idx <- 0L

  for (ri in seq_len(nrow(se))) {
    for (ci in seq_len(ncol(se))) {
      avg_dt <- as.data.table(avg_bm[ri, ci][[1L]])
      met_dt <- as.data.table(met_bm[ri, ci][[1L]])

      if (NROW(avg_dt) == 0L || NROW(met_dt) == 0L) next

      row_nm <- rownames(se)[ri]
      col_nm <- colnames(se)[ci]

      for (norm_type in normalization_types) {
        avg_sub <- avg_dt[avg_dt$normalization_type == norm_type, ]
        met_sub <- met_dt[met_dt$normalization_type == norm_type, ]
        if (NROW(avg_sub) == 0L || NROW(met_sub) == 0L) next

        # Smooth SA predictions for every point in the averaged combo data
        # drug_1 params: SA along conc1 for each cotrt conc2 value
        drug1_params <- met_sub[met_sub$dilution_drug == "drug_1", ]
        drug2_params <- met_sub[met_sub$dilution_drug == "drug_2", ]

        if (NROW(drug1_params) == 0L || NROW(drug2_params) == 0L) next

        # col_values: predict drug-1 response at conc1 for each cotrt conc2
        avg_sub$col_values <- vapply(seq_len(NROW(avg_sub)), function(i) {
          cotrt <- avg_sub[[conc2_col]][i]
          conc <- avg_sub[[conc1_col]][i]
          p <- drug1_params[abs(drug1_params$cotrt_value - cotrt) ==
                              min(abs(drug1_params$cotrt_value - cotrt)), ][1L, ]
          if (NROW(p) == 0L || any(is.na(c(p$ec50, p$h, p$x_inf, p$x_0)))) return(NA_real_)
          gDRutils::predict_efficacy_from_conc(conc, p$x_inf, p$x_0, p$ec50, p$h)
        }, numeric(1))

        # row_values: predict drug-2 response at conc2 for each cotrt conc1
        avg_sub$row_values <- vapply(seq_len(NROW(avg_sub)), function(i) {
          cotrt <- avg_sub[[conc1_col]][i]
          conc <- avg_sub[[conc2_col]][i]
          p <- drug2_params[abs(drug2_params$cotrt_value - cotrt) ==
                              min(abs(drug2_params$cotrt_value - cotrt)), ][1L, ]
          if (NROW(p) == 0L || any(is.na(c(p$ec50, p$h, p$x_inf, p$x_0)))) return(NA_real_)
          gDRutils::predict_efficacy_from_conc(conc, p$x_inf, p$x_0, p$ec50, p$h)
        }, numeric(1))

        avg_sub$smooth <- rowMeans(
          avg_sub[, c("col_values", "row_values"), with = FALSE],
          na.rm = TRUE
        )
        avg_sub[avg_sub[[conc1_col]] == 0 & avg_sub[[conc2_col]] == 0, smooth := 1]

        # SA edges for calculate_HSA / calculate_Bliss
        sa1 <- avg_sub[avg_sub[[conc2_col]] == 0, c(series_identifiers, "smooth"), with = FALSE]
        sa2 <- avg_sub[avg_sub[[conc1_col]] == 0, c(series_identifiers, "smooth"), with = FALSE]

        if (NROW(sa1) == 0L || NROW(sa2) == 0L) next

        # HSA score
        hsa_expected <- calculate_HSA(sa1, conc1_col, sa2, conc2_col, norm_type)
        h_excess <- calculate_excess(
          hsa_expected, avg_sub,
          series_identifiers = series_identifiers,
          metric_col = "metric", measured_col = "smooth"
        )
        hsa_score <- if (all(is.na(h_excess$x))) NA_real_ else
          calculate_score(h_excess$x)

        # Bliss score
        bliss_expected <- calculate_Bliss(sa1, conc1_col, sa2, conc2_col, norm_type)
        bliss_excess <- calculate_excess(
          bliss_expected, avg_sub,
          series_identifiers = series_identifiers,
          metric_col = "metric", measured_col = "smooth"
        )
        bliss_score <- if (all(is.na(bliss_excess$x))) NA_real_ else
          calculate_score(bliss_excess$x)

        idx <- idx + 1L
        score_rows[[idx]] <- data.table::data.table(
          row = row_nm,
          column = col_nm,
          normalization_type = norm_type,
          fit_source = fit_source,
          bliss_score = bliss_score,
          hsa_score = hsa_score
        )
      }
    }
  }

  if (idx == 0L) {
    return(se)
  }

  scores_dt <- data.table::rbindlist(score_rows[seq_len(idx)])
  se <- .persist_assay(se, scores_dt, "replace", scores_assay, "row", "column",
                       upsert_key = c("fit_source", "normalization_type"))
  se
}


#' bliss_fit_fn
#'
#' Reference fit function for Bliss independence synergy scoring on combination
#' data.  Computes Bliss-expected response from single-agent edges and derives
#' excess and score from raw averaged data.
#'
#' Intended for use with \code{\link{apply_fit}} on combination SEs:
#'
#' \preformatted{
#'   apply_fit(
#'     combo_se, bliss_fit_fn, "combination",
#'     output_assay = "custom_bliss", fit_source = "bliss"
#'   )
#' }
#'
#' @param avg_dt \code{data.table} for one (drug1 x drug2 x cell line x
#'   normalization_type) combination.  Must contain columns \code{Concentration},
#'   \code{Concentration_2}, \code{x}, and \code{normalization_type}.
#'
#' @return Named list with \code{bliss_score}, \code{bliss_excess_mean},
#'   \code{n_combo_points}, and \code{normalization_type}.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_out <- apply_fit(
#'   combo_se, bliss_fit_fn, "combination",
#'   output_assay = "custom_bliss", fit_source = "bliss"
#' )
#'
#' @export
bliss_fit_fn <- function(avg_dt) {
  norm_type <- avg_dt$normalization_type[1]
  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")

  sa1 <- avg_dt[avg_dt[[conc2_col]] == 0 & avg_dt[[conc1_col]] > 0, ]
  sa2 <- avg_dt[avg_dt[[conc1_col]] == 0 & avg_dt[[conc2_col]] > 0, ]
  combo <- avg_dt[avg_dt[[conc1_col]] > 0 & avg_dt[[conc2_col]] > 0, ]

  n_combo <- NROW(combo)
  if (n_combo == 0L || NROW(sa1) == 0L || NROW(sa2) == 0L) {
    return(list(
      normalization_type = norm_type,
      bliss_score        = NA_real_,
      bliss_excess_mean  = NA_real_,
      n_combo_points     = n_combo
    ))
  }

  # Match SA responses to each combo point by concentration (per-point comparison)
  combo_ordered <- combo[order(combo[[conc1_col]], combo[[conc2_col]]), ]
  sa1_match <- sa1[["x"]][match(combo_ordered[[conc1_col]], sa1[[conc1_col]])]
  sa2_match <- sa2[["x"]][match(combo_ordered[[conc2_col]], sa2[[conc2_col]])]

  expected <- if (norm_type == "RV") {
    sa1_match * sa2_match
  } else {
    # GR adaptation: Holbeck et al., Cancer Res 2017
    2^(log2(sa1_match + 1) * log2(sa2_match + 1)) - 1
  }

  excess <- expected - combo_ordered[["x"]]
  q90 <- stats::quantile(excess, 0.9, na.rm = TRUE)

  list(
    normalization_type = norm_type,
    bliss_score = mean(excess[excess >= q90], na.rm = TRUE),
    bliss_excess_mean = mean(excess, na.rm = TRUE),
    n_combo_points = n_combo
  )
}


#' hss_fit_fn
#'
#' Reference fit function for Highest Single Agent (HSA) synergy scoring on
#' combination data.  Uses raw averaged single-agent edges.
#'
#' Intended for use with \code{\link{apply_fit}} on combination SEs:
#'
#' \preformatted{
#'   apply_fit(
#'     combo_se, hss_fit_fn, "combination",
#'     output_assay = "custom_hss", fit_source = "hss"
#'   )
#' }
#'
#' @param avg_dt \code{data.table} for one (drug1 x drug2 x cell line x
#'   normalization_type) combination.  Must contain columns \code{Concentration},
#'   \code{Concentration_2}, \code{x}, and \code{normalization_type}.
#'
#' @return Named list with \code{hss_score}, \code{hss_excess_mean},
#'   \code{n_combo_points}, and \code{normalization_type}.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_out <- apply_fit(
#'   combo_se, hss_fit_fn, "combination",
#'   output_assay = "custom_hss", fit_source = "hss"
#' )
#'
#' @export
hss_fit_fn <- function(avg_dt) {
  norm_type <- avg_dt$normalization_type[1]
  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")

  sa1 <- avg_dt[avg_dt[[conc2_col]] == 0 & avg_dt[[conc1_col]] > 0, ]
  sa2 <- avg_dt[avg_dt[[conc1_col]] == 0 & avg_dt[[conc2_col]] > 0, ]
  combo <- avg_dt[avg_dt[[conc1_col]] > 0 & avg_dt[[conc2_col]] > 0, ]

  n_combo <- NROW(combo)
  if (n_combo == 0L || NROW(sa1) == 0L || NROW(sa2) == 0L) {
    return(list(
      normalization_type = norm_type,
      hss_score = NA_real_,
      hss_excess_mean = NA_real_,
      n_combo_points = n_combo
    ))
  }

  # Match SA responses to each combo point by concentration (per-point comparison)
  combo_ordered <- combo[order(combo[[conc1_col]], combo[[conc2_col]]), ]
  sa1_match <- sa1[["x"]][match(combo_ordered[[conc1_col]], sa1[[conc1_col]])]
  sa2_match <- sa2[["x"]][match(combo_ordered[[conc2_col]], sa2[[conc2_col]])]

  # HSA expected: the more potent single agent at each combo grid point
  expected <- pmin(sa1_match, sa2_match)
  excess <- expected - combo_ordered[["x"]]
  q90 <- stats::quantile(excess, 0.9, na.rm = TRUE)

  list(
    normalization_type = norm_type,
    hss_score = mean(excess[excess >= q90], na.rm = TRUE),
    hss_excess_mean = mean(excess, na.rm = TRUE),
    n_combo_points = n_combo
  )
}


####
# Internal helpers
####

#' Build a per-cell wrapper for apply_bumpy_function
#'
#' Iterates over each unique combination of slicing_col values, filters the
#' cell data.table, calls fit_fn, stamps fit_source and slicing column values,
#' and rbinds all results.  Respects on_error semantics.
#'
#' @param fit_fn user-supplied fit function
#' @param slicing_cols character vector of column names to slice by
#' @param slicing_values character vector of values to iterate, or NULL for all
#' @param on_error \code{"warn"} or \code{"stop"}
#' @param fit_source fit_source label
#'
#' @return function(DFrame) -> data.table
#'
#' @keywords internal
.make_fit_wrapper <- function(fit_fn, slicing_cols, slicing_values,
                              on_error, fit_source) {
  function(avg_df) {
    avg_dt <- data.table::as.data.table(avg_df)
    results <- list()

    # Determine the actual values to iterate for the primary slicing column.
    # Multi-column slicing is deferred to a future iteration; for the POC the
    # profiles always supply a single slicing column.
    slice_col <- slicing_cols[1L]

    vals <- if (!is.null(slicing_values)) {
      slicing_values
    } else if (slice_col %in% names(avg_dt)) {
      unique(avg_dt[[slice_col]])
    } else {
      NA_character_   # no slicing column present → one call on whole cell
    }

    for (val in vals) {
      sub_dt <- if (slice_col %in% names(avg_dt) && !is.na(val)) {
        avg_dt[avg_dt[[slice_col]] == val, ]
      } else {
        avg_dt
      }

      if (NROW(sub_dt) == 0L) next

      label <- if (!is.na(val)) sprintf("%s=%s", slice_col, val) else "all"
      result_dt <- tryCatch({
        res <- fit_fn(sub_dt)
        data.table::as.data.table(as.list(res))
      }, error = function(e) {
        if (on_error == "stop") stop(e)
        warning(sprintf("fit_fn failed for %s: %s", label, conditionMessage(e)))
        NULL
      })

      if (is.null(result_dt)) next

      # Stamp mandatory identifier columns (caller's value always wins)
      result_dt[["fit_source"]] <- fit_source
      if (!is.na(val) && !slice_col %in% names(result_dt)) {
        result_dt[[slice_col]] <- val
      }

      results[[length(results) + 1L]] <- result_dt
    }

    if (length(results) == 0L) return(data.table::data.table())
    data.table::rbindlist(results, fill = TRUE)
  }
}


#' Apply summary_fn across cells and persist into summary_assay
#'
#' @param se SummarizedExperiment
#' @param fit_dt flat data.table with "row", "column", and all fit columns
#' @param summary_fn function(data.table) -> named list
#' @param summary_assay assay name for summary output
#' @param merge merge strategy
#' @param fit_source fit_source label
#' @param on_error error handling strategy
#'
#' @return updated SummarizedExperiment
#'
#' @keywords internal
.apply_summary_fn <- function(se, fit_dt, summary_fn, summary_assay,
                              merge, fit_source, on_error) {
  cell_list <- split(fit_dt, by = c("row", "column"), sorted = FALSE)

  summary_rows <- lapply(cell_list, function(cell_fit) {
    r <- cell_fit[["row"]][1L]
    cc <- cell_fit[["column"]][1L]

    result <- tryCatch({
      res <- summary_fn(cell_fit)
      data.table::as.data.table(as.list(res))
    }, error = function(e) {
      if (on_error == "stop") stop(e)
      warning(sprintf(
        "summary_fn failed for row=%s col=%s: %s", r, cc, conditionMessage(e)
      ))
      NULL
    })

    if (is.null(result)) return(NULL)
    result[["row"]] <- r
    result[["column"]] <- cc
    result[["fit_source"]] <- fit_source
    result
  })

  summary_rows <- Filter(Negate(is.null), summary_rows)
  if (length(summary_rows) == 0L) return(se)

  summary_dt <- data.table::rbindlist(summary_rows, fill = TRUE)
  .persist_assay(se, summary_dt, merge, summary_assay, "row", "column", "fit_source")
}


#' Persist a data.table into a BumpyMatrix assay
#'
#' Performs an idempotent upsert (merge mode) keyed by \code{upsert_key_cols},
#' or a full overwrite (replace mode).
#'
#' @param se SummarizedExperiment
#' @param new_dt data.table including "row" and "column" index fields
#' @param merge \code{"merge"} or \code{"replace"}
#' @param assay_name assay name to write into
#' @param row row index column name in new_dt
#' @param col column index column name in new_dt
#' @param upsert_key_cols character vector of column(s) forming the upsert key
#'
#' @return updated SummarizedExperiment
#'
#' @keywords internal
.persist_assay <- function(se, new_dt, merge, assay_name, row, col,
                           upsert_key_cols) {
  has_assay <- assay_name %in% SummarizedExperiment::assayNames(se)

  if (merge == "merge" && has_assay) {
    existing_df <- BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se, assay_name),
      row.field = row, column.field = col
    )
    existing_dt <- data.table::as.data.table(existing_df)

    # Ensure all key columns exist in existing data (fill with NA if absent)
    for (kc in upsert_key_cols) {
      if (!kc %in% names(existing_dt)) existing_dt[[kc]] <- NA_character_
    }

    key_cols_present <- intersect(upsert_key_cols, names(new_dt))
    # data.table anti-join: keep existing rows whose key is not in new_dt
    existing_pruned <- existing_dt[!new_dt, on = key_cols_present]

    merged_dt <- data.table::rbindlist(list(existing_pruned, new_dt), fill = TRUE)
  } else {
    merged_dt <- new_dt
  }

  data_cols <- names(merged_dt)[!names(merged_dt) %in% c(row, col)]
  mx <- BumpyMatrix::splitAsBumpyMatrix(
    merged_dt[, data_cols, with = FALSE],
    row = merged_dt[[row]], col = merged_dt[[col]]
  )
  SummarizedExperiment::assay(se, assay_name) <- mx
  se
}


#' @keywords internal
.empty_fit_result <- function(norm_type) {
  list(
    fit_source = "custom",
    normalization_type = norm_type,
    x_mean = NA_real_,
    x_AOC = NA_real_,
    N_conc = 0L,
    maxlog10Concentration = NA_real_,
    xc50 = NA_real_,
    h = NA_real_,
    r2 = NA_real_,
    x_0 = NA_real_,
    x_inf = NA_real_,
    fit_type = "DRCInvalidFitResult"
  )
}


#' @keywords internal
.estimate_xc50_fallback <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0L) {
    NA_real_
  } else if (all(x_clean > 0.5)) {
    Inf
  } else if (all(x_clean <= 0.5)) {
    -Inf
  } else {
    NA_real_
  }
}


####
# Deprecated aliases
####

#' @rdname apply_fit
#' @export
#' @keywords internal
apply_custom_fit <- function(...) {
  .Deprecated("apply_fit")
  apply_fit(...)
}

#' @rdname apply_fits
#' @export
#' @keywords internal
apply_custom_fits <- function(...) {
  .Deprecated("apply_fits")
  apply_fits(...)
}

#' @keywords internal
.persist_metrics <- function(se, new_metrics, merge, metrics_assay, row, col) {
  .persist_assay(
    se, new_metrics, merge, metrics_assay, row, col,
    upsert_key_cols = c("fit_source", "normalization_type")
  )
}
