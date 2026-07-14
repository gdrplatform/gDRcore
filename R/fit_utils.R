####
# Data-type profiles
####

# Default slicing configuration per experiment type.
# slicing_cols:   column(s) in the Averaged assay that define sub-experiment slices.
# slicing_values: default values to iterate; NULL means "all unique values found in data".
# input_assay:    default source assay name.
#
# Override any field by passing explicit arguments to apply_custom_fit().
.CUSTOM_FIT_PROFILES <- list(
  "single-agent" = list(
    slicing_cols   = "normalization_type",
    slicing_values = c("GR", "RV"),
    input_assay    = "Averaged"
  ),
  "combination" = list(
    slicing_cols   = "normalization_type",
    slicing_values = c("GR", "RV"),
    input_assay    = "Averaged"
  ),
  "time-course" = list(
    slicing_cols   = "normalization_type",
    slicing_values = c("GR", "RV"),
    input_assay    = "Averaged"
  )
)


####
# Public API
####

#' apply_custom_fit
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
#'     apply_custom_fit(bliss_fn, "combination",
#'                      output_assay = "custom_bliss", ...) |>
#'     apply_custom_fit(musyc_fn,  "combination",
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
#' se_out <- apply_custom_fit(
#'   se, mean_fn, "single-agent",
#'   output_assay = "custom_mean", fit_source = "demo"
#' )
#'
#' @keywords metrics
#' @export
#'
apply_custom_fit <- function(se,
                             fit_fn,
                             data_type      = c("single-agent",
                                                "combination",
                                                "time-course"),
                             slicing_cols   = NULL,
                             slicing_values = NULL,
                             input_assay    = NULL,
                             output_assay,
                             summary_fn     = NULL,
                             summary_assay  = NULL,
                             merge          = "merge",
                             on_error       = "warn",
                             fit_source) {

  data_type <- match.arg(data_type)
  profile <- .CUSTOM_FIT_PROFILES[[data_type]]

  if (is.null(slicing_cols))   slicing_cols   <- profile$slicing_cols
  if (is.null(slicing_values)) slicing_values <- profile$slicing_values
  if (is.null(input_assay))    input_assay    <- profile$input_assay

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

  se_out <- gDRutils::apply_bumpy_function(
    se             = se,
    FUN            = wrapper_fn,
    req_assay_name = input_assay,
    out_assay_name = tmp_assay
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
#' This is a convenience wrapper around \code{\link{apply_custom_fit}} for the
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
  # (apply_custom_fit does not impose column expectations on custom assays)
  result <- apply_custom_fit(
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
# Reference fit functions
####

#' fit_drug_response_metrics
#'
#' Reference fit function replicating standard gDR Hill curve fitting.
#' For use with \code{\link{apply_fit_to_se}} or \code{\link{apply_custom_fit}}
#' on single-agent data.
#'
#' @param avg_dt \code{data.table} of averaged data for one
#'   (drug x cell line x normalization_type) triplet
#' @param capping_fold numeric capping fold passed to
#'   \code{\link[gDRutils]{cap_xc50}}
#'
#' @return Named list of fit metrics
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

  x_mean <- mean(x, na.rm = TRUE)
  x_AOC <- 1 - x_mean
  N_conc <- length(unique(conc))
  maxlog10Conc <- log10(max(conc, na.rm = TRUE))

  x_inf_prior <- if (norm_type == "GR") 0.1 else 0.4
  controls <- drc::drmc(relTol = 1e-06, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)
  fit <- tryCatch(
    drc::drm(x ~ conc,
             data = data.table::data.table(x = x, conc = conc),
             fct = drc::LL.4(),
             start = c(2, x_inf_prior, 1, stats::median(conc)),
             control = controls),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    coefs <- stats::coef(fit)
    r2 <- 1 - sum(stats::residuals(fit)^2) / sum((x - mean(x))^2)
    list(
      normalization_type    = norm_type,
      x_mean                = x_mean,
      x_AOC                 = x_AOC,
      N_conc                = N_conc,
      maxlog10Concentration = maxlog10Conc,
      xc50    = gDRutils::cap_xc50(coefs[4], max(conc), capping_fold = capping_fold),
      h       = coefs[1],
      r2      = r2,
      x_0     = coefs[3],
      x_inf   = coefs[2],
      fit_type = "DRC4pHillFitModel"
    )
  } else {
    list(
      normalization_type    = norm_type,
      x_mean                = x_mean,
      x_AOC                 = x_AOC,
      N_conc                = N_conc,
      maxlog10Concentration = maxlog10Conc,
      xc50     = .estimate_xc50_fallback(x),
      h        = NA_real_,
      r2       = NA_real_,
      x_0      = NA_real_,
      x_inf    = NA_real_,
      fit_type = "DRCInvalidFitResult"
    )
  }
}


#' bliss_fit_fn
#'
#' Reference fit function for Bliss independence synergy scoring on combination
#' data.  Computes Bliss-expected response from single-agent edges and derives
#' excess and score from raw averaged data.
#'
#' Intended for use with \code{\link{apply_custom_fit}} on combination SEs:
#'
#' \preformatted{
#'   apply_custom_fit(
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
#' @export
bliss_fit_fn <- function(avg_dt) {
  norm_type   <- avg_dt$normalization_type[1]
  conc1_col   <- gDRutils::get_env_identifiers("concentration")
  conc2_col   <- gDRutils::get_env_identifiers("concentration2")

  sa1   <- avg_dt[avg_dt[[conc2_col]] == 0 & avg_dt[[conc1_col]] > 0, ]
  sa2   <- avg_dt[avg_dt[[conc1_col]] == 0 & avg_dt[[conc2_col]] > 0, ]
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

  sa1_x <- sa1[order(sa1[[conc1_col]]), ][["x"]]
  sa2_x <- sa2[order(sa2[[conc2_col]]), ][["x"]]

  expected <- if (norm_type == "RV") {
    as.vector(outer(sa1_x, sa2_x, `*`))
  } else {
    # GR adaptation: Holbeck et al., Cancer Res 2017
    as.vector(outer(log2(sa1_x + 1), log2(sa2_x + 1),
                    function(a, b) 2^(a * b) - 1))
  }

  combo_x <- combo[["x"]]
  excess <- expected - mean(combo_x, na.rm = TRUE)
  q90 <- stats::quantile(excess, 0.9, na.rm = TRUE)

  list(
    normalization_type = norm_type,
    bliss_score        = mean(excess[excess >= q90], na.rm = TRUE),
    bliss_excess_mean  = mean(excess, na.rm = TRUE),
    n_combo_points     = n_combo
  )
}


#' hss_fit_fn
#'
#' Reference fit function for Highest Single Agent (HSA) synergy scoring on
#' combination data.  Uses raw averaged single-agent edges.
#'
#' Intended for use with \code{\link{apply_custom_fit}} on combination SEs:
#'
#' \preformatted{
#'   apply_custom_fit(
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
#' @export
hss_fit_fn <- function(avg_dt) {
  norm_type <- avg_dt$normalization_type[1]
  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")

  sa1   <- avg_dt[avg_dt[[conc2_col]] == 0 & avg_dt[[conc1_col]] > 0, ]
  sa2   <- avg_dt[avg_dt[[conc1_col]] == 0 & avg_dt[[conc2_col]] > 0, ]
  combo <- avg_dt[avg_dt[[conc1_col]] > 0 & avg_dt[[conc2_col]] > 0, ]

  n_combo <- NROW(combo)
  if (n_combo == 0L || NROW(sa1) == 0L || NROW(sa2) == 0L) {
    return(list(
      normalization_type = norm_type,
      hss_score          = NA_real_,
      hss_excess_mean    = NA_real_,
      n_combo_points     = n_combo
    ))
  }

  sa1_x <- sa1[order(sa1[[conc1_col]]), ][["x"]]
  sa2_x <- sa2[order(sa2[[conc2_col]]), ][["x"]]

  # HSA expected: the more potent single agent at each grid point
  expected <- as.vector(outer(sa1_x, sa2_x, pmin))
  combo_x <- combo[["x"]]
  excess <- expected - mean(combo_x, na.rm = TRUE)
  q90 <- stats::quantile(excess, 0.9, na.rm = TRUE)

  list(
    normalization_type = norm_type,
    hss_score          = mean(excess[excess >= q90], na.rm = TRUE),
    hss_excess_mean    = mean(excess, na.rm = TRUE),
    n_combo_points     = n_combo
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
  cell_pairs <- unique(fit_dt[, c("row", "column")])

  summary_rows <- lapply(seq_len(NROW(cell_pairs)), function(i) {
    r <- cell_pairs[["row"]][i]
    cc <- cell_pairs[["column"]][i]
    cell_fit <- fit_dt[fit_dt[["row"]] == r & fit_dt[["column"]] == cc, ]

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
    new_key <- do.call(paste, new_dt[, key_cols_present, with = FALSE])
    existing_key <- do.call(paste, existing_dt[, key_cols_present, with = FALSE])
    existing_pruned <- existing_dt[!existing_key %in% unique(new_key), ]

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
    normalization_type    = norm_type,
    x_mean                = NA_real_,
    x_AOC                 = NA_real_,
    N_conc                = 0L,
    maxlog10Concentration = NA_real_,
    xc50                  = NA_real_,
    h                     = NA_real_,
    r2                    = NA_real_,
    x_0                   = NA_real_,
    x_inf                 = NA_real_,
    fit_type              = "DRCInvalidFitResult"
  )
}


#' @keywords internal
.estimate_xc50_fallback <- function(x) {
  if (all(x > 0.5, na.rm = TRUE)) {
    Inf
  } else if (all(x <= 0.5, na.rm = TRUE)) {
    -Inf
  } else {
    NA
  }
}


####
# Deprecated alias — kept only until .persist_metrics callers are updated
####

#' @keywords internal
.persist_metrics <- function(se, new_metrics, merge, metrics_assay, row, col) {
  .persist_assay(
    se, new_metrics, merge, metrics_assay, row, col,
    upsert_key_cols = c("fit_source", "normalization_type")
  )
}
