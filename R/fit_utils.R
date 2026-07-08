#' apply_fit_to_se
#'
#' Apply a user-supplied fit function once per (drug \eqn{\times} cell line
#' \eqn{\times} normalization_type) triplet of the Averaged
#' \linkS4class{BumpyMatrix} assay and persist the results into the Metrics
#' assay with idempotent merge semantics and configurable error handling.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   response data.  Must contain the assay named by \code{averaged_assay}.
#' @param fit_fn A function accepting a single \code{\link[data.table]{data.table}}
#'   (the Averaged assay slice pre-filtered to one drug \eqn{\times} cell line
#'   \eqn{\times} normalization_type triplet) and returning a named list of
#'   scalar metric values.
#' @param normalization_types Character vector of normalization types to iterate.
#'   Defaults to \code{c("GR", "RV")}.
#' @param averaged_assay String name of the averaged assay in \code{se}.
#'   Defaults to \code{"Averaged"}.
#' @param metrics_assay String name of the metrics assay in \code{se}.
#'   Defaults to \code{"Metrics"}.
#' @param merge String controlling how new results are combined with any
#'   existing Metrics assay data.  \code{"merge"} (default) performs an
#'   idempotent upsert keyed by (\code{fit_source} +
#'   \code{normalization_type}): rows whose key already exists are replaced,
#'   all other existing rows are preserved.  \code{"replace"} overwrites the
#'   entire Metrics assay.
#' @param on_error String controlling behaviour when \code{fit_fn} throws an
#'   error for a given triplet.  \code{"warn"} (default) emits a warning and
#'   skips that entry.  \code{"stop"} propagates the error immediately without
#'   any partial write.
#' @param fit_source String identifier recorded as the \code{fit_source}
#'   column in the Metrics assay output.  Used as an upsert key together with
#'   \code{normalization_type} when \code{merge = "merge"}.
#'   Defaults to \code{"custom"}.
#'
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} with the
#'   Metrics assay updated according to the \code{merge} strategy.
#'
#' @details
#' For each (drug \eqn{\times} cell line \eqn{\times} normalization_type)
#' triplet the function:
#' \enumerate{
#'   \item Extracts the Averaged assay slice for the given drug/cell line pair
#'         (via \code{\link[gDRutils]{apply_bumpy_function}}).
#'   \item Filters rows where the \code{normalization_type} column matches the
#'         current normalization type (all rows are passed through unchanged
#'         when the \code{normalization_type} column is absent).
#'   \item Calls \code{fit_fn} on the resulting \code{data.table}.
#'   \item Records \code{fit_source} and \code{normalization_type} metadata on
#'         the returned list before converting it to a data.table row.
#' }
#' If \code{fit_fn} output is missing any column listed in
#' \code{gDRutils::get_header("response_metrics")} a warning (not an error)
#' is emitted once after all triplets are processed.
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

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_function(fit_fn)
  checkmate::assert_character(normalization_types, min.len = 1L)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_choice(merge, c("merge", "replace"))
  checkmate::assert_choice(on_error, c("warn", "stop"))
  checkmate::assert_string(fit_source)

  if (!averaged_assay %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf(
      "%s assay is required for %s()",
      averaged_assay,
      as.character(match.call()[[1]])
    ))
  }

  norm_col <- "normalization_type"

  wrapper_fn <- function(avg_df) {
    avg_dt <- data.table::as.data.table(avg_df)
    results <- list()

    for (nt in normalization_types) {
      if (norm_col %in% names(avg_dt)) {
        sub_dt <- avg_dt[avg_dt[[norm_col]] == nt, ]
      } else {
        sub_dt <- avg_dt
      }

      if (NROW(sub_dt) == 0L) {
        next
      }

      result_dt <- tryCatch({
        res <- fit_fn(sub_dt)
        data.table::as.data.table(as.list(res))
      }, error = function(e) {
        if (on_error == "stop") {
          stop(e)
        }
        warning(sprintf(
          "fit_fn failed for normalization_type=%s: %s",
          nt, conditionMessage(e)
        ))
        NULL
      })

      if (is.null(result_dt)) {
        next
      }

      result_dt$fit_source <- fit_source
      result_dt$normalization_type <- nt
      results[[length(results) + 1L]] <- result_dt
    }

    if (length(results) == 0L) {
      return(data.table::data.table())
    }
    data.table::rbindlist(results, fill = TRUE)
  }

  se_with_new <- gDRutils::apply_bumpy_function(
    se = se,
    FUN = wrapper_fn,
    req_assay_name = averaged_assay,
    out_assay_name = "__custom_fit_tmp__"
  )

  tmp_assay_name <- "__custom_fit_tmp__"
  if (!tmp_assay_name %in% SummarizedExperiment::assayNames(se_with_new)) {
    return(se)
  }

  new_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_with_new, tmp_assay_name),
    row.field = "row",
    column.field = "column"
  )
  new_metrics <- data.table::as.data.table(new_df)

  if (NROW(new_metrics) == 0L) {
    return(se)
  }

  expected_metrics <- gDRutils::get_header("response_metrics")
  missing_cols <- setdiff(expected_metrics, names(new_metrics))
  if (length(missing_cols) > 0L) {
    warning(sprintf(
      "fit_fn output is missing response_metrics columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  .persist_metrics(se, new_metrics, merge, metrics_assay, "row", "column")
}


#' fit_drug_response_metrics
#'
#' Reference fit function equivalent to current gDR RV/GR fitting.
#' Replicates what \code{\link[gDRutils]{logisticFit}} + \code{fit_curves()}
#' produce for one entry.
#'
#' @param avg_dt \code{data.table} of averaged data for one triplet.
#' @param capping_fold numeric capping fold. Defaults to \code{5}.
#'
#' @return Named list of metrics.
#'
#' @export
fit_drug_response_metrics <- function(avg_dt, capping_fold = 5) {
  norm_type <- avg_dt$normalization_type[1]
  conc_col  <- gDRutils::get_env_identifiers("concentration")
  conc      <- avg_dt[[conc_col]]
  x         <- avg_dt$x

  keep <- !is.na(x) & !is.na(conc)
  x <- x[keep]
  conc <- conc[keep]

  if (length(x) == 0L) {
    return(list(
      fit_source            = "custom",
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
    ))
  }

  x_mean                <- mean(x, na.rm = TRUE)
  x_AOC                 <- 1 - x_mean
  N_conc                <- length(unique(conc))
  maxlog10Concentration <- log10(max(conc, na.rm = TRUE))

  fit <- tryCatch(
    drc::drm(x ~ conc, data = data.frame(x = x, conc = conc), # nolint
             fct = drc::LL.4(),
             start   = c(2, 0.4, 1, stats::median(conc)),
             control = drc::drmc(relTol = 1e-06, errorm = FALSE,
                                 noMessage = TRUE, rmNA = TRUE)),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    p        <- stats::coef(fit)
    h <- p[1]
    x_inf <- p[2]
    x_0 <- p[3]
    ec50 <- p[4]
    xc50     <- gDRutils::cap_xc50(ec50, max(conc), capping_fold = capping_fold)
    r2       <- 1 - sum(stats::residuals(fit)^2) / sum((x - mean(x))^2)
    fit_type <- "DRC4pHillFitModel"
  } else {
    h <- x_inf <- x_0 <- r2 <- NA_real_
    if (all(x > 0.5, na.rm = TRUE)) {
      xc50 <- Inf
    } else if (all(x <= 0.5, na.rm = TRUE)) {
      xc50 <- -Inf
    } else {
      xc50 <- NA
    }
    fit_type <- "DRCInvalidFitResult"
  }

  list(
    fit_source            = "custom",
    normalization_type    = norm_type,
    x_mean                = x_mean,
    x_AOC                 = x_AOC,
    N_conc                = N_conc,
    maxlog10Concentration = maxlog10Concentration,
    xc50                  = xc50,
    h                     = h,
    r2                    = r2,
    x_0                   = x_0,
    x_inf                 = x_inf,
    fit_type              = fit_type
  )
}


#' Persist new metrics into the Metrics assay of a SummarizedExperiment
#'
#' @param se \code{SummarizedExperiment}.
#' @param new_metrics \code{data.table} of new metrics.
#' @param merge One of \code{"merge"} or \code{"replace"}.
#' @param metrics_assay String naming the metrics assay.
#' @param row Name of the row-index column.
#' @param col Name of the column-index column.
#'
#' @return Updated \code{SummarizedExperiment}.
#'
#' @keywords internal
.persist_metrics <- function(se, new_metrics, merge, metrics_assay, row, col) {
  has_metrics <- metrics_assay %in% SummarizedExperiment::assayNames(se)

  if (merge == "merge" && has_metrics) {
    existing_df <- BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se, metrics_assay),
      row.field = row,
      column.field = col
    )
    existing_dt <- data.table::as.data.table(existing_df)

    if (!"fit_source" %in% names(existing_dt)) {
      existing_dt[, fit_source := "gDR"]
    }
    if (!"normalization_type" %in% names(existing_dt)) {
      existing_dt[, normalization_type := NA_character_]
    }

    new_key <- paste(
      new_metrics[["fit_source"]],
      new_metrics[["normalization_type"]]
    )
    existing_key <- paste(
      existing_dt[["fit_source"]],
      existing_dt[["normalization_type"]]
    )
    existing_pruned <- existing_dt[!existing_key %in% unique(new_key), ]

    metrics_merged <- data.table::rbindlist(
      list(existing_pruned, new_metrics),
      fill = TRUE
    )
  } else {
    metrics_merged <- new_metrics
  }

  data_cols <- names(metrics_merged)[!names(metrics_merged) %in% c(row, col)]
  metrics_mx <- BumpyMatrix::splitAsBumpyMatrix(
    metrics_merged[, data_cols, with = FALSE],
    row = metrics_merged[[row]],
    col = metrics_merged[[col]]
  )
  SummarizedExperiment::assay(se, metrics_assay) <- metrics_mx
  se
}
