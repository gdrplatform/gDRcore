#' apply_fit_to_se
#'
#' Apply a user-supplied fit function per (drug x cell line x normalization_type)
#' triplet and persist results into the Metrics assay.
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
    stop(sprintf("%s assay is required", averaged_assay))
  }

  norm_col <- "normalization_type"
  tmp_assay <- "__custom_fit_tmp__"

  wrapper_fn <- .make_fit_wrapper(
    fit_fn, normalization_types, norm_col, on_error, fit_source
  )

  se_out <- gDRutils::apply_bumpy_function(
    se = se,
    FUN = wrapper_fn,
    req_assay_name = averaged_assay,
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

  # Warn once for missing standard metric columns
  expected <- gDRutils::get_header("response_metrics")
  missing <- setdiff(expected, names(new_metrics))
  if (length(missing) > 0L) {
    warning(sprintf("fit_fn output missing columns: %s", paste(missing, collapse = ", ")))
  }

  .persist_metrics(se, new_metrics, merge, metrics_assay, "row", "column")
}


####
# Internal helpers
####

#' Build wrapper function for apply_bumpy_function
#'
#' @param fit_fn user-supplied fit function
#' @param normalization_types character vector
#' @param norm_col column name for normalization type
#' @param on_error error handling strategy
#' @param fit_source label for fit_source column
#'
#' @return function suitable for apply_bumpy_function FUN
#'
#' @keywords internal
.make_fit_wrapper <- function(fit_fn, normalization_types, norm_col,
                              on_error, fit_source) {
  function(avg_df) {
    avg_dt <- data.table::as.data.table(avg_df)
    results <- list()

    for (nt in normalization_types) {
      sub_dt <- if (norm_col %in% names(avg_dt)) {
        avg_dt[avg_dt[[norm_col]] == nt, ]
      } else {
        avg_dt
      }

      if (NROW(sub_dt) == 0L) next

      result_dt <- tryCatch({
        res <- fit_fn(sub_dt)
        data.table::as.data.table(as.list(res))
      }, error = function(e) {
        if (on_error == "stop") stop(e)
        warning(sprintf("fit_fn failed for normalization_type=%s: %s",
                        nt, conditionMessage(e)))
        NULL
      })

      if (is.null(result_dt)) next

      result_dt$fit_source <- fit_source
      result_dt$normalization_type <- nt
      results[[length(results) + 1L]] <- result_dt
    }

    if (length(results) == 0L) return(data.table::data.table())
    data.table::rbindlist(results, fill = TRUE)
  }
}


####
# Reference fit function
####

#' fit_drug_response_metrics
#'
#' Reference fit function replicating standard gDR Hill curve fitting.
#'
#' @param avg_dt \code{data.table} of averaged data for one triplet
#' @param capping_fold numeric capping fold
#'
#' @return Named list of metrics
#'
#' @export
fit_drug_response_metrics <- function(avg_dt, capping_fold = 5) {
  norm_type <- avg_dt$normalization_type[1]
  conc_col <- gDRutils::get_env_identifiers("concentration")
  conc <- avg_dt[[conc_col]]
  x <- avg_dt$x

  # Filter NAs to ensure alignment for r2 calculation
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

  controls <- drc::drmc(relTol = 1e-06, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)
  fit <- tryCatch(
    drc::drm(x ~ conc,
             data = data.table::data.table(x = x, conc = conc),
             fct = drc::LL.4(),
             start = c(2, 0.4, 1, stats::median(conc)),
             control = controls),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    coefs <- stats::coef(fit)
    r2 <- 1 - sum(stats::residuals(fit)^2) / sum((x - mean(x))^2)
    out <- list(
      fit_source = "custom", normalization_type = norm_type,
      x_mean = x_mean, x_AOC = x_AOC, N_conc = N_conc,
      maxlog10Concentration = maxlog10Conc,
      xc50 = gDRutils::cap_xc50(coefs[4], max(conc), capping_fold = capping_fold),
      h = coefs[1], r2 = r2, x_0 = coefs[3], x_inf = coefs[2],
      fit_type = "DRC4pHillFitModel"
    )
  } else {
    xc50 <- .estimate_xc50_fallback(x)
    out <- list(
      fit_source = "custom", normalization_type = norm_type,
      x_mean = x_mean, x_AOC = x_AOC, N_conc = N_conc,
      maxlog10Concentration = maxlog10Conc,
      xc50 = xc50, h = NA_real_, r2 = NA_real_,
      x_0 = NA_real_, x_inf = NA_real_,
      fit_type = "DRCInvalidFitResult"
    )
  }

  out
}


#' @keywords internal
.empty_fit_result <- function(norm_type) {
  list(
    fit_source = "custom", normalization_type = norm_type,
    x_mean = NA_real_, x_AOC = NA_real_, N_conc = 0L,
    maxlog10Concentration = NA_real_, xc50 = NA_real_,
    h = NA_real_, r2 = NA_real_, x_0 = NA_real_, x_inf = NA_real_,
    fit_type = "DRCInvalidFitResult"
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
# Metrics persistence
####

#' Persist new metrics into the Metrics assay
#'
#' @param se SummarizedExperiment
#' @param new_metrics data.table of new metrics
#' @param merge \code{"merge"} or \code{"replace"}
#' @param metrics_assay assay name
#' @param row row index column name
#' @param col column index column name
#'
#' @return updated SummarizedExperiment
#'
#' @keywords internal
.persist_metrics <- function(se, new_metrics, merge, metrics_assay, row, col) {
  has_metrics <- metrics_assay %in% SummarizedExperiment::assayNames(se)

  if (merge == "merge" && has_metrics) {
    existing_df <- BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se, metrics_assay),
      row.field = row, column.field = col
    )
    existing_dt <- data.table::as.data.table(existing_df)

    # Ensure upsert keys exist in existing data
    if (!"fit_source" %in% names(existing_dt)) {
      existing_dt[, fit_source := "gDR"]
    }
    if (!"normalization_type" %in% names(existing_dt)) {
      existing_dt[, normalization_type := NA_character_]
    }

    new_key <- paste(new_metrics[["fit_source"]], new_metrics[["normalization_type"]])
    existing_key <- paste(existing_dt[["fit_source"]], existing_dt[["normalization_type"]])
    existing_pruned <- existing_dt[!existing_key %in% unique(new_key), ]

    metrics_merged <- data.table::rbindlist(list(existing_pruned, new_metrics), fill = TRUE)
  } else {
    metrics_merged <- new_metrics
  }

  data_cols <- names(metrics_merged)[!names(metrics_merged) %in% c(row, col)]
  metrics_mx <- BumpyMatrix::splitAsBumpyMatrix(
    metrics_merged[, data_cols, with = FALSE],
    row = metrics_merged[[row]], col = metrics_merged[[col]]
  )
  SummarizedExperiment::assay(se, metrics_assay) <- metrics_mx
  se
}
