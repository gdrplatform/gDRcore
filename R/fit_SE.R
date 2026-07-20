#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
fit_SE <- function(se,
                   data_type = "single-agent",
                   nested_identifiers = NULL,
                   averaged_assay = "Averaged",
                   metrics_assay = "Metrics",
                   n_point_cutoff = 4,
                   range_conc = c(5e-3, 5),
                   force_fit = FALSE,
                   pcutoff = 0.05,
                   cap = 0.1,
                   curve_type = c("GR", "RV")) {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(data_type)
  checkmate::assert_choice(data_type, choices = c(gDRutils::get_supported_experiments("sa"),
                                                  gDRutils::get_supported_experiments("cd")))
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_number(n_point_cutoff)
  checkmate::assert_numeric(range_conc)
  checkmate::assert_logical(force_fit)
  checkmate::assert_number(pcutoff)
  checkmate::assert_number(cap)
  checkmate::assert_character(curve_type)
  checkmate::assert_true(all(curve_type %in% c("GR", "RV")))
  gDRutils::validate_se_assay_name(se, averaged_assay)

  if (is.null(nested_identifiers)) {
    nested_identifiers <- get_default_nested_identifiers(
      se,
      data_model(data_type)
    )
  }

  # Build a fit_fn closure that forwards all fit_SE parameters.
  # slicing_values = curve_type filters to only the requested normalization types.
  fit_fn <- function(avg_dt) {
    fit_drug_response_metrics(
      avg_dt,
      capping_fold   = 5,
      range_conc     = range_conc,
      pcutoff        = pcutoff,
      n_point_cutoff = n_point_cutoff,
      force_fit      = force_fit,
      cap            = cap
    )
  }

  # slicing_cols: use nested_identifiers only when they extend beyond the default
  # "normalization_type" (e.g. co-dilution adds "ratio"). For standard SA data
  # the profile default ("normalization_type") is used.
  default_nested <- get_default_nested_identifiers(se, data_model(data_type))
  slicing_cols <- if (identical(nested_identifiers, default_nested)) {
    NULL  # let apply_fit use the profile default (normalization_type)
  } else {
    nested_identifiers
  }

  se <- apply_fit(
    se             = se,
    fit_fn         = fit_fn,
    data_type      = data_type,
    slicing_cols   = slicing_cols,
    slicing_values = curve_type,
    input_assay    = averaged_assay,
    output_assay   = metrics_assay,
    merge          = "replace",
    on_error       = "warn",
    fit_source     = "gDR"
  )

  se <- gDRutils::set_SE_fit_parameters(se,
    value = list(
      n_point_cutoff = n_point_cutoff,
      range_conc     = range_conc,
      force_fit      = force_fit,
      pcutoff        = pcutoff,
      cap            = cap)
  )
  se <- gDRutils::set_SE_processing_metadata(
    se,
    value = list(
      date_processed = Sys.Date(),
      session_info   = utils::sessionInfo()
    )
  )

  se
}


#' @keywords internal
fit_FUN <- function(x,
                    metric_cols = gDRutils::get_header("response_metrics"),
                    conc = gDRutils::get_env_identifiers("concentration"),
                    nested_identifiers,
                    n_point_cutoff,
                    range_conc,
                    force_fit,
                    pcutoff,
                    cap,
                    curve_type) {
  fit_df <- S4Vectors::DataFrame(matrix(NA, length(curve_type), length(metric_cols)))
  colnames(fit_df) <- metric_cols
  rownames(fit_df) <- c("RV", "GR")[c("RV", "GR") %in% curve_type]

  if (!is.null(x) && all(dim(x) > 0)) {
    if (!all(is.na(x[[conc]]))) {
      x <- x[x[[conc]] != 0, ]
    }

    fit_df <- S4Vectors::DataFrame(gDRutils::fit_curves(
      data.table::as.data.table(x),
      series_identifiers = nested_identifiers,
      e_0 = 1,
      GR_0 = 1,
      n_point_cutoff = n_point_cutoff,
      range_conc = range_conc,
      force_fit = force_fit,
      pcutoff = pcutoff,
      cap = cap,
      normalization_type = curve_type))
  }
  fit_df
}
