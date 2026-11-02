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
  if (data_type == gDRutils::get_supported_experiments("time-course")) {
    stop("fit_SE() does not fit time-course data, use fit_SE.timecourse()")
  }
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

  # Pass the settings through fit_args rather than a closure, so the recipe that
  # apply_fit() records names gDRutils::fit_drug_response_metrics and its arguments.
  fit_args <- list(
    capping_fold = 5,
    range_conc = range_conc,
    pcutoff = pcutoff,
    n_point_cutoff = n_point_cutoff,
    force_fit = force_fit,
    cap = cap
  )

  # A co-dilution cell holds both concentrations; fit_codilution_metrics() fits each
  # ratio along the total concentration and returns one row per ratio.
  fit_fn <- if (data_type == gDRutils::get_supported_experiments("cd")) {
    gDRutils::fit_codilution_metrics
  } else {
    gDRutils::fit_drug_response_metrics
  }

  # The fit profile of data_type slices by normalization_type; slicing_values =
  # curve_type keeps only the requested normalization types.
  se <- gDRutils::apply_fit(
    se = se,
    fit_fn = fit_fn,
    data_type = data_type,
    slicing_values = curve_type,
    input_assay = averaged_assay,
    output_assay = metrics_assay,
    merge = "replace",
    on_error = "warn",
    fit_source = "gDR",
    fit_args = fit_args
  )

  # Ensure the Metrics assay always exists after fitting, even when apply_fit
  # returns the original SE without it (all-empty results scenario where
  # apply_bumpy_function errors on dimnames mismatch).
  if (!metrics_assay %in% SummarizedExperiment::assayNames(se)) {
    # Build an empty BumpyDataFrameMatrix with the full Metrics column schema
    # so that validate_SE() can verify normalization_type and fit_source are
    # present even when no fits succeeded.
    nr <- NROW(se)
    nc <- NCOL(se)
    metric_cols <- c(gDRutils::get_header("response_metrics"), "normalization_type", "fit_source")
    empty_df <- S4Vectors::DataFrame(
      matrix(character(0), nrow = 0, ncol = length(metric_cols)),
      check.names = FALSE
    )
    names(empty_df) <- metric_cols
    empty_bm <- BumpyMatrix::splitAsBumpyMatrix(
      empty_df,
      row = rep(seq_len(nr), nc),
      column = rep(seq_len(nc), each = nr)
    )
    dimnames(empty_bm) <- dimnames(se)
    SummarizedExperiment::assay(se, metrics_assay, withDimnames = FALSE) <- empty_bm
  }

  se <- gDRutils::set_SE_fit_parameters(se,
    value = list(
      n_point_cutoff = n_point_cutoff,
      range_conc = range_conc,
      force_fit = force_fit,
      pcutoff = pcutoff,
      cap = cap)
  )
  # the versions that produced the fit, as in the provenance run block; the whole
  # sessionInfo() also listed every package loaded at that moment, so two fits of
  # the same data gave different objects
  se <- gDRutils::set_SE_processing_metadata(
    se,
    value = list(
      date_processed = Sys.Date(),
      versions = list(
        R = paste(R.version$major, R.version$minor, sep = "."),
        gDRcore = as.character(utils::packageVersion("gDRcore")),
        gDRutils = as.character(utils::packageVersion("gDRutils"))
      )
    )
  )

  se
}
