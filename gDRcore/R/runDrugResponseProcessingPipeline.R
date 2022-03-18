#' Run drug response processing pipeline
#'
#' Run different components of the gDR drug response processing pipeline.
#' Either: create a SummarizedExperiment and normalize raw treated and control data (create_and_normalize_SE),
#' average data (average_SE), or fit the processed data (fit_SE). See details for more in-depth explanations.
#'
#' @param df_ data.frame of raw drug response data containing both treated and untreated values. 
#' @param se \code{SummarizedExperiment} object.
#' @param readout string of the name containing the cell viability readout values.
#' @param control_mean_fxn function indicating how to average controls.
#' Defaults to \code{mean(x, trim = 0.25)}.
#' @param nested_keys character vector of column names to include in the data.frames
#' in the assays of the resulting \code{SummarizedExperiment} object.
#' @param Defaults to the \code{nested_identifiers} and \code{nested_confounders} if passed through
#' \code{create_and_normalize_SE} or \code{runDrugResponseProcessingPipeline}.
#' @param override_untrt_controls named list containing defining factors in the treatments.
#' Defaults to \code{NULL}.
#' @param override_masked boolean indicating whether or not to override the masked wells
#' in the averaging and include all wells. 
#' Defaults to \code{FALSE}.
#' @param ndigit_rounding integer indicating number of digits to round to in calculations.
#' Defaults to \code{4}.
#' @param n_point_cutoff integer of how many points should be considered the minimum required to try to fit a curve.
#' Defaults to \code{4}.
#' @param control_assay string containing the name of the assay representing the controls in the \code{se}.
#' Defaults to \code{"Controls"}.
#' @param raw_treated_assay string containing the name of the assay representing the raw treated data in the \code{se}.
#' Defaults to \code{"RawTreated"}.
#' @param normalized_assay string of the assay name containing the normalized data.
#' Defaults to \code{"Normalized"}.
#' @param averaged_assay string of the name of the averaged assay in the \linkS4class{SummarizedExperiment}.
#' Defaults to \code{"Averaged"}.
#' @param metrics_assay string of the name of the metrics assay to output
#' in the returned \linkS4class{SummarizedExperiment}
#' Defaults to \code{"Metrics"}.
#' @param add_raw_data  boolean indicating whether or not to include raw data into experiment metadata.
#'
#' @details
#' \code{runDrugResponseProcessingPipeline} is made up of 3 separate steps:
#' \itemize{
#'  \item{"create_and_normalize_SE"}
#'  \item{"average_SE"}
#'  \item{"fit_SE"}
#'}
#'
#' For create_and_normalize_SE, this creates a SummarizedExperiment object from a data.frame, 
#' where the data.frame contains treatments on rows, and conditions on columns. 
#' A \linkS4class{SummarizedExperiment} object containing two asssays is created:
#' treated readouts will live in an assay called \code{"RawTreated"},
#' and reference readouts live in an assay called \code{"Controls"}. Subsequently, the treated
#' and control elements will be normalized to output two metrics: 
#'
#' For average_SE, take the normalized assay and average the nested \code{DataFrame}s across unique
#' \code{nested_identifiers}. 
#'
#' For fit_SE, take the averaged assay and fit curves to obtain metrics, one set of metrics for each
#' normalization type set.
#'
#' @name runDrugResponseProcessingPipelineFxns
#'
NULL


#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
create_and_normalize_SE <- function(df_,
                                    readout = "ReadoutValue",
                                    control_mean_fxn = function(x) {
                                      mean(x, trim = 0.25)
                                    },
                                    nested_identifiers = NULL,
                                    nested_confounders = gDRutils::get_env_identifiers("barcode"),
                                    override_untrt_controls = NULL,
                                    ndigit_rounding = 4,
                                    control_assay = "Controls",
                                    raw_treated_assay = "RawTreated",
                                    normalized_assay = "Normalized") {
  
  checkmate::assert_data_frame(df_)
  checkmate::assert_string(readout)
  checkmate::assert_function(control_mean_fxn)
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  checkmate::assert_numeric(ndigit_rounding)
  checkmate::assert_string(control_assay)
  checkmate::assert_string(raw_treated_assay)
  checkmate::assert_string(normalized_assay)
  
  se <- create_SE(df_ = df_, 
    readout = readout, 
    control_mean_fxn = control_mean_fxn, 
    nested_identifiers = nested_identifiers,
    nested_confounders = nested_confounders,
    override_untrt_controls = override_untrt_controls)
  se <- normalize_SE(se = se, 
    nested_identifiers = nested_identifiers,
    nested_confounders = nested_confounders,
    control_assay = control_assay, 
    raw_treated_assay = raw_treated_assay, 
    normalized_assay = normalized_assay,
    ndigit_rounding = ndigit_rounding)
  se
}


#' @rdname runDrugResponseProcessingPipelineFxns 
#' @export
runDrugResponseProcessingPipeline <- function(df_,
                                              readout = "ReadoutValue",
                                              control_mean_fxn = function(x) {
                                                mean(x, trim = 0.25)
                                              },
                                              nested_identifiers = NULL,
                                              nested_confounders = gDRutils::get_env_identifiers("barcode"),
                                              override_untrt_controls = NULL,
                                              override_masked = FALSE,
                                              ndigit_rounding = 4,
                                              n_point_cutoff = 4,
                                              control_assay = "Controls",
                                              raw_treated_assay = "RawTreated",
                                              normalized_assay = "Normalized",
                                              averaged_assay = "Averaged",
                                              metrics_assay = "Metrics",
                                              add_raw_data = FALSE) {
  
  checkmate::assert_data_frame(df_)
  checkmate::assert_string(readout)
  checkmate::assert_function(control_mean_fxn)
  checkmate::assert_multi_class(nested_identifiers, c("character", "list"), null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  checkmate::assert_flag(override_masked)
  checkmate::assert_numeric(ndigit_rounding)
  checkmate::assert_number(n_point_cutoff)
  checkmate::assert_string(control_assay)
  checkmate::assert_string(raw_treated_assay)
  checkmate::assert_string(normalized_assay)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  
  df_ <- identify_data_type(df_)
  df_list <- split_raw_data(df_)
  mae <- MultiAssayExperiment::MultiAssayExperiment()
  if (!is.list(nested_identifiers)) {
    nested_identifiers <- if (is.null(nested_identifiers)) {
      list(`matrix` = .get_default_combo_identifiers(),
           `single-agent` = .get_default_single_agent_identifiers())
    } else if (length(nested_identifiers) == 1 || names(df_list) == 1) {
        list(`single-agent` = nested_identifiers)
    } else {
      stop("Number of detected data types is greater that 1.
           Please provide a named list of nested_identifiers")
    }
  }
  for (experiment in names(df_list)) {
    message("Processing ", experiment)
    experiment_identifier <- if (experiment %in% names(nested_identifiers)) {
      nested_identifiers[[experiment]]
    } else {
      nested_identifiers[["single-agent"]]
    }
    se <- purrr::quietly(create_and_normalize_SE)(df_ = df_list[[experiment]],
                                                     readout = readout,
                                                     control_mean_fxn = control_mean_fxn,
                                                     nested_identifiers = experiment_identifier,
                                                     nested_confounders = nested_confounders,
                                                     override_untrt_controls = override_untrt_controls,
                                                     control_assay = control_assay, 
                                                     raw_treated_assay = raw_treated_assay, 
                                                     normalized_assay = normalized_assay,
                                                     ndigit_rounding = ndigit_rounding)
    
    paste_warnings(se$warnings)
    se <- purrr::quietly(average_SE)(se = se$result, 
                                        series_identifiers = experiment_identifier,
                                        override_masked = override_masked, 
                                        normalized_assay = normalized_assay, 
                                        averaged_assay = averaged_assay)
    paste_warnings(se$warnings)
    se <- if (experiment == "matrix") {
      purrr::quietly(fit_SE.combinations)(se = se$result,
                          series_identifiers = experiment_identifier,
                          averaged_assay = averaged_assay)
    } else {
      purrr::quietly(fit_SE)(se = se$result, 
             nested_identifiers = experiment_identifier,
             averaged_assay = averaged_assay, 
             metrics_assay = metrics_assay, 
             n_point_cutoff = n_point_cutoff)
    }
    paste_warnings(se$warnings)
    if (add_raw_data) {
      se$result <- gDRutils::set_SE_experiment_raw_data(se$result, df_list[[experiment]])
    }
    mae <- c(mae, MultiAssayExperiment::MultiAssayExperiment(experiments = list(experiment = se$result)))
    names(mae)[length(names(mae))] <- experiment
  }
  mae
}

#' @keywords internal
paste_warnings <- function(list, sep = "\n") {
  warning(paste0(list, sep = sep), call. = FALSE)
}
