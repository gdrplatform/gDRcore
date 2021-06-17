#' Run drug response processing pipeline
#'
#' Create a SummarizedExperiment of raw data and proceed to 
#' normalize, average, and fit the processed data. 
#'
#' @param df_ data.frame of raw drug response data containing both treated and untreated values. 
#' @param readout string of the name containing the cell viability readout values.
#' @param control_mean_fxn function indicating how to average controls.
#' Defaults to \code{mean(x, trim = 0.25)}.
#' @param nested_keys character vector of column names to include in the data.frames in the assays of the resulting \code{SummarizedExperiment} object.
#' Defaults to \code{c("Barcode", gDRutils::get_identifier("masked_tag"))}.
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
#' @param ref_GR_assay string of the name of the reference GR assay in the \linkS4class{SummarizedExperiment}.
#' Defaults to \code{"RefGRvalue"}.
#' @param ref_RV_assay string of the name of the reference Relative Viability assay in the \linkS4class{SummarizedExperiment}.
#' Defaults to \code{"RefRelativeViability"}.
#' @param metrics_assay string of the name of the metrics assay to output in the returned \linkS4class{SummarizedExperiment}
#' Defaults to \code{"Metrics"}.
#'
#'
#' @family runDrugResponseProcessingPipelineFxns
#' @export
#'
runDrugResponseProcessingPipeline <- function(df_, 
                                               readout = "ReadoutValue",
                                               control_mean_fxn = function(x) {mean(x, trim = 0.25)},
                                               nested_keys = c("Barcode", gDRutils::get_identifier("masked_tag")),
                                               override_untrt_controls = NULL,
                                               override_masked = FALSE,
                                               ndigit_rounding = 4,
                                               n_point_cutoff = 4,
                                               control_assay = "Controls",
                                               raw_treated_assay = "RawTreated",
                                               normalized_assay = "Normalized",
                                               ref_GR_assay = "RefGRvalue",
                                               ref_RV_assay = "RefRelativeViability",
                                               averaged_assay = "Averaged",
                                               metrics_assay = "Metrics") {
  se <- create_SE(df_ = df_, 
                   readout = readout, 
                   control_mean_fxn = control_mean_fxn, 
                   nested_keys = nested_keys, 
                   override_untrt_controls = override_untrt_controls)
  se <- normalize_SE(se = se, 
                      control_assay = control_assay, 
                      raw_treated_assay = raw_treated_assay, 
                      normalized_assay = normalized_assay,
                      ref_GR_assay = ref_GR_assay, 
                      ref_RV_assay = ref_RV_assay, 
                      ndigit_rounding = ndigit_rounding)
  se <- average_SE(se = se, 
                    override_masked = override_masked, 
                    normalized_assay = normalized_assay, 
                    averaged_assay = averaged_assay)
  se <- fit_SE(se = se, 
                averaged_assay = averaged_assay, 
                ref_GR_assay = ref_GR_assay, 
                metrics_assay = metrics_assay, 
                n_point_cutoff = n_point_cutoff)
  se <- add_codrug_group_SE(se)
  se
}
