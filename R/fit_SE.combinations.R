#' fit_SE for combination screens
#'
#' Perform fittings for combination screens.
#' @param se \code{SummarizedExperiment} object with a BumpyMatrix assay
#' containing averaged data.
#' @param data_type single-agent vs combination
#' @param series_identifiers character vector of the column names in the
#' nested \code{DFrame} corresponding to nested identifiers.
#' @param normalization_types character vector of normalization types used for
#' calculating combo matrix.
#' @param averaged_assay string of the name of the averaged assay to use as
#' input. in the \code{se}.
#' @param metrics_assay string of the name of the metrics assay to output
#' in the returned \linkS4class{SummarizedExperiment}.
#' whose combination represents a unique series for which to fit curves.
#' @param score_FUN function used to calculate score for HSA and Bliss
#'
#' @details
#' This function assumes that the combination is set up with both
#' concentrations nested in the assay.
#'
#' @examples
#' fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#'
#' se1 <- fmae_cms[[gDRutils::get_supported_experiments("combo")]]
#' SummarizedExperiment::assays(se1) <-
#'   SummarizedExperiment::assays(se1)["Averaged"]
#' fit_SE.combinations(se1[1, 1])
#'
#' @keywords runDrugResponseProcessingPipeline
#' @return A \code{SummarizedExperiment} object with an additional assay
#' containing the combination metrics.
#'
#' @export
#'
fit_SE.combinations <- function(se,
                                data_type = gDRutils::get_supported_experiments("combo"),
                                series_identifiers = NULL,
                                normalization_types = c("GR", "RV"),
                                averaged_assay = "Averaged",
                                metrics_assay = "Metrics",
                                score_FUN = calculate_score) {

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE)
  checkmate::assert_character(normalization_types)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_function(score_FUN)

  if (is.null(series_identifiers)) {
    series_identifiers <- get_default_nested_identifiers(
      se,
      data_model(data_type)
    )
  }

  if (length(series_identifiers) != 2L) {
    stop(
      "gDR only supports 'series_identifiers' arguments with length '2'
      for the 'fit_SE.combinations' function"
    )
  }

  # Delegate to composable apply_combo_*() steps:
  #   Step 1+2: SA fits → Metrics
  #   Step 3:   smooth → excess
  #   Step 4:   Loewe CI → isobolograms + all_iso_points
  #   Step 5:   scores (read from excess for numerical identity)
  excess_assay <- "excess"
  se <- apply_combo_sa_fits(
    se,
    series_identifiers = series_identifiers,
    normalization_types = normalization_types,
    averaged_assay = averaged_assay,
    metrics_assay = metrics_assay,
    fit_source = "gDR"
  )
  se <- apply_combo_excess(
    se,
    series_identifiers = series_identifiers,
    normalization_types = normalization_types,
    averaged_assay = averaged_assay,
    metrics_assay = metrics_assay,
    excess_assay = excess_assay
  )
  se <- apply_combo_isobolograms(
    se,
    series_identifiers = series_identifiers,
    normalization_types = normalization_types,
    averaged_assay = averaged_assay,
    metrics_assay = metrics_assay,
    isobolograms_assay = "isobolograms",
    iso_points_assay = "all_iso_points"
  )
  se <- apply_combo_scores(
    se,
    scores_assay = "scores",
    averaged_assay = averaged_assay,
    metrics_assay = metrics_assay,
    excess_assay = excess_assay,
    normalization_types = normalization_types,
    fit_source = "gDR",
    score_FUN = score_FUN
  )
  se
}
