#' Testing synthetic data form gDRtestData package
#'
#' @param original original MAE assay
#' @param data datase MAE or data.table
#' @param dataName dataset name
#' @param override_untrt_controls named list containing defining factors in 
#' the treatments
#' @param tolerance tolerance factor
#' 
#' @examples
#' set.seed(2)
#' cell_lines <- gDRtestData::create_synthetic_cell_lines()
#' drugs <- gDRtestData::create_synthetic_drugs()
#' data <- "finalMAE_small.RDS"
#' original <- gDRutils::get_synthetic_data(data)
#' test_synthetic_data(original, original, "test")
#' 
#' @return `NULL`
#' @export
test_synthetic_data <- function(original,
                                data,
                                dataName,
                                override_untrt_controls = NULL,
                                tolerance = 10e-4) {
  if (inherits(data, "MultiAssayExperiment")) {
    reprocessed <- data
  } else {
    reprocessed <- runDrugResponseProcessingPipeline(
      data,
      override_untrt_controls = override_untrt_controls
    )
  }
  
  normalized <- 
    gDRutils::convert_mae_assay_to_dt(original, "Normalized")
  averaged <-
    gDRutils::convert_mae_assay_to_dt(original, "Averaged")
  metrics <-
    gDRutils::convert_mae_assay_to_dt(original, "Metrics")
  normalized_new <-
    gDRutils::convert_mae_assay_to_dt(reprocessed, "Normalized")
  averaged_new <-
    gDRutils::convert_mae_assay_to_dt(reprocessed, "Averaged")
  metrics_new <-
    gDRutils::convert_mae_assay_to_dt(reprocessed, "Metrics")

  
  testthat::test_that(
    sprintf(
      "reprocessed data %s is identical to data stored in gDRtestData", 
      dataName
    ), {
    testthat::expect_equal(
      normalized_new, normalized, tolerance = tolerance
    ) # nolint
    testthat::expect_equal(
      averaged_new, averaged, tolerance = tolerance
    ) # nolint
    testthat::expect_equal(
      metrics_new, metrics, tolerance = tolerance
    )
  })
}
