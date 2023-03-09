#' Testing synthetic data form gDRtestData package
#'
#' @param original original MAE assay
#' @param data datase MAE or data frame
#' @param dataName dataset name
#' @param override_untrt_controls named list containing defining factors in the treatments
#' @param tolerance tolerance factor
#' 
#' @export
test_synthetic_data <- function(original,
                                data,
                                dataName,
                                override_untrt_controls = NULL,
                                tolerance = 10e-5) {
  if (inherits(data, "MultiAssayExperiment")) {
    reprocessed <- data
  } else {
    reprocessed <- gDRcore::runDrugResponseProcessingPipeline(data,
                                                              override_untrt_controls = override_untrt_controls)
  }
  
  normalized <- as.data.frame(gDRutils::convert_mae_assay_to_dt(original, "Normalized"))
  averaged <- as.data.frame(gDRutils::convert_mae_assay_to_dt(original, "Averaged"))
  metrics <- as.data.frame(gDRutils::convert_mae_assay_to_dt(original, "Metrics"))
  normalized_new <- as.data.frame(gDRutils::convert_mae_assay_to_dt(reprocessed, "Normalized"))
  averaged_new <- as.data.frame(gDRutils::convert_mae_assay_to_dt(reprocessed, "Averaged"))
  metrics_new <- as.data.frame(gDRutils::convert_mae_assay_to_dt(reprocessed, "Metrics"))
  
  
  testthat::test_that(sprintf("reprocessed data %s is identical to data stored in gDRtestData", dataName), {
    testthat::expect_equal(normalized_new, normalized, tolerance = tolerance) # nolint
    testthat::expect_equal(averaged_new, averaged, tolerance = tolerance) # nolint
    testthat::expect_equal(metrics_new, metrics, tolerance = tolerance)
  })
}

#' Get synthetic data from gDRtestData package
#'
#' @param rds RDS filename
#'
#' @return loaded data
#' 
#' @export
get_synthetic_data <- function(rds) {
  readRDS(system.file("testdata", rds, package = "gDRtestData"))
}
