#' Testing synthetic data form gDRtestData package
#'
#' @param original original MAE assay
#' @param data datase MAE or data.table
#' @param dataName dataset name
#' @param override_untrt_controls named list containing defining factors in 
#' the treatments
#' @param assays assays to test
#' @param tolerance tolerance factor
#' 
#' @examples
#' set.seed(2)
#' cell_lines <- gDRtestData::create_synthetic_cell_lines()
#' drugs <- gDRtestData::create_synthetic_drugs()
#' data <- "finalMAE_small"
#' original <- gDRutils::get_synthetic_data(data)
#' test_synthetic_data(original, original, "test")
#' 
#' @return `NULL`
#' @keywords test_utils
#' @export
test_synthetic_data <- function(original,
                                data,
                                dataName,
                                override_untrt_controls = NULL,
                                assays = c("Normalized", "Averaged", "Metrics"),
                                tolerance = 10e-4) {
  if (inherits(data, "MultiAssayExperiment")) {
    reprocessed <- data
  } else {
    reprocessed <- runDrugResponseProcessingPipeline(
      data,
      override_untrt_controls = override_untrt_controls
    )
  }
  
  masked_tag <- gDRutils::get_env_identifiers("masked_tag")
  
  purrr::walk(assays, function(x) {
              dt_original <- gDRutils::convert_mae_assay_to_dt(original, x)
              if (masked_tag %in% names(dt_original)) {
                dt_original[, (masked_tag) := NULL]
              }
              
              dt_reprocessed <- gDRutils::convert_mae_assay_to_dt(data, x)
              if (masked_tag %in% names(dt_reprocessed)) {
                dt_reprocessed[, (masked_tag) := NULL]
              }

              data.table::setcolorder(dt_original, names(dt_reprocessed))
              data.table::setorderv(dt_original, names(dt_original))
              data.table::setorderv(dt_reprocessed, names(dt_reprocessed))
              
              testthat::test_that(
                sprintf(
                  "reprocessed data %s is identical to data stored in gDRtestData", 
                  dataName
                ), {
                  testthat::expect_equal(
                    dt_original, dt_reprocessed, tolerance = tolerance
                  )
                })
  })
}
