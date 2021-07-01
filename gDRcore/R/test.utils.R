#' @export
#'
test_synthetic_data <- function(original, reprocessed, dataName) {
  normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
  averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
  metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")
  normalized_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Normalized")
  averaged_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Averaged")
  metrics_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Metrics")
  
  tolerance <- 10e-8
  test_that(sprintf("Original data %s and recreated data are identical", dataName), {
    expect_equal(normalized_new, normalized, tolerance = tolerance)
    expect_equal(averaged_new, averaged, tolerance = tolerance)
    expect_equal(metrics_new, metrics, tolerance = tolerance)
  })
}

#' @export
#'
get_synthetic_data <- function(rds) {
  readRDS(system.file("testdata", rds, package = "gDRtestData"))
}
