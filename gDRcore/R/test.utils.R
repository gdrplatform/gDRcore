#' @export
#'
test_synthetic_data <- function(original, reprocessed, dataName) {
  normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
  averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
  metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")
  normalized_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Normalized")
  averaged_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Averaged")
  metrics_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Metrics")
  
  tolerance <- 10e-4
  test_that(sprintf("Original data %s and recreated data are identical", dataName), {
    expect_equal(normalized_new, normalized)
    expect_equal(averaged_new, averaged)
    expect_equal(metrics_new, metrics, tolerance = tolerance)
  })
}

#' @export
#'
test_synthetic_data2 <- function(original, reprocessed, dataName, asys = c("Normalized", "Averaged", "Metrics")) {
  # TODO: remove row_id, col_id
  # TODO: order rows appropriately
  for (asy in asys) {
    o_df <- S4Vectors::DataFrame(gDRutils::convert_se_assay_to_dt(original, asy))
    n_df <- S4Vectors::DataFrame(gDRutils::convert_se_assay_to_dt(reprocessed, asy))[names(o_dt)]
    
    tolerance <- 10e-4
    test_that(sprintf("Original data %s and recreated data are identical for assay %s", dataName, asy), {
      expect_equal(o_df, n_df, tolerance = tolerance)
    })
  }
}

#' @export
#'
get_synthetic_data <- function(rds) {
  readRDS(system.file("testdata", rds, package = "gDRtestData"))
}
