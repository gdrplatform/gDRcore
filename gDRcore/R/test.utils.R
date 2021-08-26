#' @export
#'
test_synthetic_data <- function(original, reprocessed, dataName) {
  normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
  averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
  metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")
  normalized_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Normalized")
  averaged_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Averaged")
  metrics_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Metrics")
  
  test_that(sprintf("Original data %s and recreated data are identical", dataName), {
    expect_equal(normalized_new, normalized)
    expect_equal(averaged_new, averaged)
    expect_equal(metrics_new, metrics, tolerance = 10e-4)
  })
}

#' @export
#'
test_synthetic_data2 <- function(original, reprocessed, dataName, asys = c("Normalized", "Averaged", "Metrics")) {
  for (asy in asys) {
    o_df <- S4Vectors::DataFrame(gDRutils::convert_se_assay_to_dt(original, asy))
    n_df <- S4Vectors::DataFrame(gDRutils::convert_se_assay_to_dt(reprocessed, asy))
    
    o_df$rId <- o_df$cId <- NULL
    n_df$rId <- n_df$cId <- NULL
   
    sort_cols <- c("Gnumber", "Gnumber_2", "clid", "Concentration", "Concentration_2")
    order_cols <- c(intersect(sort_cols, colnames(o_df)), setdiff(colnames(o_df), sort_cols))
    
    o_df <- unique(o_df)[order_cols]
    n_df <- unique(n_df)[order_cols]
    
    o_df <- o_df[S4Vectors::order(o_df), ]
    n_df <- n_df[S4Vectors::order(n_df), ]

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
