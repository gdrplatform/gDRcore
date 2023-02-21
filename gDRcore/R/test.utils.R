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
    testthat::expect_equal(normalized_new, normalized, tolerance = tolerance)
    testthat::expect_equal(averaged_new, averaged, tolerance = tolerance)
    testthat::expect_equal(metrics_new, metrics, tolerance = tolerance)
  })
}


#' Test that the data is consistent after moving the RefReadout out of the create_and_normalize_SE logic.
#'
#' @param original original MAE assay
#' @param long_df dataset in data frame format
#' @param dataName dataset name
#' @param asys names of assays
#'
#' @export
test_synthetic_data2 <- function(original, long_df, dataName, asys = c("Normalized", "Averaged")) {
  nested_ids <- intersect(c("Concentration", "Concentration_2"), colnames(long_df))
  reprocessed <- create_and_normalize_SE(long_df, nested_identifiers = nested_ids, nested_confounders = "Barcode")
  reprocessed <- average_SE(reprocessed, series_identifiers = nested_ids)
  for (asy in asys) {
    o_df <- S4Vectors::DataFrame(gDRutils::convert_se_assay_to_dt(original, asy))
    n_df <- S4Vectors::DataFrame(gDRutils::convert_se_assay_to_dt(reprocessed, asy))
    
    o_df$rId <- o_df$cId <- NULL
    n_df$rId <- n_df$cId <- NULL
   
    sort_cols <- c("Gnumber", "Gnumber_2", "clid", "Concentration", "Concentration_2")
    order_cols <- c(intersect(sort_cols, colnames(o_df)), setdiff(colnames(o_df), sort_cols))
    
    # Get rid of all untreated from old and references from new.    
    o_df <- o_df[o_df$Concentration == 0 & o_df$Concentration_2 == 0, ]
    n_df <- n_df[n_df$Concentration == 0 & n_df$Concentration_2 == 0, ]

    o_df <- unique(o_df)[order_cols]
    n_df <- unique(n_df)[order_cols]

    o_df <- o_df[S4Vectors::order(o_df), ]
    n_df <- n_df[S4Vectors::order(n_df), ]

    tolerance <- 10e-4

    testthat::test_that(sprintf("Original data %s and recreated data are identical for assay %s", dataName, asy), {
      testthat::expect_equal(o_df, n_df, tolerance = tolerance)
    })
  }
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
