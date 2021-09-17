#' @export
#'
test_synthetic_data <- function(original,
                                data,
                                dataName,
                                additional_columns = 0,
                                override_untrt_controls = NULL,
                                tolerance = 10e-7,
                                combo = FALSE) {
  if (inherits(data, "SummarizedExperiment")) {
    reprocessed <- data
  } else {
    reprocessed <- gDRcore::runDrugResponseProcessingPipeline(data, override_untrt_controls = override_untrt_controls)
  }
  
  if (!is.null(override_untrt_controls)) {
    original <- original[SummarizedExperiment::rowData(original)[[names(override_untrt_controls)]]
                         == override_untrt_controls, ]
    reprocessed <- reprocessed[SummarizedExperiment::rowData(reprocessed)[[names(override_untrt_controls)]]
                               == override_untrt_controls, ]
  }
  normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
  averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
  metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")
  normalized_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Normalized")
  averaged_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Averaged")
  metrics_new <- gDRutils::convert_se_assay_to_dt(reprocessed, "Metrics")
  
  OMITTED_COLUMNS_TO_TEST_NORMALIZED <- c(
  )
  
  OMITTED_COLUMNS_TO_TEST_AVERAGED <- c(
  )
  
  OMITTED_COLUMNS_TO_TEST_METRICS <- c(
  )
  
  if (combo) {
    for (assay in c("normalized", "averaged", "metrics")) {
      refColNames <- intersect(unname(unlist(get_env_identifiers())), names(get(assay)))
      concCols <- grep("Concentration", refColNames, value = TRUE)
      original <- unique(get(assay)[get(assay)[, apply(.SD != 0, 1, all), .SDcols = refColNames], ])
      new <- unique(get(paste0(assay, "_new"))[get(paste0(assay, "_new"))[, apply(.SD != 0, 1, all),
                                                                          .SDcols = refColNames], ])
      original$Concentration_2 <- round(original$Concentration_2, 7)
      new$Concentration_2 <- round(new$Concentration_2, 7)
      
      if (assay == "averaged") {
        original <- data.frame(lapply(original, function(x) if(is.numeric(x)) round(x, 4) else x))
        new <- data.frame(lapply(new, function(x) if(is.numeric(x)) round(x, 4) else x))
      }
      
      d1 <- do.call("paste", original)
      d2 <- do.call("paste", new)
      expect_equal(dim(unique(original[d1%in%d2, ])), dim(original))
    }
  } else {
    test_that(sprintf("Original data %s and recreated data are identical", dataName), {
      expect_equal(ncol(normalized), 14 + additional_columns)
      expect_equal(ncol(averaged), 15 + additional_columns)
      expect_equal(ncol(metrics), 26 + additional_columns)
      
      expect_equivalent(
        subset(normalized_new, select = which(!colnames(normalized_new) %in% OMITTED_COLUMNS_TO_TEST_NORMALIZED)),
        subset(normalized, select = which(!colnames(normalized) %in% OMITTED_COLUMNS_TO_TEST_NORMALIZED))
      )
      
      expect_equivalent(
        subset(averaged_new, select = which(!colnames(averaged_new) %in% OMITTED_COLUMNS_TO_TEST_AVERAGED)),
        subset(averaged, select = which(!colnames(averaged) %in% OMITTED_COLUMNS_TO_TEST_AVERAGED))
      )
      
      expect_equivalent(
        subset(metrics_new, select = which(!colnames(metrics_new) %in% OMITTED_COLUMNS_TO_TEST_METRICS)),
        subset(metrics, select = which(!colnames(metrics) %in% OMITTED_COLUMNS_TO_TEST_METRICS)), 
        tolerance = tolerance
      )
    })
  }
}


# Test that the data is consistent after moving the RefReadout out of the create_and_normalize_SE logic.
#' @export
test_synthetic_data2 <- function(original, long_df, dataName, asys = c("Normalized", "Averaged")) {
  nested_ids <- intersect(c("Concentration", "Concentration_2"), colnames(long_df))
  reprocessed <- create_and_normalize_SE(long_df, nested_identifiers = nested_ids, nested_confounders = "Barcode")
  reprocessed <- average_SE(reprocessed, nested_identifiers = nested_ids)
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
