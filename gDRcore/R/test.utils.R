#' @export
#'
test_synthetic_data <- function(original,
                                data,
                                dataName,
                                additional_columns = 0,
                                override_untrt_controls = NULL,
                                tolerance = 10e-7,
                                combo = FALSE,
                                OMITTED_COLUMNS_TO_TEST_NORMALIZED = c("rId", "cId"),
                                OMITTED_COLUMNS_TO_TEST_AVERAGED = c("rId", "cId"),
                                OMITTED_COLUMNS_TO_TEST_METRICS = c("rId", "cId")
                                ) {
  if (inherits(data, "MultiAssayExperiment")) {
    reprocessed <- data
  } else {
    reprocessed <- gDRcore::runDrugResponseProcessingPipeline(data,
                                                              override_untrt_controls = override_untrt_controls)
  }
  
  if (!is.null(override_untrt_controls)) {
    original <- original[SummarizedExperiment::rowData(original)[[names(override_untrt_controls)]]
                         == override_untrt_controls, ]
    reprocessed[[1]] <- reprocessed[[1]][SummarizedExperiment::rowData(reprocessed[[1]])
                                         [[names(override_untrt_controls)]]
                               == override_untrt_controls, ]
  }
  reprocessed <- reprocessed[["single-agent"]]
  if (gDRutils::get_env_identifiers("drug_name2") %in% names(SummarizedExperiment::rowData(original))) {
    original <- original[SummarizedExperiment::rowData(original)[[gDRutils::get_env_identifiers("drug_name2")]]
                         %in% gDRutils::get_env_identifiers("untreated_tag")]
  }
  
  normalized <- as.data.frame(gDRutils::convert_se_assay_to_dt(original, "Normalized"))
  averaged <- as.data.frame(gDRutils::convert_se_assay_to_dt(original, "Averaged"))
  normalized_new <- as.data.frame(gDRutils::convert_se_assay_to_dt(reprocessed, "Normalized"))
  averaged_new <- as.data.frame(gDRutils::convert_se_assay_to_dt(reprocessed, "Averaged"))

  if (!combo) {
    metrics <- as.data.frame(gDRutils::convert_se_assay_to_dt(original, "Metrics"))
    metrics_new <- as.data.frame(gDRutils::convert_se_assay_to_dt(reprocessed, "Metrics"))
  }
  
  if (combo) {
    for (assay in c("averaged")) {
      data.table::setDT(normalized)
      data.table::setDT(averaged)
      data.table::setDT(normalized_new)
      data.table::setDT(averaged_new)
      refColNames <- intersect(unname(unlist(gDRutils::get_env_identifiers())), names(get(assay)))
      concCols <- grep("Concentration", refColNames, value = TRUE)
      original <- unique(get(assay)[get(assay)[, apply(.SD != 0, 1, all), .SDcols = concCols], ])
      new <- unique(get(paste0(assay, "_new"))[
        get(paste0(assay, "_new"))[, apply(.SD != 0, 1, all), .SDcols = concCols], ])

      original$Concentration_2 <- round(original$Concentration_2, 7)
      new$Concentration_2 <- round(new$Concentration_2, 7)
      if (assay == "averaged") {
        original <- data.frame(lapply(original, function(x) if (is.numeric(x)) round(x, 4) else x))
        new <- data.frame(lapply(new, function(x) if (is.numeric(x)) round(x, 4) else x))
        data.table::setDT(original)
        data.table::setDT(new)
      }
      colsCompare <- setdiff(colnames(new), c(refColNames, "rId", "cId"))
      
      data.table::setorderv(new, colsCompare)
      data.table::setorderv(original, colsCompare)
      test_that(sprintf("Original data %s and recreated data are identical", dataName), {
      expect_equal(new[, ..colsCompare], original[, ..colsCompare])
      })
    }
  } else if (all(c("cotreatment", "single-agent") %in% names(reprocessed))) {
    
    normalized_new$Concentration_2[is.na(normalized_new$Concentration_2)] <- 0
    averaged_new$Concentration_2[is.na(averaged_new$Concentration_2)] <- 0
    normalized_new$Concentration_2[is.na(normalized_new$Concentration_2)] <- 0
    metrics_new$Concentration_2[is.na(metrics_new$Concentration_2)] <- 0
    order_cols <- unlist(gDRutils::get_env_identifiers(c("drug", "cellline", "concentration"), simplify = FALSE))
    
    
    test_that(sprintf("Original data %s and recreated data are identical", dataName), {
      expect_equal(ncol(normalized), ncol(normalized_new) + length(cotrt_cols))
      expect_equal(ncol(averaged), ncol(averaged_new) + length(cotrt_cols))
      expect_equal(ncol(metrics), ncol(metrics_new) + length(cotrt_cols))
      
      cotrt_cols <- grep("_2", names(normalized), value = TRUE)
      expect_equivalent(
        subset(normalized_new[do.call(order, normalized_new[order_cols]), ],
               select = which(!colnames(normalized_new) %in% c(OMITTED_COLUMNS_TO_TEST_NORMALIZED, cotrt_cols))),
        subset(normalized[do.call(order, normalized[order_cols]), ],
               select = which(!colnames(normalized) %in% c(OMITTED_COLUMNS_TO_TEST_NORMALIZED, cotrt_cols)))
      )
      expect_equivalent(
      subset(averaged_new[do.call(order, averaged_new[order_cols]), ],
             select = which(!colnames(averaged_new) %in% c(OMITTED_COLUMNS_TO_TEST_AVERAGED, cotrt_cols))),
      subset(averaged[do.call(order, averaged[order_cols]), ],
             select = which(!colnames(averaged) %in% c(OMITTED_COLUMNS_TO_TEST_AVERAGED, cotrt_cols)))
    )
    })
  } else {
    cotrt_cols_norm <- grep("_2", names(normalized), value = TRUE)
    cotrt_cols_avg <- grep("_2", names(averaged), value = TRUE)
    
    test_that(sprintf("Original data %s and recreated data are identical", dataName), {
      expect_equal(ncol(normalized), ncol(normalized_new) + length(cotrt_cols_norm))
      expect_equal(ncol(averaged), ncol(averaged_new) + length(cotrt_cols_avg))
      expect_equal(ncol(metrics), ncol(metrics_new) + length(cotrt_cols_avg))
      
      expect_equivalent(
        subset(normalized_new, select = which(!colnames(normalized_new) %in% OMITTED_COLUMNS_TO_TEST_NORMALIZED)),
        subset(normalized, select = which(!colnames(normalized) %in%
                                            c(OMITTED_COLUMNS_TO_TEST_NORMALIZED, cotrt_cols_norm)))
      )
      
      expect_equivalent(
        subset(averaged_new, select = which(!colnames(averaged_new) %in% OMITTED_COLUMNS_TO_TEST_AVERAGED)),
        subset(averaged, select = which(!colnames(averaged) %in% c(OMITTED_COLUMNS_TO_TEST_AVERAGED, cotrt_cols_avg)))
      )
      
      expect_equivalent(
        subset(metrics_new, select = which(!colnames(metrics_new) %in% OMITTED_COLUMNS_TO_TEST_METRICS)),
        subset(metrics, select = which(!colnames(metrics) %in% c(OMITTED_COLUMNS_TO_TEST_METRICS, cotrt_cols_avg))), 
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
