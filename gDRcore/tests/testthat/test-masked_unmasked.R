test_that("masked and unmasked values are processed properly", {
  data <- "finalSE_small.RDS"
  original <- get_synthetic_data(data)
  
  df_layout <- merge(cell_lines[2:11, ], drugs[2:11, ], by = NULL)
  df_layout <- gDRtestData::add_data_replicates(df_layout)
  df_layout <- gDRtestData::add_concentration(df_layout)
  
  df_merged_data <- gDRtestData::generate_response_data(df_layout)
  
  df_merged_data$masked <- FALSE
  df_merged_data[df_merged_data$clid == df_merged_data$clid[[3]] &
                   df_merged_data$Gnumber ==  unique(df_merged_data$Gnumber)[[3]],
                 "masked"] <- TRUE
  
  
  finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)
  testthat::expect_s4_class(finalSE, "SummarizedExperiment")
})

