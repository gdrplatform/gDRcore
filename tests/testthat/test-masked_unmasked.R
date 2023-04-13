test_that("masked and unmasked values are processed properly", {
  data <- "finalMAE_small.RDS"
  original <- gDRutils::get_synthetic_data(data)
  
  df_layout <- merge(cell_lines[2:11, ], drugs[2:11, ], by = NULL)
  df_layout <- gDRtestData::add_data_replicates(df_layout)
  df_layout <- gDRtestData::add_concentration(df_layout)
  
  df_merged_data <- gDRtestData::generate_response_data(df_layout)
  
  df_merged_data$masked <- FALSE
  df_merged_data[df_merged_data$clid == df_merged_data$clid[[3]] &
                   df_merged_data$Gnumber ==  unique(df_merged_data$Gnumber)[[3]],
                 "masked"] <- TRUE
  df_merged_data <- data.table::setDT(df_merged_data)
  
  finalMAE <- purrr::quietly(runDrugResponseProcessingPipeline)(
    df_merged_data, 
    override_untrt_controls = NULL,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  expect_length(finalMAE$warnings, 3)
  
  testthat::expect_s4_class(finalMAE$result, "MultiAssayExperiment")
})

