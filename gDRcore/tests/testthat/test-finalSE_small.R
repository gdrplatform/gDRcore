library(testthat); library(gDRcore)
data <- "finalSE_small.RDS"
original <- get_synthetic_data(data)

df_layout <- merge(cell_lines[2:11, ], drugs[2:11, ], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_merged_data <- gDRtestData::generate_response_data(df_layout)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, data)
