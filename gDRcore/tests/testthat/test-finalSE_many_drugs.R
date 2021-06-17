original <- readRDS(system.file("testdata", "finalSE_many_drugs.RDS", package = "gDRtestData"))

df_layout <- merge(gDRtestData::create_synthetic_cell_lines()[1:10,], gDRtestData::create_synthetic_drugs(), by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_merged_data <- gDRtestData::generate_response_data(df_layout)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, "finalSE_many_drugs.RDS")
