data <- "finalMAE_many_drugs.RDS"
original <- get_synthetic_data(data)

df_layout <- merge(cell_lines[1:10, ], drugs, by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_merged_data <- gDRtestData::generate_response_data(df_layout)

finalMAE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalMAE, data)
