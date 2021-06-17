original <- readRDS(system.file("testdata", "finalSE_combo_1dose_many_drugs.RDS", package = "gDRtestData"))

df_layout <- merge(gDRtestData::create_synthetic_cell_lines()[2:4,], gDRtestData::create_synthetic_drugs()[-1,], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_2 <- cbind(gDRtestData::create_synthetic_drugs()[c(1,1),], Concentration = c(0, 2))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- gDRtestData::generate_response_data(df_layout_2)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, "finalSE_combo_1dose_many_drugs.RDS")

