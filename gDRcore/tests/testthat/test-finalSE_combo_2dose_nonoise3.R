original <- readRDS(system.file("testdata", "finalSE_combo_2dose_nonoise3.RDS", package = "gDRtestData"))

df_layout <- merge(gDRtestData::create_synthetic_cell_lines()[2:4,], gDRtestData::create_synthetic_drugs()[2:4,], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_2 <- cbind(gDRtestData::create_synthetic_drugs()[c(26,26,26), ], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)
df_layout_2 = df_layout_2[!(df_layout_2$Concentration == 0 & df_layout_2$Concentration_2 > 0), ]

df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, "finalSE_combo_2dose_nonoise3.RDS")
