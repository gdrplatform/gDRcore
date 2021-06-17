original <- readRDS(system.file("testdata", "finalSE_combo_matrix.RDS", package = "gDRtestData"))

df_layout <- merge(gDRtestData::create_synthetic_cell_lines()[seq(1,30,4),], gDRtestData::create_synthetic_drugs()[c(1,2,11,12,16,17),], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_2 <- merge(gDRtestData::create_synthetic_cell_lines()[gDRtestData::create_synthetic_cell_lines()$clid %in% df_layout$clid,], gDRtestData::create_synthetic_drugs()[c(21,26,31),], by = NULL)
df_2 <- gDRtestData::add_data_replicates(df_2)
df_2 <- gDRtestData::add_concentration(df_2)
colnames(df_2)[colnames(df_2) %in% c(colnames(gDRtestData::create_synthetic_drugs()),'Concentration')] <- 
  paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(gDRtestData::create_synthetic_drugs()),'Concentration')], '_2')

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- gDRtestData::generate_response_data(df_layout_2)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, "finalSE_combo_matrix.RDS")

