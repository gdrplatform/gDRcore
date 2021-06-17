original <- readRDS(system.file("testdata", "finalSE_combo_matrix_small.RDS", package = "gDRtestData"))
normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")

df_layout <- merge(gDRtestData::create_synthetic_cell_lines()[7:8,], gDRtestData::create_synthetic_drugs()[c(4:6),], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout, Concentrations = 10**(seq(-3,.5,.5)))

df_2 <- merge(gDRtestData::create_synthetic_cell_lines()[gDRtestData::create_synthetic_cell_lines()$clid %in% df_layout$clid,], gDRtestData::create_synthetic_drugs()[c(21,26),], by = NULL)
df_2 <- gDRtestData::add_data_replicates(df_2)
df_2 <- gDRtestData::add_concentration(df_2, Concentrations = 10**(seq(-3,.5,.5)))
colnames(df_2)[colnames(df_2) %in% c(colnames(gDRtestData::create_synthetic_drugs()),'Concentration')] <- 
  paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(gDRtestData::create_synthetic_drugs()),'Concentration')], '_2')

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, "finalSE_combo_matrix_small.RDS")

