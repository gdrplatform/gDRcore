original <- readRDS(system.file("testdata", "finalSE_combo_matrix_small.RDS", package = "gDRtestData"))
normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[7:8,], Drugs[c(4:6),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout, Concentrations = 10**(seq(-3,.5,.5)))

df_2 <- merge(CellLines[CellLines$clid %in% df_layout$clid,], Drugs[c(21,26),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2, Concentrations = 10**(seq(-3,.5,.5)))
colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')] <- 
  paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')], '_2')

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- generate_response_data(df_layout_2, 0)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE)

