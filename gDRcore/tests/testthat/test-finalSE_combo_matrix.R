original <- readRDS(system.file("testdata", "finalSE_combo_matrix.RDS", package = "gDRtestData"))

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[seq(1,30,4),], Drugs[c(1,2,11,12,16,17),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- merge(CellLines[CellLines$clid %in% df_layout$clid,], Drugs[c(21,26,31),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2)
colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')] <- 
  paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')], '_2')

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- generate_response_data(df_layout_2)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE)

