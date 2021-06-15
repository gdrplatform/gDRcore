original <- readRDS(system.file("testdata", "finalSE_combo_codilution_small.RDS", package = "gDRtestData"))

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[1:2,], Drugs[1:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[1,,drop = FALSE], df_layout[,'Concentration',drop = FALSE])
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- cbind(df_layout, df_2)
df_layout_2 = df_layout_2[df_layout_2$DrugName != df_layout_2$DrugName_2,]
df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] <- 
  df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] / 2

df_merged_data <- generate_response_data(df_layout_2, 0)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE)
