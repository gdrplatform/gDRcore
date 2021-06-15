original <- readRDS(system.file("testdata", "finalSE_combo_codilution.RDS", package = "gDRtestData"))

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[seq(1,15,2),], Drugs[1:12,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[c(1,1),], df_layout[,'Concentration',drop=F])
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- cbind(df_layout, df_2)
df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] <- 
  df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] / 2

df_merged_data <- generate_response_data(df_layout_2)
finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE)