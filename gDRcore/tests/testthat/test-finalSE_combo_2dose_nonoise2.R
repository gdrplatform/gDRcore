original <- readRDS(system.file("testdata", "finalSE_combo_2dose_nonoise2.RDS", package = "gDRtestData"))

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[2:4,], Drugs[c(2:4,26),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[26 * c(1,1,1),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- generate_response_data(df_layout_2, 0)
df_merged_data <- df_merged_data[!(df_merged_data$Gnumber %in% c('vehicle', Drugs$Gnumber[26]) & 
                                     df_merged_data$Gnumber_2 == Drugs$Gnumber[26]),]

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE)
