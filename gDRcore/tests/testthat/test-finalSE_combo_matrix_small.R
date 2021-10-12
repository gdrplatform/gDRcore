data <- "finalSE_combo_matrix_small.RDS"
original <- get_synthetic_data(data)

df_layout <- merge(cell_lines[7:8, ], drugs[c(4:6), ], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout, concentrations = 10 ^ (seq(-3, .5, .5)))

df_2 <- merge(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26), ], by = NULL)
df_2 <- gDRtestData::add_data_replicates(df_2)
df_2 <- gDRtestData::add_concentration(df_2, concentrations = 10 ^ (seq(-3, .5, .5)))
colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")] <- 
  paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")], "_2")

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, data, combo = TRUE)
