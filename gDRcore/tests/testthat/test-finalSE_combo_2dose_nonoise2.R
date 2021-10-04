data <- "finalSE_combo_2dose_nonoise2.RDS"
original <- get_synthetic_data(data)

df_layout <- merge(cell_lines[2:4, ], drugs[c(2:4, 26), ], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_2 <- cbind(drugs[26 * c(1, 1, 1), ], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)
df_merged_data <- df_merged_data[!(df_merged_data$Gnumber %in% c("vehicle", drugs$Gnumber[26]) & 
                                     df_merged_data$Gnumber_2 == drugs$Gnumber[26]), ]

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, data, combo = TRUE)
