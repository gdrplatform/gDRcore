data <- "finalSE_medium.RDS"
original <- get_synthetic_data(data)

df_layout <- merge(cell_lines[1:15, ], drugs[1:40, ], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_merged_data <- gDRtestData::generate_response_data(df_layout)

test_synthetic_data(original, df_merged_data, data)
