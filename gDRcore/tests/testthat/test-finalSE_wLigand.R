
skip(message = "Ligand")

data <- "finalSE_wLigand.RDS"
original <- get_synthetic_data(data)

df_layout <- merge(cell_lines[2:6, ], drugs[2:5, ], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_merged_data <- gDRtestData::generate_response_data(df_layout, 0)
df_merged_data$Ligand <- 0.1
df_merged_data2 <- df_merged_data[df_merged_data$Gnumber %in% c("vehicle", "G00002", "G00003"), ]
df_merged_data2$Ligand <- 0
df_merged_data2$ReadoutValue <- 105 - pmax(0, pmin(104, (105 - df_merged_data2$ReadoutValue) ** 1.1))
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 11:12)] <- 
  0.8 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 11:12)]
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 13:14)] <- 
  0.5 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 13:14)]
df_merged_data2$ReadoutValue <- round(df_merged_data2$ReadoutValue, 1)
df_merged_data2$Barcode <- paste0(df_merged_data2$Barcode, "1")
df_merged_data <- rbind(df_merged_data, df_merged_data2)


finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = c(Ligand = 0.1))

test_synthetic_data(original, finalSE, data, additional_columns = 1)
