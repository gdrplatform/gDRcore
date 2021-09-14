
skip(message = "Will be fixed in https://jira.gene.com/jira/browse/GDR-1090")

skip("[FIX IT] skipped due to invalid IC50 values on R 4.1/BioC 3.13")
data <- "finalSE_combo_2dose_nonoise3.RDS"
original <- get_synthetic_data(data)

df_layout <- merge(cell_lines[2:4, ], drugs[2:4, ], by = NULL)
df_layout <- gDRtestData::add_data_replicates(df_layout)
df_layout <- gDRtestData::add_concentration(df_layout)

df_2 <- cbind(drugs[c(26, 26, 26), ], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- merge(df_layout, df_2, by = NULL)
df_layout_2 <- df_layout_2[!(df_layout_2$Concentration == 0 & df_layout_2$Concentration_2 > 0), ]

df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE, data)
