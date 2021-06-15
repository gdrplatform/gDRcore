original <- readRDS(system.file("testdata", "finalSE_small_no_noise.RDS", package = "gDRtestData"))

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)

finalSE <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = NULL)

test_synthetic_data(original, finalSE)
