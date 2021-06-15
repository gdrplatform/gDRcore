original <- readRDS(system.file("testdata", "finalSE_combo_2dose_nonoise3.RDS", package = "gDRtestData"))
normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[2:4,], Drugs[2:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[c(26,26,26),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)
df_layout_2 = df_layout_2[!(df_layout_2$Concentration == 0 & df_layout_2$Concentration_2 > 0), ]

df_merged_data <- generate_response_data(df_layout_2, 0)

se <- gDRcore::create_SE(df_merged_data, override_untrt_controls = NULL)
normSE <- gDRcore::normalize_SE(se)  
avgSE <- gDRcore::average_SE(normSE)
metricsSE <- gDRcore::fit_SE(avgSE)
finalSE <- gDRcore::add_codrug_group_SE(metricsSE)

normalized_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Normalized")
averaged_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Averaged")
metrics_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Metrics")

test_that("Original finalSE_combo_2dose_nonoise3 data and recreated data are identical", {
  expect_identical(normalized, normalized_new)
  expect_identical(averaged, averaged_new)
  expect_identical(metrics, metrics_new)
})

