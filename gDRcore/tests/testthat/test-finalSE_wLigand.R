original <- readRDS(system.file("testdata", "finalSE_wLigand.RDS", package = "gDRtestData"))
normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")

source(system.file("scripts", "functions_generate_data.R", package = "gDRtestData"))
df_layout <- merge(CellLines[2:6,], Drugs[2:5,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data$Ligand <- 0.1
df_merged_data2 <- df_merged_data[df_merged_data$Gnumber %in% c('vehicle', 'G00002', 'G00003'),]
df_merged_data2$Ligand <- 0
df_merged_data2$ReadoutValue <- 105 - pmax(0, pmin(104, (105-df_merged_data2$ReadoutValue) ** 1.1))
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 11:12)] <- 
  0.8 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 11:12)]
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 13:14)] <- 
  0.5 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 13:14)]
df_merged_data2$ReadoutValue <- round(df_merged_data2$ReadoutValue,1)
df_merged_data2$Barcode <- paste0(df_merged_data2$Barcode, '1')
df_merged_data <- rbind(df_merged_data, df_merged_data2)

se <- gDRcore::create_SE(df_merged_data, override_untrt_controls = c(Ligand = 0.1))
normSE <- gDRcore::normalize_SE(se)  
avgSE <- gDRcore::average_SE(normSE)
metricsSE <- gDRcore::fit_SE(avgSE)
finalSE <- gDRcore::add_codrug_group_SE(metricsSE)

normalized_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Normalized")
averaged_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Averaged")
metrics_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Metrics")

test_that("Original finalSE_wLigand data and recreated data are identical", {
  expect_identical(normalized, normalized_new)
  expect_identical(averaged, averaged_new)
  expect_identical(metrics, metrics_new)
})

