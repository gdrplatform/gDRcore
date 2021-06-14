context("Test finalSE_medium")
original <- readRDS(system.file("testdata", "finalSE_medium.RDS", package = "gDRtestData"))
normalized <- gDRutils::convert_se_assay_to_dt(original, "Normalized")
averaged <- gDRutils::convert_se_assay_to_dt(original, "Averaged")
metrics <- gDRutils::convert_se_assay_to_dt(original, "Metrics")

source('gDRtestData/gDRtestData/inst/scripts/functions_generate_data.R')
df_layout <- merge(CellLines[1:15,], Drugs[1:40,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)

se <- gDRcore::create_SE(df_merged_data, override_untrt_controls = NULL)
normSE <- gDRcore::normalize_SE(se)  
avgSE <- gDRcore::average_SE(normSE)
metricsSE <- gDRcore::fit_SE(avgSE)
finalSE <- gDRcore::add_codrug_group_SE(metricsSE)

normalized_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Normalized")
averaged_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Averaged")
metrics_new <- gDRutils::convert_se_assay_to_dt(finalSE, "Metrics")

test_that("Original finalSE_medium data and recreated data are identical", {
  expect_identical(normalized, normalized_new)
  expect_identical(averaged, averaged_new)
  expect_identical(metrics, metrics_new)
})

