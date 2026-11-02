test_that("convert_mae_to_raw_data works as expected with sa data", {
  data <- "finalMAE_small"
  original <- gDRutils::get_synthetic_data(data)
  set.seed(2)
  mae <- purrr::quietly(generateNoiseRawData)(
    cell_lines, drugs, FALSE
  )
  input_df <- convert_mae_to_raw_data(mae$result)
  mae2 <- purrr::quietly(runDrugResponseProcessingPipeline)(
    input_df
  )
  test_synthetic_data(original, mae$result, data)
  test_synthetic_data(mae$result, mae2$result, data)
})

test_that("convert_mae_to_raw_data works as expected with matrix data", {
  data <- "finalMAE_combo_2dose_nonoise"
  original <- gDRutils::get_synthetic_data(data)
  set.seed(2)
  mae <- purrr::quietly(generateComboNoNoiseData)(
    cell_lines, drugs, FALSE
  )
  input_df <- convert_mae_to_raw_data(mae$result)
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  input_df <- data.table::as.data.table(lapply(input_df, function(x) {
    ifelse(x %in% untreated_tag, untreated_tag[2], x)
  }))
  mae2 <- purrr::quietly(runDrugResponseProcessingPipeline)(
    input_df
  )
  test_synthetic_data(original, mae$result, data)
  test_synthetic_data(mae$result, mae2$result, data)

  trt <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(original[[1]], "RawTreated"))
  trt <- trt[!duplicated(trt$record_id), ]
  dt_raw1 <- convert_mae_to_raw_data(original)
  SummarizedExperiment::assay(original[[1]], "RawTreated") <- BumpyMatrix::splitAsBumpyMatrix(
    trt[!colnames(trt) %in% c("row", "column")],
    row = trt$row,
    column = trt$column)
  dt_raw2 <- convert_mae_to_raw_data(original)
  expect_true("G00026" %in% dt_raw1$Gnumber)
  expect_true(!"G00026" %in% dt_raw2$Gnumber)
})


test_that("convert_mae_to_raw_data keeps the data of an MAE with one experiment", {
  cd_name <- gDRutils::get_supported_experiments("cd")
  set.seed(2)
  mae <- purrr::quietly(generateCodilutionSmall)(cell_lines, drugs, save = FALSE)$result
  raw <- convert_mae_to_raw_data(mae)
  expect_gt(NROW(raw), 0L)
  # reprocessing the raw data gives the same co-dilution fits
  mae2 <- purrr::quietly(runDrugResponseProcessingPipeline)(raw, nested_confounders = "Barcode")$result
  cols <- c("DrugName", "CellLineName", "normalization_type", "ratio", "fit_type", "xc50")
  m1 <- gDRutils::convert_se_assay_to_dt(mae[[cd_name]], "Metrics")[, cols, with = FALSE]
  m2 <- gDRutils::convert_se_assay_to_dt(mae2[[cd_name]], "Metrics")[, cols, with = FALSE]
  data.table::setorder(m1)
  data.table::setorder(m2)
  expect_equal(m2, m1)

  combo <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  combo_only <- purrr::quietly(function() combo[, , gDRutils::get_supported_experiments("combo")])()$result
  expect_gt(NROW(convert_mae_to_raw_data(combo_only)), 0L)
})

test_that("convert_mae_to_raw_data works for co-dilution and single-agent data", {
  set.seed(2)
  df <- gDRtestData::prepareCodilutionData(cell_lines[seq_len(2), ], drugs[seq_len(3), ], noise = 0)
  conc2 <- gDRutils::get_env_identifiers("concentration2")
  # single-agent arm of drug 1 on extra plates
  sa <- df[df[[conc2]] > 0]
  sa[[conc2]] <- 0
  sa$Gnumber_2 <- sa$DrugName_2 <- sa$drug_moa_2 <- gDRutils::get_env_identifiers("untreated_tag")[1]
  sa$Barcode <- paste0(sa$Barcode, "_sa")
  mae <- purrr::quietly(runDrugResponseProcessingPipeline)(rbind(df, sa), nested_confounders = "Barcode")$result
  expect_setequal(names(mae), c(gDRutils::get_supported_experiments("cd"), gDRutils::get_supported_experiments("sa")))

  raw <- convert_mae_to_raw_data(mae)
  expect_gt(NROW(raw), 0L)
  mae2 <- purrr::quietly(runDrugResponseProcessingPipeline)(raw, nested_confounders = "Barcode")$result
  expect_setequal(names(mae2), names(mae))
})
