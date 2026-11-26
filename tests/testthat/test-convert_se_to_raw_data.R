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

test_that("convert_se_to_raw_data reads an SE with custom identifiers and keeps the session", {
  se <- gDRutils::get_synthetic_data("finalMAE_small")[[1]]
  ids <- gDRutils::get_SE_identifiers(se)
  ids[c("drug", "cellline", "concentration", "duration")] <- list("compound_id", "cell_id", "dose_uM", "time_h")
  S4Vectors::metadata(se)$identifiers <- ids
  se_custom <- gDRutils::standardize_se(se, use_default = FALSE)
  expect_true("time_h" %in% names(SummarizedExperiment::rowData(se_custom)))

  gDRutils::reset_env_identifiers()
  on.exit(gDRutils::reset_env_identifiers(), add = TRUE)
  raw_default <- convert_se_to_raw_data(gDRutils::standardize_se(se_custom))

  # a session with its own names is the same after the call
  gDRutils::set_env_identifier("drug", "my_drug")
  session <- gDRutils::get_env_identifiers()
  raw_custom <- convert_se_to_raw_data(se_custom)
  expect_identical(gDRutils::get_env_identifiers(), session)
  expect_equal(raw_custom, raw_default)
  expect_true(all(c("Gnumber", "clid", "Concentration", "Duration") %in% names(raw_custom)))
})
