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
})
