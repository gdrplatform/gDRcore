test_that("convert_se_to_raw_data works as expected with sa data", {
  data <- "finalMAE_small.RDS"
  original <- gDRutils::get_synthetic_data(data)
  set.seed(2)
  mae <- purrr::quietly(generateNoiseRawData)(
    cell_lines, drugs, FALSE
  )
  input_df <- convert_mae_to_raw_data(mae$result)
  mae2 <- purrr::quietly(runDrugResponseProcessingPipeline)(
    as.data.frame(input_df)
  )
  test_synthetic_data(original, mae$result, data)
  test_synthetic_data(mae$result, mae2$result, data)
})

test_that("convert_se_to_raw_data works as expected with matrix data", {
  data <- "finalMAE_combo_matrix_small.RDS"
  original <- gDRutils::get_synthetic_data(data)
  set.seed(2)
  mae <- purrr::quietly(generateComboMatrixSmall)(
    cell_lines, drugs, FALSE
  )
  input_df <- convert_mae_to_raw_data(mae$result)
  mae2 <- purrr::quietly(runDrugResponseProcessingPipeline)(
    as.data.frame(input_df)
  )
  test_synthetic_data(original, mae$result, data)
  test_synthetic_data(mae$result, mae2$result, data)
})
