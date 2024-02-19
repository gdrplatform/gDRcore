test_that("small_no_noise: test_synthetic_data", {
  data <- "finalMAE_small_no_noise"
  original <- gDRutils::get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(generateNoNoiseRawData)(
    cell_lines, drugs, FALSE
  )
  expect_lte(length(mae$warnings), 2)

  test_synthetic_data(original, mae$result, data)
})
