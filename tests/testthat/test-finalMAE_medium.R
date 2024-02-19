test_that("medium: test_synthetic_data", {
  data <- "finalMAE_medium"
  original <- gDRutils::get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(generateMediumData)(
    cell_lines, drugs, FALSE
  )
  expect_lte(length(mae$warnings), 2)
  
  test_synthetic_data(original, mae$result, data)
})
