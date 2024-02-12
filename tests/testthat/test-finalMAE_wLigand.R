test_that("wLigand: test_synthetic_data", {
  data <- "finalMAE_wLigand"
  original <- gDRutils::get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(generateLigandData)(
    cell_lines, drugs, FALSE
  )
  expect_lte(length(mae$warnings), 2)
  
  test_synthetic_data(original, mae$result, data)
})
