combo_assays <- c("Normalized", "Averaged", "Metrics", "excess", "scores")

test_that("combo_2dose_nonoise: test_synthetic_data", {
  data <- "finalMAE_combo_2dose_nonoise"
  original <- gDRutils::get_synthetic_data(data)

  set.seed(2)
  mae <- purrr::quietly(generateComboNoNoiseData)(
    cell_lines, drugs, FALSE
  )

  test_synthetic_data(original, mae$result, data, assays = combo_assays)
})

test_that("combo_matrix_small: test_synthetic_data", {
  data <- "finalMAE_combo_matrix_small"
  original <- gDRutils::get_synthetic_data(data)

  set.seed(2)
  mae <- purrr::quietly(generateComboMatrixSmall)(
    cell_lines, drugs, FALSE
  )

  test_synthetic_data(original, mae$result, data, assays = combo_assays)
})

test_that("combo_triple: test_synthetic_data", {
  data <- "finalMAE_combo_triple"
  original <- gDRutils::get_synthetic_data(data)

  set.seed(2)
  mae <- purrr::quietly(generateTripleComboMatrix)(
    cell_lines, drugs, FALSE
  )

  test_synthetic_data(original, mae$result, data, assays = combo_assays)
})
