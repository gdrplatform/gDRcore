combo_assays <- c("Normalized", "Averaged", "Metrics", "excess", "scores")
# The combination fits differ between platforms by up to about 0.6% (x_AOC_range
# on Linux vs macOS), more than the default tolerance; the regressions these tests
# guard against (the arm fit, CIScore) change smooth and the scores by 1% to 200%.
combo_tolerance <- 1e-2

test_that("combo_2dose_nonoise: test_synthetic_data", {
  data <- "finalMAE_combo_2dose_nonoise"
  original <- gDRutils::get_synthetic_data(data)

  set.seed(2)
  mae <- purrr::quietly(generateComboNoNoiseData)(
    cell_lines, drugs, FALSE
  )

  test_synthetic_data(original, mae$result, data, assays = combo_assays, tolerance = combo_tolerance)
})

test_that("combo_matrix_small: test_synthetic_data", {
  data <- "finalMAE_combo_matrix_small"
  original <- gDRutils::get_synthetic_data(data)

  set.seed(2)
  mae <- purrr::quietly(generateComboMatrixSmall)(
    cell_lines, drugs, FALSE
  )

  test_synthetic_data(original, mae$result, data, assays = combo_assays, tolerance = combo_tolerance)
})

test_that("combo_triple: test_synthetic_data", {
  data <- "finalMAE_combo_triple"
  original <- gDRutils::get_synthetic_data(data)

  set.seed(2)
  mae <- purrr::quietly(generateTripleComboMatrix)(
    cell_lines, drugs, FALSE
  )

  test_synthetic_data(original, mae$result, data, assays = combo_assays, tolerance = combo_tolerance)
})
