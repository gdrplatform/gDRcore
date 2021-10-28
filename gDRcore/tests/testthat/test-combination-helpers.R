test_that(".calculate_dilution_ratio works as expected", {
  ratio <- 0.5
  concs <- 10 ^ (seq(-3, 1, ratio))
  obs <- gDRcore:::.calculate_dilution_ratio(concs)
  expect_equal(obs, ratio)
})


test_that("map_conc_to_standardized_conc works as expected", {
  ratio <- 0.5
  conc1 <- c(0, 10 ^ (seq(-3, 1, ratio)))

  shorter_range <- conc1[-1]
  noise <- runif(length(shorter_range), 1e-12, 1e-11)
  conc2 <- shorter_range + noise

  obs <- map_conc_to_standardized_conc(conc1, conc2)
  expect_true(is(obs, "data.frame"))
})


test_that("replace_conc_with_standardized_conc works as expected", {
  conc_map <- data.frame(orig = c(0.99, 0.6, 0.456, 0.4), std = c(1, 0.6, 0.46, 0.4))
  original_concs <- c(0.456, 0.456, 0.4, 0.99)
  exp <- c(0.46, 0.46, 0.4, 1)
  obs <- replace_conc_with_standardized_conc(original_concs, conc_map, original_conc_col = "orig", standardized_conc_col = "std")
  expect_equal(obs, exp)
})
