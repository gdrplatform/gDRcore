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
