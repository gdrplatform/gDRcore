test_that("calculate_isobolograms caps co-dilution points on the co-dilution line", {
  fit <- data.table::data.table(x_inf = 0.1, x_0 = 1, ec50 = 1, h = 2)
  drug_1 <- cbind(fit, cotrt_value = c(0, 1))
  drug_2 <- cbind(fit, cotrt_value = c(0, 1))
  # at 50% the summed concentration is about 3.35, so both drugs are above their caps of 0.8
  codilution <- data.table::data.table(x_inf = c(0.1, 0.1), x_0 = 1, ec50 = 3, h = 2, ratio = c(2, 0.5))

  iso <- calculate_isobolograms(drug_2, drug_1, codilution, isobol_value = 0.5, max1_cap = 0.8, max2_cap = 0.8)
  codil <- iso[iso$fit_type == "by_codil"]

  expect_false(anyNA(codil$conc_1))
  expect_false(anyNA(codil$conc_2))
  # ratio = conc_2 / conc_1, and the point stops at whichever cap it reaches first
  expect_equal(codil$conc_2 / codil$conc_1, codilution$ratio)
  expect_equal(codil$conc_1, c(0.4, 0.8))
  expect_equal(codil$conc_2, c(0.8, 0.4))
})
