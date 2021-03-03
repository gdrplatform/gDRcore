library(testthat); library(gDR);

test_that("calculate_GR_value throws expected errors", {
  readouts <- c(rep(10000, 5))
  expect_error(calculate_GR_value(rel_viability = readouts, 
                     corrected_readout = readouts, 
                     day0_readout = readouts[-1],
                     untrt_readout = readouts,
                     ndigit_rounding = 4,
                     duration = duration, 
                     ref_div_time = duration/2)
  )
})


test_that("calculate_GR_value works as expected", {
  duration <- 144
  rv <- seq(0.1, 1, 0.1)
  corrected <- seq(41000, 50000, 1000)
  day0 <- seq(91000, 95500, 500)
  untrt <- rep(c(115000, 118000), 5)

  # No cell line is specified.
  gr1 <- calculate_GR_value(rel_viability = rv, 
    corrected_readout = corrected, 
    day0_readout = day0,
    untrt_readout = untrt,
    ndigit_rounding = 4,
    duration = duration, 
    ref_div_time = duration/2)
  expect_true(is.numeric(gr1))

  # day0_readout is NA, reference division time missing.
  expect_warning(gr2 <- calculate_GR_value(rel_viability = rv, 
    corrected_readout = corrected, 
    day0_readout = NA, 
    untrt_readout = untrt, 
    ndigit_rounding = 4,
    duration = duration, 
    ref_div_time = NA,
    cl_name = "CELLLINE_NAME"), regexp = "GR values are NA")
  expect_true(all(is.na(gr2)))

  # day0_readout is NA, reference division time is present and valid.
  expect_warning(gr3 <- calculate_GR_value(rel_viability = rv, 
    corrected_readout = corrected, 
    day0_readout = NA, 
    untrt_readout = untrt, 
    ndigit_rounding = 4,
    duration = duration, 
    ref_div_time = (duration * 1.5) - 1,
    cl_name = "CELLLINE_NAME"), regexp = "calculating GR value based on reference doubling time")
  expect_true(all(gr1 != gr3))

  # day0_readout is NA, reference division time is present but not valid.
  expect_warning(gr4 <- calculate_GR_value(rel_viability = rv, 
    corrected_readout = corrected, 
    day0_readout = NA, 
    untrt_readout = untrt, 
    ndigit_rounding = 4,
    duration = duration, 
    ref_div_time = (duration * 1.5) + 1,
    cl_name = "CELLLINE_NAME"), regexp = "too long for GR calculation")
  expect_true(all(is.na(gr4)))
})
