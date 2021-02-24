library(testthat); library(gDR);

test_that("calculate_GR_value throws expected warnings", {
  calculate_GR_value()
})


test_that("calculate_GR_value works as expected", {
  duration = 144

  # No cell line is specified.
  calculate_GR_value(rel_viability = , 
                     corrected_readout = , 
                     day0_readout = ,
                     untrt_readout = ,
                     ndigit_rounding = 4,
                     duration = duration, 
                     ref_div_time = )

  # day0_readout is NA, reference division time missing for some cell lines.
  calculate_GR_value(rel_viability =, 
                     day0_readout = NA, 
                     untrt_readout = , 
                     ndigit_rounding = 4,
                     duration = duration, 
                     ref_div_time = NA,
                     cl_name = "CELLLINE_NAME")

  # day0_readout is NA, reference division time missing for all.
  calculate_GR_value(rel_viability =, 
                     day0_readout = NA, 
                     untrt_readout = , 
                     ndigit_rounding = 4,
                     duration = duration, 
                     ref_div_time = NULL,
                     cl_name = "CELLLINE_NAME")

  # day0_readout is NA, reference division time is present and valid.
  calculate_GR_value(rel_viability =, 
                     day0_readout = NA, 
                     untrt_readout = , 
                     ndigit_rounding = 4,
                     duration = duration, 
                     ref_div_time = (duration * 1.5) - 1,
                     cl_name = "CELLLINE_NAME")

  # day0_readout is NA, reference division time is present but not valid.
  calculate_GR_value(rel_viability =, 
                     day0_readout = NA, 
                     untrt_readout = , 
                     ndigit_rounding = 4,
                     duration = duration, 
                     ref_div_time = (duration * 1.5) + 1,
                     cl_name = "CELLLINE_NAME")
})
