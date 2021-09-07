
test_that("check printed output", {
  se <- get_synthetic_data("finalSE_combo_codilution_small.RDS")
  se <- get_synthetic_data("finalSE_combo_2dose_nonoise.RDS")
  
  
  n <- 4
  drug_names <- paste0("drug_00", seq(n))
  printed_output <- cat(paste0("[1] \"", drug_names[2:4], " x ", drug_names[1], " (", 
                               seq_along(drug_names[2:4]), "/", length(drug_names[2:4]), 
                               ")\"\n"), 
                        sep = "")
  expect_output(calculate_combo_cotrt(se), printed_output)
})
