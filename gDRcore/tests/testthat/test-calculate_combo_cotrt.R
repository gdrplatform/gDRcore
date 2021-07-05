library(testthat)
library(gDRcore)
se <- readRDS(system.file("testdata", 
                          "finalSE_combo_codilution_small.RDS", 
                          package = "gDRtestData"))

n <- 4
drug_names <- paste0("drug_00", seq(n))
printed_output <- cat(paste0("[1] \"", drug_names[2:4], " x ", drug_names[1], " (", 
                             seq_along(drug_names[2:4]), "/", length(drug_names[2:4]), 
                             ")\"\n"), 
                      sep = "")
test_that("check printed output",
          expect_output(calculate_combo_cotrt(se), printed_output))
