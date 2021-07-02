library(testthat)

se <- readRDS(system.file("testdata", 
                          "finalSE_combo_codilution_small.RDS", 
                          package = "gDRtestData"))
calculate_combo_cotrt(se)

n <- 4
drug_names <- paste0("drug_00", seq(n))
# printed_output <- paste0(drug_names[1], " x ", drug_names[2:n], " (", seq(n-1), "/", n - 1, ")")
printed_output <- "[1] \"drug_002 x drug_001 (1/3)\"\n[1] \"drug_003 x drug_001 (2/3)\"\n[1] \"drug_004 x drug_001 (3/3)"
test_that("check output",
          expect_output(calculate_combo_cotrt(se), printed_output))
