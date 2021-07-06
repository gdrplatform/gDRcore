library(testthat)
library(gDRcore)
se <- readRDS(system.file("testdata", 
                          "finalSE_combo_matrix_small.RDS", 
                          package = "gDRtestData"))

cell_line_names <- rep(c("cellLine_GB", "cellline_HB"), 12)
drug_names1 <- rep(rep(paste0("drug_00", seq(4, 6)), each = 4), 2)
drug_names2 <- rep(rep(c("drug_021", "drug_026"), each = 2), 6)
printed_output <- cat(paste0("[1] \"Calculation for cell line ", cell_line_names,
                             " treated with ", drug_names1,
                             " x ", drug_names2, "\"\n"), sep = "")
test_that("check printed output",
          expect_output(suppressWarnings(calculate_combo_matrix(se)), printed_output))
