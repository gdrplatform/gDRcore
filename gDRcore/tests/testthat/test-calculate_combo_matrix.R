test_that("check printed output", {
  se <- get_synthetic_data("finalSE_combo_matrix_small.RDS")
  cell_line_names <- rep(c("cellLine_GB", "cellline_HB"), 12)
  drug_names1 <- rep(paste0("drug_00", seq(4, 6)), each = 2)
  drug_names2 <- rep(c("drug_021", "drug_026"), 3)
  
  ccm <- suppressWarnings(calculate_combo_matrix(se))
  
  exp_data_names <-
    c("bliss_q10", "CI_100x_80", "CI_100x_50", "hsa_q10")
  expect_equal(names(ccm$agg_results_norm),
               c("RelativeViability", "GRvalue"))
  
  exp_data_cols <- c("cellline_GB CL00016", "cellline_HB CL00017")
  expect_equal(names(ccm$agg_results_norm[[1]]), exp_data_names)
  expect_equal(colnames(ccm$agg_results_norm[[1]][[1]]), exp_data_cols)
  expect_equal(dim(ccm$agg_results_norm[[1]][[1]]), c(6, 2))
  
  rnames <- paste0(drug_names1, " x ", drug_names2, " (T=72)")
  expect_equal(rownames(ccm$agg_results_norm[[1]][[1]]), rnames)
  
})
