#library(testthat); library(gDR)
source("setUp.R")

test_that("splitSEComponents splits the correct columns", {
  # Standard case.
  md <- splitSEComponents(test_df)
  expect_true(all(c("Gnumber", "DrugName", "replicates", "drug_moa") %in% colnames(md$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md$condition_md)))
  expect_equal(ncol(test_df), 
    sum(ncol(md$treatment_md), ncol(md$condition_md), length(md$data_fields), ncol(md$experiment_md))) 

  # Check that discard_keys argument works as expected.
  md2 <- splitSEComponents(test_df, discard_keys = c("replicates"))
  expect_true(all(c("Gnumber", "DrugName", "drug_moa") %in% colnames(md2$treatment_md)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md2$condition_md)))
  expect_true(all(c("WellRow", "WellColumn", "Concentration", "replicates") %in% md2$data_fields))
  expect_equal(ncol(test_df), 
    sum(ncol(md2$treatment_md), ncol(md2$condition_md), length(md2$data_fields), ncol(md2$experiment_md)))
})


test_that("splitSEComponents throws a warning for bad cell line metadata", {
  df2 <- test_df
  df2$CellLineName[7] <- "Abnormality"
  expect_warning(splitSEComponents(df2), 
    regexp = "'CellLineName' not metadata for unique cell line identifier column")
})
