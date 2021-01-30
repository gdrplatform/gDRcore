#library(testthat); library(gDR)
source("setUp.R")

test_that("getMetaData splits the correct columns", {
  # Standard case.
  md <- getMetaData(test_df)
  expect_true(all(c("Gnumber", "DrugName", "replicates", "drug_moa", "row_id", "name_") %in% colnames(md$rowData)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md$colData)))
  expect_equal(ncol(test_df), 
    sum(ncol(md$rowData), ncol(md$colData), length(md$dataCols), ncol(md$csteData)) - 4) 
    # (- 4) because of appended 'name_', 'row_id', and 'col_id'.

  # Check that discard_keys argument works as expected.
  md2 <- getMetaData(test_df, discard_keys = c("replicates"))
  expect_true(all(c("Gnumber", "DrugName", "drug_moa", "row_id", "name_") %in% colnames(md2$rowData)))
  expect_true(all(c("clid", "CellLineName", "Tissue", "ReferenceDivisionTime") %in% colnames(md2$colData)))
  expect_true(all(c("WellRow", "WellColumn", "Concentration", "replicates") %in% md2$dataCols))
  expect_equal(ncol(test_df), 
    sum(ncol(md2$rowData), ncol(md2$colData), length(md2$dataCols), ncol(md2$csteData)) - 4)
})

test_that("getMetaData throws a warning for bad cell line metadata", {
  df2 <- test_df
  df2$CellLineName[7] <- "Abnormality"
  expect_warning(getMetaData(df2), 
    regexp = "'CellLineName' not metadata for unique cell line identifier column")
})
