test_that("add_CellLine_annotation works", {
  df <- data.table::data.table(ReadouValue = runif(5),
                               clid = paste0("CL", 1:5))
  df_annotated <- add_CellLine_annotation(df)
  expect_true(all(is.na(df_annotated$ReferenceDivisionTime)))
  expect_equal(df_annotated$clid, df_annotated$CellLineName)
})

test_that("add_Drug_annotation works", {
  df <- data.table::data.table(Gnumber = "drug_id")
  df_annotated <- add_Drug_annotation(df)
  expect_equal(df_annotated$drug_moa, "unknown")
  expect_equal(df_annotated$Gnumber, df_annotated$DrugName)
})

test_that("remove_drug_batch works", {
  df <- data.table::data.table(Gnumber = "DRUG.123")
  gnumber_without_batch <- remove_drug_batch(df$Gnumber)
  expect_equal(gnumber_without_batch, "DRUG")
})
