test_that("add_CellLine_annotation works", {
  df <- data.frame(ReadouValue = runif(5),
                   CLID = paste0("CL", 1:5))
  df_annotated <- add_CellLine_annotation(df)
  expect_true(all(is.na(df_annotated$ReferenceDivisionTime)))
  expect_equal(df_annotated$CLID, df_annotated$CellLineName)
})
