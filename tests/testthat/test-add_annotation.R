test_that("add_CellLine_annotation works", {
  dt_unknown <- data.table::data.table(ReadouValue = runif(5),
                                       clid = paste0("CL", 1:5))
  dt_unknown_annotated <- add_CellLine_annotation(dt_unknown)
  expect_true(all(is.na(dt_unknown_annotated$ReferenceDivisionTime)))
  expect_equal(dt_unknown_annotated$CellLineName, dt_unknown$clid)
  
  dt <- data.table::data.table(ReadouValue = runif(5),
                               clid = paste0("CL000", 11:15))
  CLs_info <- data.table::fread(
    system.file("annotation_data", "cell_lines.csv", package = "gDRtestData")
  )[cell_line_identifier %in% dt$clid]
  dt_annotated <- add_CellLine_annotation(dt, annotationPackage = "gDRtestData")
  
  expect_equal(dt_annotated$CellLineName, CLs_info$cell_line_name)
  expect_equal(dt_annotated$Tissue, CLs_info$primary_tissue)
  expect_true(all(dt_annotated$subtype == "unknown"))
})

test_that("add_Drug_annotation works", {
  dt_unknown <- data.table::data.table(Gnumber = "drug_id")
  dt_unknown_annotated <- add_Drug_annotation(dt_unknown)
  expect_equal(dt_unknown_annotated$drug_moa, "unknown")
  expect_equal(dt_unknown_annotated$Gnumber, dt_unknown_annotated$DrugName)
  
  dt <- data.table::data.table(Gnumber = sprintf("G0000%s", 4:8),
                               DrugName = sprintf("drug_00%s", 4:8))
  Drug_info <- data.table::fread(
    system.file("annotation_data", "drugs.csv", package = "gDRtestData"), header = TRUE,
    col.names = c("drug", "drug_name", "drug_moa"))[drug %in% dt$Gnumber]
  dt_result <- data.table::setnames(Drug_info, c("drug", "drug_name"), c("Gnumber", "DrugName"))
  dt_annotated <- add_Drug_annotation(dt, annotationPackage = "gDRtestData")
  
  expect_identical(dt_annotated, dt_result)
})

test_that("remove_drug_batch works", {
  dt <- data.table::data.table(Gnumber = "DRUG.123")
  gnumber_without_batch <- remove_drug_batch(dt$Gnumber)
  expect_equal(gnumber_without_batch, "DRUG")
})
