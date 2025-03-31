test_that("get_cell_line_annotation works correctly", {
  data <- data.table::data.table(clid = c("CL1", "CL2", "CL3"))
  # Assuming the annotation file "cell_lines.csv" is available in the package "gDRtestData"
  result <- get_cell_line_annotation(data, fill = "unknown", annotation_package = "gDRtestData")
  expect_true("data.table" %in% class(result))
  expect_equal(nrow(result), 3)
  expect_equal(result$clid, c("CL1", "CL2", "CL3"))
  # Check if the fill value is correctly applied to missing annotations
  expect_equal(result$Tissue[result$clid == "CL3"], "unknown")
})

test_that("annotate_dt_with_cell_line works correctly", {
  data <- data.table::data.table(clid = c("CL1", "CL2", "CL3"))
  cell_line_annotation <- data.table::data.table(
    clid = c("CL1", "CL2"),
    CellLineName = c("Cell Line 1", "Cell Line 2"),
    Tissue = c("Tissue 1", "Tissue 2"),
    ReferenceDivisionTime = c(24, 48),
    parental_identifier = c("Parent 1", "Parent 2"),
    subtype = c("Subtype 1", "Subtype 2")
  )
  result <- annotate_dt_with_cell_line(data, cell_line_annotation, fill = "unknown")
  expect_true(all(c(24, 48) %in% result$ReferenceDivisionTime))
  expect_true("data.table" %in% class(result))
  expect_equal(ncol(result), ncol(data) + ncol(cell_line_annotation) - 1)
  expect_equal(result$CellLineName, c("Cell Line 1", "Cell Line 2", NA))
})

test_that("get_drug_annotation works correctly", {
  data <- data.table::data.table(Gnumber = c("drug1", "drug2", "drug3"))
  # Assuming the annotation file "drugs.csv" is available in the package "gDRtestData"
  result <- get_drug_annotation(data, fill = "unknown", annotation_package = "gDRtestData")
  expect_true("data.table" %in% class(result))
  expect_equal(nrow(result), 3)
  expect_equal(result$Gnumber, c("drug1", "drug2", "drug3"))
  # Check if the fill value is correctly applied to missing annotations
  expect_equal(result$drug_moa[result$Gnumber == "drug3"], "unknown")
  
  complex_data <- data.table::data.table(Gnumber = c("D1", "D2", "D3"), Gnumber_2 = c("D4", "D5", "D6"))
  result <- get_drug_annotation(complex_data, fill = "unknown", annotation_package = "gDRtestData")
  expect_true("data.table" %in% class(result))
  expect_equal(nrow(result), 6)
})

test_that("annotate_dt_with_drug works correctly", {
  data <- data.table::data.table(Gnumber = c("D1", "D2", "D3"))
  drug_annotation <- data.table::data.table(
    Gnumber = c("D1", "D2"),
    DrugName = c("Drug 1", "Drug 2"),
    drug_moa = c("MOA 1", "MOA 2")
  )
  result <- annotate_dt_with_drug(data, drug_annotation, fill = "unknown")
  expect_true("data.table" %in% class(result))
  expect_equal(ncol(result), ncol(data) + ncol(drug_annotation) - 1)
  expect_equal(result$DrugName, c("Drug 1", "Drug 2", "D3"))
  
  complex_data <- data.table::data.table(Gnumber = c("D1", "D2", "D3"), Gnumber_2 = c("D4", "D5", "D6"))
  complex_drug_annotation <- data.table::data.table(
    Gnumber = c("D1", "D2", "D4", "D5"),
    DrugName = c("Drug 1", "Drug 2", "Drug 4", "Drug 5"),
    drug_moa = c("MOA 1", "MOA 2", "MOA 4", "MOA 5")
  )
  result <- annotate_dt_with_drug(complex_data, complex_drug_annotation, fill = "unknown")
  expect_true("data.table" %in% class(result))
  expect_equal(result$DrugName, c("Drug 1", "Drug 2", "D3"))
  expect_equal(result$DrugName_2, c("Drug 4", "Drug 5", "D6"))
})

test_that("get_drug_annotation_from_dt works as expected", {
  dt_example <- data.table::data.table(Gnumber = "drug_id",
                                       DrugName = "DrugName",
                                       drug_moa = "drug_moa",
                                       some_col = "value")
  annotation <- get_drug_annotation_from_dt(dt_example)
  testthat::expect_true(data.table::is.data.table(annotation))
  testthat::expect_equal(dim(annotation), c(1, 3))
  dt_example$drug_moa <- NULL
  annotation2 <- get_drug_annotation_from_dt(dt_example)
  expect_equal(annotation2$drug_moa, "unknown")
})

test_that("get_cellline_annotation_from_dt works as expected", {
  dt_example <- data.table::data.table(ReadoutValue = runif(5),
                                       clid = paste0("CL", 1:5),
                                       CellLineName = paste0("RandomName", 1:5),
                                       Tissue =  paste0("Tissue", 1:5),
                                       ReferenceDivisionTime = 1:5,
                                       parental_identifier = 1:5,
                                       subtype = "subtype",
                                       some_col = "value")
  annotation <- get_cellline_annotation_from_dt(dt_example)
  testthat::expect_true(data.table::is.data.table(annotation))
  testthat::expect_equal(dim(annotation), c(5, 6))
  dt_example$subtype <- NULL
  annotation2 <- get_cellline_annotation_from_dt(dt_example)
  expect_true(all(annotation2$subtype == "unknown"))
})

