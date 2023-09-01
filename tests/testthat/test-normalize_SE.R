test_that("normalize_SE works as expected", {
  # Set up. 
  conc <- rep(seq(0.1, 0.3, 0.1), 2)
  ctrl_df <- S4Vectors::DataFrame(Barcode = c(1, 2),
                                  CorrectedReadout = c(rep(1, 2),
                                                       rep(6, 2)),
                                  control_type = c(rep("UntrtReadout", 2),
                                                   rep("Day0Readout", 2)),
                                  isDay0 = c(rep(FALSE, 2),
                                             rep(TRUE, 2)))
  
  trt_df <- S4Vectors::DataFrame(CorrectedReadout = rep(seq(1, 3, 1), 2),
                                 Concentration = conc,
                                 Barcode = rep(c(1, 2), each = 3),
                                 masked = rep(FALSE, 6))

  coldata <- S4Vectors::DataFrame(CellLineName = "Mickey Mouse",
                                  ReferenceDivisionTime = 1)
  rowdata <- S4Vectors::DataFrame(Duration = 2)

  ctrl <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = ctrl_df)
  trted <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = trt_df)
  
  keys <- list("nested_keys" = "Barcode", 
               "Trt" = "Concentration")
  
  metadata <- list(identifiers = list("cellline_name" = "CellLineName", 
                                      "cellline_ref_div_time" = "ReferenceDivisionTime", 
                                      "duration" = "Duration", 
                                      "masked_tag" = "masked",
                                      "barcode" = "Barcode",
                                      "concentration" = "Concentration",
                                      "concentration2" = "Concentration2"),
                   Keys = keys)
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = list("RawTreated" = trted, "Controls" = ctrl), 
                                                   colData = coldata, 
                                                   rowData = rowdata,
                                                   metadata = metadata)


  se <- normalize_SE(se, data_type = "single-agent")
  normalized <- SummarizedExperiment::assays(se)[["Normalized"]][1, 1][[1]]

  expect_true(methods::is(normalized, "DataFrame"))
  expect_equal(dim(normalized), c(12, 4))
  expect_true(all(colnames(normalized) %in% c("Concentration", "masked", "normalization_type", "x")))
  expect_equal(unique(normalized$Concentration), unique(conc))
  
  se2 <- normalize_SE(se, "single-agent", nested_confounders = 
                        c(gDRutils::get_SE_identifiers(se, "barcode", simplify = TRUE), "masked"))
  expect_s4_class(se2, "SummarizedExperiment")
  
  
  ctrl_df$Barcode <- NULL
  trt_df$Barcode <- NULL
  ctrl2 <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = ctrl_df)
  trted2 <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = trt_df)
  se3 <- SummarizedExperiment::SummarizedExperiment(assays = list("RawTreated" = trted2, "Controls" = ctrl2), 
                                                   colData = coldata, 
                                                   rowData = rowdata,
                                                   metadata = metadata)
  
  se3 <- normalize_SE(se3, "single-agent", nested_confounders = NULL, "masked")
  expect_s4_class(se3, "SummarizedExperiment")
})


test_that("normalize_SE works as expected with Day0data", {
  
  td <- gDRimport::get_test_data()
  l_tbl <- purrr::quietly(gDRimport::load_data)(gDRimport::manifest_path(td), 
                                                gDRimport::template_path(td), 
                                                gDRimport::result_path(td))
  imported_data <- purrr::quietly(merge_data)(
    data.table::setDT(l_tbl$result$manifest),
    data.table::setDT(l_tbl$result$treatments),
    data.table::setDT(l_tbl$result$data)
  )
  
  data_type <- "single-agent"
  se <- purrr::quietly(create_SE)(imported_data$result, data_type = data_type, nested_confounders = NULL)
  se <- purrr::quietly(normalize_SE)(se$result, data_type = data_type)
  
  
  # Check Day0 data
  norm <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se[[1]], "Normalized"))
  expect_true(all(is.na(controls$CorrectedReadout[
    controls$control_type == "Day0Readout"])))
})