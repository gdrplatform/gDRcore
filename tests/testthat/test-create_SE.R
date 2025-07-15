test_that("validate_mapping works as expected", {
  ref_df <- data.table::data.table(Gnumber = paste0("DRUG_", c(rep(10, 3), rep(11, 7))),
                       Gnumber2 = "vehicle", value = runif(10))
  trt_df <- data.table::data.table(Gnumber = paste0("DRUG_", 10),
                       Gnumber2 = paste0("DRUG_", 1:10),
                       value = runif(10))
  merged_df <- validate_mapping(trt_df, ref_df, nested_confounders = NULL)
  expect_equal(dim(merged_df), c(13, 3))  
  
  ref_df$Barcode <- c("A", rep("B", 9))
  trt_df$Barcode <- "A"
  merged_df2 <- validate_mapping(trt_df, ref_df, nested_confounders = "Barcode")
  expect_equal(dim(merged_df2), c(11, 4))
})

test_that("validate_mapping catches reverse single-agent data", {
  trt_df <- data.table::data.table(Gnumber = "DRUG_10",
                       Gnumber_2 = "DRUG_2", value = runif(1),
                       Concentration = runif(1),
                       Concentration_2 = runif(1))
  ref_df <- data.table::data.table(Gnumber = c("DRUG_10", "DRUG_2", "DRUG_3"),
                       Gnumber_2 = rep("vehicle", 3), value = runif(3),
                       Concentration = runif(3), Concentration_2 = rep(0, 3))
  merged_df <- validate_mapping(trt_df, ref_df, nested_confounders = NULL)
  expect_equal(dim(merged_df), c(3, 6))  
})



test_that("create_SE works as expected", {
  
  td <- gDRimport::get_test_data()
  l_tbl <- purrr::quietly(gDRimport::load_data)(gDRimport::manifest_path(td), 
                                                gDRimport::template_path(td), 
                                                gDRimport::result_path(td))
  imported_data <- purrr::quietly(merge_data)(
    data.table::setDT(l_tbl$result$manifest),
    data.table::setDT(l_tbl$result$treatments),
    data.table::setDT(l_tbl$result$data)
  )

  se <- purrr::quietly(create_SE)(imported_data$result, data_type = "single-agent", nested_confounders = NULL)

  testthat::expect_s4_class(se$result, "SummarizedExperiment")
  expect_equal(dim(se$result), c(12, 6))
  
  # Check Day0 data
  controls <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se[[1]], "Controls"))
  expect_true(all(is.na(controls$CorrectedReadout[
    controls$control_type == "Day0Readout"])))
})


test_that("create_SE works with empty nested confounder", {
  
  td <- gDRimport::get_test_data()
  l_tbl <- purrr::quietly(gDRimport::load_data)(gDRimport::manifest_path(td), 
                                                gDRimport::template_path(td), 
                                                gDRimport::result_path(td))
  imported_data <- purrr::quietly(merge_data)(
    data.table::setDT(l_tbl$result$manifest),
    data.table::setDT(l_tbl$result$treatments),
    data.table::setDT(l_tbl$result$data)
  )
  
  imported_data$result$Barcode <- NULL
  
  se <- purrr::quietly(create_SE)(imported_data$result, data_type = "single-agent", nested_confounders = NULL)
  
  testthat::expect_s4_class(se$result, "SummarizedExperiment")
  expect_equal(dim(se$result), c(2, 6))
})

test_that("create_SE swap drugs properly", {
  
  data_mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  data_raw <- gDRcore::convert_mae_to_raw_data(data_mae)
  data_new <- data_raw[data_raw$Gnumber == "G00005" & data_raw$Gnumber_2 == "G00026"]
  
  data_new$Gnumber_2 <- "G00005"
  data_new$DrugName_2 <- "drug_005"
  data_new$drug_moa_2 <- "moa_A"
  
  data_new$Gnumber <- "G00026"
  data_new$DrugName <- "drug_026"
  data_new$drug_moa <- "moa_E"
  data_final <- rbind(data_raw, data_new)
  mae <- gDRcore::create_SE(data_final, data_type = "combination")
  mae2 <- gDRcore::create_SE(data_raw, data_type = "combination")
  expect_equal(dim(mae), dim(mae2))
})

