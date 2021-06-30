library(testthat)
library(gDRcore)

test_that("identify_keys works", {
  d_id1 <- c("Gnumber",
    "DrugName",
    "Concentration")
  d_ids <- c("Gnumber_2",
    "DrugName_2",
    "Concentration_2",
    "Gnumber_3",
    "DrugName_3",
    "Concentration_3")
  duration <- "Duration"
  cl <- c("clid",
    "CellLineName",
    "UserCellLineDesignation",
    "Tissue")
  misc <- c("Medium",
    "E2",
    "Replicate",
    "Barcode")

  cols <- c(
    cl,
    d_id1,
    d_ids,
    duration,
    misc,
    "masked",
    "ReadoutValue",
    "BackgroundValue",
    "ReferenceDivisionTime",
    "Template",
    "WellRow",
    "WellColumn",
    "CorrectedReadout"
  )

  df_ <- data.frame(matrix(0, nrow = 1, ncol = length(cols)))
  colnames(df_) <- cols
  df_$DrugName <- gDRutils::get_identifier("untreated_tag")[1] 

  k1 <- identify_keys(df_, nested_keys = NULL)
  expect_equal(sort(k1[["Trt"]]), sort(c(cl, misc, d_id1, d_ids, duration)))
  expect_equal(sort(k1[["ref_Endpoint"]]), sort(c(cl, misc, d_ids, duration)))
  expect_equal(sort(k1[["untrt_Endpoint"]]), sort(c(cl, misc, duration)))
  expect_equal(sort(k1[["Day0"]]), sort(c(cl, misc)))

  # nested_keys argument works.
  nested_key <- "Barcode"
  k2 <- identify_keys(df_, nested_keys = nested_key)
  expect_equal(k1[!names(k1) %in% c("Trt", "nested_keys")], k2[!names(k2) %in% c("Trt", "nested_keys")])
  expect_equal(setdiff(k1$Trt, k2$Trt), nested_key)

  # Remove NA keys.
  df_$E2 <- NA
  k3 <- identify_keys(df_, nested_keys = NULL)
  expect_equal(setdiff(k1[["untrt_Endpoint"]], k3[["untrt_Endpoint"]]), "E2")
})
