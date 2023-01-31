test_that("fit_SE errors as expected", {
  se <- SummarizedExperiment::SummarizedExperiment()
  expect_error(fit_SE(se))
})

test_that("fit_SE works as expected", {
 
  # combo data 
  fmae_cms_path <-
    system.file(package = "gDRtestData", "testdata", "finalMAE_combo_matrix_small.RDS")
  fmae_cms <- readRDS(fmae_cms_path)
  se1 <- fmae_cms[["matrix"]]
  assays(se1) <- assays(se1)["Averaged"]
  
  new_se1 <- purrr::quietly(gDRcore::fit_SE.combinations)(se1[1, 1])
  exp_as <-
    c(
      "Averaged",
      "SmoothMatrix",
      "BlissExcess",
      "HSAExcess",
      "all_iso_points",
      "isobolograms",
      "BlissScore",
      "HSAScore",
      "CIScore_50",
      "CIScore_80",
      "Metrics"
    )
  expect_equal(assayNames(new_se1$result), exp_as)
 expect_length(new_se1$warnings, 57)
  
  aip_df <-
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(new_se1$result, "all_iso_points"))
  expect_true(all(dim(aip_df) > 0))
})
