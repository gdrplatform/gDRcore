test_that("fit_SE errors as expected", {
  se <- SummarizedExperiment::SummarizedExperiment()
  expect_error(fit_SE(se))
})

test_that("fit_SE works as expected", {
  maeReal <-
    readRDS(
      system.file("testdata", "finalMAE_combo_2dose_nonoise2.RDS", package = "gDRtestData")
    )
  
  nrm_se <- normalize_SE(se = maeReal[[1]], data_type = 'single-agent')
  avg_se <- average_SE(se = nrm_se, data_type = "single-agent")
  fit_se <- fit_SE(se = avg_se)
  warns <- capture_warnings(fit_SE(fit_se))
  expect_true(all(warns %in% c("overwriting existing metadata entry: 'fit_parameters'", 
                               "overwriting existing metadata entry: '.internal'")))
})

test_that("fit_SE.combinations works as expected", {
 
  # combo data 
  fmae_cms_path <-
    system.file(package = "gDRtestData", "testdata", "finalMAE_combo_matrix_small.RDS")
  fmae_cms <- readRDS(fmae_cms_path)
  se1 <- fmae_cms[["matrix"]]
  SummarizedExperiment::assays(se1) <- SummarizedExperiment::assays(se1)["Averaged"]
  
  new_se1 <- purrr::quietly(fit_SE.combinations)(se1[1, 1])
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
  expect_equal(SummarizedExperiment::assayNames(new_se1$result), exp_as)

  aip_df <-
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(new_se1$result, "all_iso_points"))
  expect_true(all(dim(aip_df) > 0))
})
