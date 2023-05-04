test_that("fit_SE errors as expected", {
  se <- SummarizedExperiment::SummarizedExperiment()
  expect_error(fit_SE(se = 1), 
               "'se' failed: Must inherit from class 'SummarizedExperiment'")
  expect_error(fit_SE(se = se, data_type = 1), 
               "'data_type' failed: Must be of type 'string', not 'double'")
  expect_error(fit_SE(se = se, data_type = "dummy"), 
               "'data_type' failed: Must be element of set {'single-agent','cotreatment','co-dilution'}, ", 
               fixed = TRUE)
  expect_error(fit_SE(se = se, data_type = "single-agent", nested_identifiers = 1), 
               "'nested_identifiers' failed: Must be of type 'character' (or 'NULL'), not 'double'", 
               fixed = TRUE)
  expect_error(fit_SE(se, averaged_assay = 1), 
               "'averaged_assay' failed: Must be of type 'string', not 'double'")
  expect_error(fit_SE(se, averaged_assay = "dummy"), 
               "'dummy' is not on of the available assays")
  expect_error(fit_SE(se, metrics_assay = 1), 
               "'metrics_assay' failed: Must be of type 'string', not 'double'.")
  expect_error(fit_SE(se, curve_type = 1), 
               "'curve_type' failed: Must be of type 'character', not 'double'.")
  expect_error(fit_SE(se, curve_type = c("GR", "dummy")), 
               "'all(curve_type %in% c(\"GR\", \"RV\"))' failed: Must be TRUE.", 
               fixed = TRUE)
  
  maeReal <- readRDS(system.file("testdata", 
                                 "finalMAE_combo_2dose_nonoise2.RDS", package = "gDRtestData")
  )
  
  nrm_se <- normalize_SE(se = maeReal[[2]], data_type = "single-agent")
  avg_se <- average_SE(se = nrm_se, data_type = "single-agent")
  
  expect_error(fit_SE(avg_se, averaged_assay = "dummy"), 
               "'dummy' is not on of the available assays")
  expect_error(fit_SE(avg_se, metrics_assay = "dummy"), 
               "'dummy' is not on of the available assays")
})

test_that("fit_SE works as expected", {
  maeReal <- readRDS(system.file("testdata", 
                                 "finalMAE_combo_2dose_nonoise2.RDS", package = "gDRtestData")
  )
  
  nrm_se <- normalize_SE(se = maeReal[[2]], data_type = "single-agent")
  avg_se <- average_SE(se = nrm_se, data_type = "single-agent")
  fit_se <- fit_SE(avg_se)
  expect_class(fit_se, "SummarizedExperiment")
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
