test_that("fit_SE errors as expected", {
  se <- SummarizedExperiment::SummarizedExperiment()
  expect_error(fit_SE(se = 1), 
               "'se' failed: Must inherit from class 'SummarizedExperiment'")
  expect_error(fit_SE(se = se, data_type = 1), 
               "'data_type' failed: Must be of type 'string', not 'double'")
  expect_error(fit_SE(se = se, data_type = "dummy"), 
               "'data_type' failed: Must be element of set {'single-agent','co-dilution'}, ", 
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
  
  ll <- names(maeReal@ExperimentList)
  se <- maeReal[[ll[which(ll == "single-agent")[1]]]]
  SummarizedExperiment::assay(se, "Metrics") <- NULL
  
  expect_error(fit_SE(se, averaged_assay = "dummy"), 
               "'dummy' is not on of the available assays")
})

test_that("fit_SE works as expected", {
  maeReal <- readRDS(system.file("testdata", 
                                 "finalMAE_combo_2dose_nonoise2.RDS", package = "gDRtestData")
  )
  
  ll <- names(maeReal@ExperimentList)
  se <- maeReal[[ll[which(ll == "single-agent")[1]]]]
  SummarizedExperiment::assay(se, "Metrics") <- NULL
  ext_ass <- SummarizedExperiment::assayNames(se)
  
  fit_se <- fit_SE(se, metrics_assay = "testing")
  expect_class(fit_se, "SummarizedExperiment")
  expect_identical(SummarizedExperiment::assayNames(fit_se), c(ext_ass, "testing"))
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
