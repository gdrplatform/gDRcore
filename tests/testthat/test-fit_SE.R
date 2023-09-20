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
  
  maeReal <- gDRutils::get_synthetic_data("finalMAE_combo_2dose_nonoise2")
  
  se <- MultiAssayExperiment::experiments(maeReal)["single-agent"][[1]]
  SummarizedExperiment::assay(se, "Metrics") <- NULL
  
  expect_error(fit_SE(se, averaged_assay = "dummy"), 
               "'dummy' is not on of the available assays")
})

test_that("fit_SE works as expected", {
  maeReal <- gDRutils::get_synthetic_data("finalMAE_combo_2dose_nonoise2")
  
  se <- MultiAssayExperiment::experiments(maeReal)["single-agent"][[1]]
  SummarizedExperiment::assay(se, "Metrics") <- NULL
  # to avoid warnings about overwriting existing metadata entry
  S4Vectors::metadata(se)[[".internal"]] <- NULL
  S4Vectors::metadata(se)[["fit_parameters"]] <- NULL
  ext_ass <- SummarizedExperiment::assayNames(se)
  
  fit_se <- fit_SE(se, metrics_assay = "testing")
  expect_class(fit_se, "SummarizedExperiment")
  expect_identical(SummarizedExperiment::assayNames(fit_se), c(ext_ass, "testing"))
  expect_true(all(vapply(
    SummarizedExperiment::assayNames(fit_se), 
    function(x) is.character(rownames(SummarizedExperiment::assays(fit_se)[[x]][1])[[1]]),
    logical(1)
  )))
})

test_that("fit_FUN works as expected", {
  elem <- S4Vectors::DataFrame(
    normalization_type = factor(rep(c("RV", "GR"), each = 9)),
    Concentration = rep(c(10 ^ (seq(-3, 1, 0.5))), 2),
    x = c(0.9999964, 0.9999640, 0.9996401, 0.9964143, 0.9653846, 0.7428571, 0.2800000, 
          0.1219512, 0.1022444, 0.9999944, 0.9999440, 0.9994402, 0.9944223, 0.9461538, 
          0.6000000, -0.1200000, -0.3658537, -0.3965087
    ),
    x_std = rep(0.1, 18)
  )
  metric_cols <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
  
  res <- fit_FUN(elem, 
                 nested_identifiers = "Concentration",
                 n_point_cutoff = 4,
                 range_conc = c(5e-3, 5),
                 force_fit = FALSE,
                 pcutoff = 0.05,
                 cap = 0.1,
                 curve_type = c("GR", "RV"))
  expect_true(all(rownames(res) %in% c("GR_gDR", "RV_gDR")))
  expect_true(all(metric_cols %in% colnames(res)))
  
  res_null <- fit_FUN(NULL, 
                      nested_identifiers = "Concentration",
                      n_point_cutoff = 4,
                      range_conc = c(5e-3, 5),
                      force_fit = FALSE,
                      pcutoff = 0.05,
                      cap = 0.1,
                      curve_type = c("GR", "RV"))
  expect_true(all(rownames(res_null) %in% c("GR", "RV")))
  expect_true(all(metric_cols %in% colnames(res_null)))
})

test_that("fit_SE.combinations works as expected", {
  
  # combo data 
  fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
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
  expect_true(all(vapply(
    exp_as, 
    function(x) is.character(rownames(SummarizedExperiment::assays(new_se1$result)[[x]][1])[[1]]),
    logical(1)
  )))
  
  aip_df <-
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(new_se1$result, "all_iso_points"))
  expect_true(all(dim(aip_df) > 0))
})
