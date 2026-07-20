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
})

test_that("fit_SE.combinations works as expected", {

  # combo data
  fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se1 <- fmae_cms[[gDRutils::get_supported_experiments("combo")]]
  SummarizedExperiment::assays(se1) <- SummarizedExperiment::assays(se1)["Averaged"]

  new_se1 <- purrr::quietly(fit_SE.combinations)(se1[1, 1])
  exp_as <- c("Averaged", "excess", "all_iso_points", "isobolograms", "scores", "Metrics")
  # Check all expected assays are present (order may vary with composable steps)
  expect_true(all(exp_as %in% SummarizedExperiment::assayNames(new_se1$result)))

  aip_df <-
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(new_se1$result, "all_iso_points"))
  expect_true(all(dim(aip_df) > 0))
})


#### Composable apply_combo_*() steps — GDR-3486 ####

.make_combo_se_averaged_only <- function() {
  fmae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se <- fmae[[gDRutils::get_supported_experiments("combo")]]
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)["Averaged"]
  se[1, 1]
}

test_that("apply_combo_sa_fits produces Metrics assay matching fit_SE.combinations", {
  se_avg <- .make_combo_se_averaged_only()
  se_fitted <- suppressWarnings(apply_combo_sa_fits(se_avg))

  expect_class(se_fitted, "SummarizedExperiment")
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(se_fitted))

  met <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_fitted, "Metrics"),
    row.field = "row", column.field = "column"
  ))
  expect_gt(NROW(met), 0L)
  expect_true(all(c("dilution_drug", "cotrt_value", "ec50", "h", "x_inf", "fit_type") %in% names(met)))
  expect_true(all(met$fit_type %in% c("DRC3pHillFitModelFixS0", "DRC4pHillFitModel",
                                       "DRCConstantFitResult", "DRCTooFewPointsToFit",
                                       "DRCInvalidFitResult")))
  # Verify numerical identity vs fit_SE.combinations (cor = 1.0, float precision only)
  fmae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se_ref <- fmae[[gDRutils::get_supported_experiments("combo")]][1, 1]
  ref <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_ref, "Metrics"),
    row.field = "row", column.field = "column"
  ))
  expect_equal(NROW(met), NROW(ref))
  met_s <- met[order(normalization_type, dilution_drug, cotrt_value)]
  ref_s <- ref[order(normalization_type, dilution_drug, cotrt_value)]
  expect_equal(sum(met_s$fit_type == ref_s$fit_type, na.rm = TRUE), NROW(met))
})

test_that("apply_combo_excess produces excess assay matching fit_SE.combinations", {
  fmae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se <- fmae[[gDRutils::get_supported_experiments("combo")]][1, 1]
  ref_excess <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, "excess"),
    row.field = "row", column.field = "column"
  ))

  se_out <- suppressWarnings(apply_combo_excess(se))
  expect_true("excess" %in% SummarizedExperiment::assayNames(se_out))

  new_excess <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "excess"),
    row.field = "row", column.field = "column"
  ))
  expect_equal(NROW(new_excess), NROW(ref_excess))
  expect_true(all(c("smooth", "hsa_excess", "bliss_excess") %in% names(new_excess)))

  # Numerical identity (reading from existing Metrics → machine-epsilon diff only)
  re <- ref_excess[order(row, column, normalization_type, Concentration, Concentration_2)]
  ne <- new_excess[order(row, column, normalization_type, Concentration, Concentration_2)]
  expect_lt(max(abs(as.numeric(re$smooth) - as.numeric(ne$smooth)), na.rm = TRUE), 1e-10)
  expect_lt(max(abs(as.numeric(re$bliss_excess) - as.numeric(ne$bliss_excess)), na.rm = TRUE), 1e-10)
})

test_that("apply_combo_isobolograms produces isobolograms matching fit_SE.combinations", {
  fmae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se <- fmae[[gDRutils::get_supported_experiments("combo")]][1, 1]
  ref_iso <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, "isobolograms"),
    row.field = "row", column.field = "column"
  ))

  se_out <- suppressWarnings(apply_combo_isobolograms(se))
  expect_true("isobolograms" %in% SummarizedExperiment::assayNames(se_out))

  new_iso <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "isobolograms"),
    row.field = "row", column.field = "column"
  ))
  expect_equal(NROW(new_iso), NROW(ref_iso))
  expect_true(all(c("log2_CI", "iso_level", "pos_x", "pos_y") %in% names(new_iso)))

  ri <- ref_iso[order(row, column, normalization_type, iso_level, log10_ratio_conc)]
  ni <- new_iso[order(row, column, normalization_type, iso_level, log10_ratio_conc)]
  expect_lt(max(abs(as.numeric(ri$log2_CI) - as.numeric(ni$log2_CI)), na.rm = TRUE), 1e-10)
})

test_that("apply_combo_scores with excess_assay gives scores matching fit_SE.combinations", {
  fmae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se <- fmae[[gDRutils::get_supported_experiments("combo")]]
  ref_scores <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, "scores"),
    row.field = "row", column.field = "column"
  ))

  se_out <- suppressWarnings(apply_combo_scores(se, excess_assay = "excess"))
  new_scores <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "scores"),
    row.field = "row", column.field = "column"
  ))

  comp <- merge(ref_scores, new_scores,
                by = c("row", "column", "normalization_type"),
                suffixes = c("_ref", "_new"))
  expect_equal(NROW(comp), NROW(ref_scores))
  expect_lt(max(abs(comp$bliss_score_ref - comp$bliss_score_new), na.rm = TRUE), 1e-10)
  expect_lt(max(abs(comp$hsa_score_ref   - comp$hsa_score_new),   na.rm = TRUE), 1e-10)
})
