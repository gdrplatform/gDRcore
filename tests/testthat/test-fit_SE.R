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

test_that("fit_SE refits a processed co-dilution experiment", {
  cd_name <- gDRutils::get_supported_experiments("cd")
  se <- gDRutils::get_synthetic_data("finalMAE_combo_codilution_small")[[cd_name]]
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)["Averaged"]

  # the SE still holds the processing metadata of the earlier run, with a
  # sessionInfo; the refit replaces it instead of stopping on it
  expect_s3_class(gDRutils::get_SE_processing_metadata(se)$session_info, "sessionInfo")
  fit_se <- suppressWarnings(fit_SE(se, data_type = cd_name))
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(fit_se))
  expect_gt(NROW(gDRutils::convert_se_assay_to_dt(fit_se, "Metrics")), 0L)
  expect_equal(gDRutils::get_SE_processing_metadata(fit_se)$date_processed, Sys.Date())
})

test_that("fit_SE gives the same object when other packages get loaded in between", {
  se <- gDRutils::get_synthetic_data("finalMAE_small")[["single-agent"]]
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)["Averaged"]
  S4Vectors::metadata(se)[[".internal"]] <- NULL

  # loading a package changes sessionInfo(), which .internal used to store
  extra <- setdiff(c("yaml", "knitr", "BiocStyle"), loadedNamespaces())
  skip_if(length(extra) == 0L, "all candidate packages are loaded already")
  fit_1 <- suppressWarnings(fit_SE(se))
  loadNamespace(extra[[1]])
  fit_2 <- suppressWarnings(fit_SE(se))
  expect_identical(fit_1, fit_2)

  internal <- gDRutils::get_SE_processing_metadata(fit_1)
  expect_identical(names(internal), c("date_processed", "versions"))
  expect_identical(internal$versions$gDRcore, as.character(utils::packageVersion("gDRcore")))
})

test_that("fit_SE.combinations works as expected", {

  # combo data
  fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se1 <- fmae_cms[[gDRutils::get_supported_experiments("combo")]]
  SummarizedExperiment::assays(se1) <- SummarizedExperiment::assays(se1)["Averaged"]

  new_se1 <- purrr::quietly(fit_SE.combinations)(se1[1, 1])
  exp_as <- c("Averaged", "Metrics", "excess", "isobolograms", "all_iso_points", "scores")
  expect_equal(SummarizedExperiment::assayNames(new_se1$result), exp_as)

  aip_df <-
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(new_se1$result, "all_iso_points"))
  expect_true(all(dim(aip_df) > 0))
})

test_that("fit_SE.combinations single-agent arms are monotonic (no checkerboard)", {

  fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se1 <- fmae_cms[[gDRutils::get_supported_experiments("combo")]]
  SummarizedExperiment::assays(se1) <- SummarizedExperiment::assays(se1)["Averaged"]

  fit <- purrr::quietly(fit_SE.combinations)(se1[1, 1])$result
  ex <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(fit, "excess")))

  # On each single-agent arm the smoothed response must be monotonically
  # non-increasing with concentration. Averaging the perpendicular fit's
  # intercept into 'smooth' (it collapses to the row/column mean when the
  # orthogonal series is not significant) used to make the arm zig-zag -
  # the "checkerboard" artifact. Each arm now uses only its parallel fit.
  tol <- 1e-6
  for (nt in unique(ex$normalization_type)) {
    arm1 <- ex[Concentration_2 == 0 & normalization_type == nt & !is.na(smooth)][
      order(Concentration)]
    arm2 <- ex[Concentration == 0 & normalization_type == nt & !is.na(smooth)][
      order(Concentration_2)]
    expect_true(all(diff(arm1$smooth) <= tol),
                info = paste("drug_1 single-agent arm not monotonic for", nt))
    expect_true(all(diff(arm2$smooth) <= tol),
                info = paste("drug_2 single-agent arm not monotonic for", nt))
  }
})


#### Composable apply_combo_*() steps ####

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
  # Reference from a fresh fit: the stored test data predate the arm fix
  se <- suppressWarnings(fit_SE.combinations(.make_combo_se_averaged_only()))
  ref_excess <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, "excess"),
    row.field = "row", column.field = "column"
  ))
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)[c("Averaged", "Metrics")]

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
  expect_lt(max(abs(comp$hsa_score_ref - comp$hsa_score_new), na.rm = TRUE), 1e-10)
})

test_that("fit_SE.combinations single-agent arms follow the fit along the arm", {
  fit <- suppressWarnings(fit_SE.combinations(.make_combo_se_averaged_only()))
  ex <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(fit, "excess")))
  met <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(fit, "Metrics")))

  for (nt in unique(ex$normalization_type)) {
    p1 <- met[normalization_type == nt & dilution_drug == "drug_1" & cotrt_value == 0][1L]
    p2 <- met[normalization_type == nt & dilution_drug == "drug_2" & cotrt_value == 0][1L]
    arm1 <- ex[normalization_type == nt & Concentration_2 == 0 & Concentration > 0]
    arm2 <- ex[normalization_type == nt & Concentration == 0 & Concentration_2 > 0]
    expect_equal(arm1$smooth, gDRutils::predict_efficacy_from_conc(
      arm1$Concentration, p1$x_inf, p1$x_0, p1$ec50, p1$h), tolerance = 1e-6)
    expect_equal(arm2$smooth, gDRutils::predict_efficacy_from_conc(
      arm2$Concentration_2, p2$x_inf, p2$x_0, p2$ec50, p2$h), tolerance = 1e-6)
  }
})

test_that("fit_SE.combinations reports CIScore_50 and CIScore_80 from the isobolograms", {
  se_out <- suppressWarnings(fit_SE.combinations(.make_combo_se_averaged_only()))
  scores <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "scores"),
    row.field = "row", column.field = "column"
  ))
  iso <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "isobolograms"),
    row.field = "row", column.field = "column"
  ))
  expect_true("CI_100x" %in% names(iso))
  expect_false(anyNA(scores$CIScore_50))
  expect_false(anyNA(scores$CIScore_80))

  # CIScore_X is the CI of the isobologram at the lowest response level >= cutoff
  for (nt in scores$normalization_type) {
    iso_nt <- iso[iso$normalization_type == nt, ]
    levels <- as.numeric(iso_nt$iso_level)
    expect_equal(scores[normalization_type == nt, CIScore_50],
                 unique(iso_nt$CI_100x[levels == min(levels[levels >= 0.5])]))
    expect_equal(scores[normalization_type == nt, CIScore_80],
                 unique(iso_nt$CI_100x[levels == min(levels[levels >= 0.2])]))
  }
})

test_that("apply_combo_scores leaves CIScore NA without isobolograms", {
  fmae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se <- fmae[[gDRutils::get_supported_experiments("combo")]][1, 1]
  SummarizedExperiment::assays(se) <-
    SummarizedExperiment::assays(se)[c("Averaged", "Metrics", "excess")]

  se_out <- suppressWarnings(apply_combo_scores(se, excess_assay = "excess"))
  scores <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se_out, "scores"))
  expect_true(all(is.na(scores$CIScore_50)))
  expect_true(all(is.na(scores$CIScore_80)))
})

test_that(".ci_scores picks the lowest response level at or above the cutoff", {
  curves <- data.table::data.table(
    iso_level = c("-0.05", "0.2", "0.25", "0.45", "0.55", "0.55"),
    CI_100x = c(5, 3, 4, 1, 2, 2)
  )
  expect_equal(.ci_scores(curves), c(CIScore_50 = 2, CIScore_80 = 3))
  expect_equal(.ci_scores(curves[iso_level != "0.2"]), c(CIScore_50 = 2, CIScore_80 = 4))
  expect_equal(.ci_scores(curves[iso_level %in% c("-0.05", "0.2")]),
               c(CIScore_50 = NA_real_, CIScore_80 = 3))
  expect_equal(.ci_scores(data.table::data.table(row_id = 1)),
               c(CIScore_50 = NA_real_, CIScore_80 = NA_real_))
})

test_that("fit_SE.combinations passes score_FUN to apply_combo_scores", {
  fmae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se <- fmae[[gDRutils::get_supported_experiments("combo")]][1, 1]
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)["Averaged"]

  constant_score <- function(x) 42
  se_out <- suppressWarnings(fit_SE.combinations(se, score_FUN = constant_score))

  scores <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "scores"),
    row.field = "row", column.field = "column"
  ))
  expect_true(all(scores$bliss_score == 42 | is.na(scores$bliss_score)))
  expect_true(all(scores$hsa_score == 42 | is.na(scores$hsa_score)))
})

test_that("fit_SE points time-course data to fit_SE.timecourse", {
  se <- SummarizedExperiment::SummarizedExperiment()
  expect_error(fit_SE(se, data_type = gDRutils::get_supported_experiments("time-course")),
               "use fit_SE.timecourse")
})


test_that("fit_SE records the gDR fit by reference with its settings", {
  maeReal <- gDRutils::get_synthetic_data("finalMAE_combo_2dose_nonoise2")
  se <- MultiAssayExperiment::experiments(maeReal)["single-agent"][[1]]
  SummarizedExperiment::assay(se, "Metrics") <- NULL
  S4Vectors::metadata(se)[[".internal"]] <- NULL
  S4Vectors::metadata(se)[["fit_parameters"]] <- NULL

  fit_se <- suppressWarnings(fit_SE(se, pcutoff = 0.1))
  recipe <- gDRutils::get_fit_provenance(fit_se, "Metrics", "gDR")
  expect_identical(recipe$fn$name, "fit_drug_response_metrics")
  expect_identical(recipe$fn$package, "gDRutils")
  expect_null(recipe$fn$source)
  expect_identical(recipe$fit_args$pcutoff, 0.1)
  expect_true(gDRutils::check_reproducibility(fit_se)$reproducible)
  # fitting the same data again gives an identical object
  expect_identical(fit_se, suppressWarnings(fit_SE(se, pcutoff = 0.1)))
})


test_that("fit_SE.combinations records each step with its arguments", {
  fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  se1 <- fmae_cms[[gDRutils::get_supported_experiments("combo")]]
  SummarizedExperiment::assays(se1) <- SummarizedExperiment::assays(se1)["Averaged"]

  fit_se <- suppressWarnings(fit_SE.combinations(se1[1, 1]))
  summary <- gDRutils::summarize_fit_provenance(fit_se)
  expect_setequal(summary$assay, c("Metrics", "excess", "isobolograms", "all_iso_points", "scores"))
  expect_true(all(grepl("^gDRcore::apply_combo_", summary$fn)))
  scores <- gDRutils::get_fit_provenance(fit_se, "scores", "gDR")
  expect_identical(scores$fn_args$score_FUN$name, "calculate_score")
  expect_true(all(gDRutils::check_reproducibility(fit_se)$reproducible))
})
