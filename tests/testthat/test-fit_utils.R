suppressWarnings(
  se_small <- qs2::qs_read(
    system.file("testdata/finalMAE_small.qs2", package = "gDRtestData")
  )[["single-agent"]]
)
# Helper for combination SE: 2 drugs x 1 cell line, 3x3 conc grid + SA edges
.build_combo_se <- function(drug1 = "DRUG_A", drug2 = "DRUG_B",
                            cl = "CL_1",
                            norm_types = c("GR", "RV"),
                            conc1 = c(0, 0.1, 1.0),
                            conc2 = c(0, 0.1, 1.0),
                            seed = 99L) {
  dt <- data.table::CJ(
    row = drug1,
    column = cl,
    normalization_type = norm_types,
    Concentration = conc1,
    Concentration_2 = conc2
  )
  set.seed(seed)
  dt[, x := runif(.N, 0.1, 1.0)]

  data_cols <- setdiff(names(dt), c("row", "column"))
  bumpy <- BumpyMatrix::splitAsBumpyMatrix(
    dt[, data_cols, with = FALSE],
    row = dt[["row"]], col = dt[["column"]]
  )
  SummarizedExperiment::SummarizedExperiment(assays = list(Averaged = bumpy))
}



####
# bliss_fit_fn
####

test_that("bliss_fit_fn returns expected columns for valid combo input", {
  se <- .build_combo_se()
  dt <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, "Averaged"),
    row.field = "row", column.field = "column"
  )
  dt <- data.table::as.data.table(dt)
  sub_rv <- dt[dt$normalization_type == "RV", ]

  result <- bliss_fit_fn(sub_rv)
  expect_type(result, "list")
  expect_named(result,
    c("normalization_type", "bliss_score", "bliss_excess_mean", "n_combo_points"),
    ignore.order = TRUE
  )
  expect_equal(result$normalization_type, "RV")
  expect_gt(result$n_combo_points, 0L)
})

test_that("bliss_fit_fn returns NA when no combo points present", {
  # Only single-agent rows (Concentration_2 == 0 or Concentration == 0)
  dt <- data.table::data.table(
    normalization_type = "RV",
    Concentration = c(0, 0.1, 1.0),
    Concentration_2 = c(0, 0, 0),
    x = c(1.0, 0.8, 0.5)
  )
  result <- bliss_fit_fn(dt)
  expect_true(is.na(result$bliss_score))
  expect_equal(result$n_combo_points, 0L)
})

test_that("bliss_fit_fn integrates with apply_fit on combination data", {
  se <- .build_combo_se()
  se_out <- apply_fit(
    se, bliss_fit_fn, "combination",
    output_assay = "custom_bliss",
    fit_source = "bliss"
  )
  expect_true("custom_bliss" %in% SummarizedExperiment::assayNames(se_out))
  df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "custom_bliss"),
    row.field = "row", column.field = "column"
  )
  expect_true("bliss_score" %in% names(df))
  # Two norm types → two rows per (drug x cell line)
  expect_equal(NROW(df), 2L)
})


####
# hss_fit_fn
####

test_that("hss_fit_fn returns expected columns for valid combo input", {
  se <- .build_combo_se()
  dt <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, "Averaged"),
    row.field = "row", column.field = "column"
  )
  dt <- data.table::as.data.table(dt)
  sub_rv <- dt[dt$normalization_type == "RV", ]

  result <- hss_fit_fn(sub_rv)
  expect_type(result, "list")
  expect_named(result,
    c("normalization_type", "hss_score", "hss_excess_mean", "n_combo_points"),
    ignore.order = TRUE
  )
  expect_gt(result$n_combo_points, 0L)
})

test_that("hss_fit_fn returns NA when no combo points present", {
  dt <- data.table::data.table(
    normalization_type = "GR",
    Concentration = c(0, 0.1, 1.0),
    Concentration_2 = c(0, 0, 0),
    x = c(1.0, 0.7, 0.4)
  )
  result <- hss_fit_fn(dt)
  expect_true(is.na(result$hss_score))
})

test_that("hss_fit_fn integrates with apply_fit on combination data", {
  se <- .build_combo_se()
  se_out <- apply_fit(
    se, hss_fit_fn, "combination",
    output_assay = "custom_hss",
    fit_source = "hss"
  )
  expect_true("custom_hss" %in% SummarizedExperiment::assayNames(se_out))
  df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "custom_hss"),
    row.field = "row", column.field = "column"
  )
  expect_true("hss_score" %in% names(df))
})

test_that("chaining bliss and hss on the same SE adds both assays independently", {
  se <- .build_combo_se(seed = 200L)
  se_out <- apply_fit(se, bliss_fit_fn, "combination",
                     output_assay = "custom_bliss", fit_source = "bliss")
  se_out <- apply_fit(se_out, hss_fit_fn, "combination",
                     output_assay = "custom_hss", fit_source = "hss")

  expect_true("custom_bliss" %in% SummarizedExperiment::assayNames(se_out))
  expect_true("custom_hss" %in% SummarizedExperiment::assayNames(se_out))
})


test_that("apply_fits on combination data: bliss + hss in one pass", {
  se <- .build_combo_se(seed = 34L)

  se_out <- apply_fits(
    se,
    fit_fns = list(custom_bliss = bliss_fit_fn, custom_hss = hss_fit_fn),
    data_type = "combination",
    fit_source = "synergy"
  )

  expect_true("custom_bliss" %in% SummarizedExperiment::assayNames(se_out))
  expect_true("custom_hss" %in% SummarizedExperiment::assayNames(se_out))

  df_bliss <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "custom_bliss"),
    row.field = "row", column.field = "column"
  )
  df_hss <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "custom_hss"),
    row.field = "row", column.field = "column"
  )
  expect_true("bliss_score" %in% names(df_bliss))
  expect_true("hss_score" %in% names(df_hss))
  # Two normalization types → two rows per (drug x cell line)
  expect_equal(NROW(df_bliss), 2L)
  expect_equal(NROW(df_hss), 2L)
})
