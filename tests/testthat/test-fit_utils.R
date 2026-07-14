suppressWarnings(
  se_small <- qs2::qs_read(
    system.file("testdata/finalMAE_small.qs2", package = "gDRtestData")
  )[["single-agent"]]
)

# Simple fit_fn that returns a fixed named list
simple_fit_fn <- function(dt) {
  list(x_mean = mean(dt[["x"]], na.rm = TRUE), ec50 = 1e-6)
}

# Helper to build a minimal SE for testing — avoids repeating CJ + split + SE pattern
.build_test_se <- function(drug_ids = "DRUG_A",
                           cl_ids = "CL_1",
                           norm_types = c("GR", "RV"),
                           conc_vals = c(0.01, 0.1, 1.0),
                           x_col = "x",
                           seed = 42L,
                           x_min = 0,
                           x_max = 1) {
  avg_rows <- data.table::CJ(
    row = drug_ids,
    column = cl_ids,
    normalization_type = norm_types,
    Concentration = conc_vals
  )
  set.seed(seed)
  avg_rows[, (x_col) := runif(.N, min = x_min, max = x_max)]

  data_cols <- setdiff(names(avg_rows), c("row", "column"))
  bumpy <- BumpyMatrix::splitAsBumpyMatrix(
    avg_rows[, data_cols, with = FALSE],
    row = avg_rows[["row"]],
    col = avg_rows[["column"]]
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Averaged = bumpy)
  )
  list(se = se, avg_rows = avg_rows)
}


test_that("apply_fit_to_se validates inputs correctly", {
  expect_error(
    apply_fit_to_se("not_an_se", simple_fit_fn),
    regexp = "se"
  )
  expect_error(
    apply_fit_to_se(se_small, "not_a_function"),
    regexp = "fit_fn"
  )
  expect_error(
    apply_fit_to_se(se_small, simple_fit_fn, normalization_types = character(0)),
    regexp = "normalization_types|slicing_values"
  )
  expect_error(
    apply_fit_to_se(se_small, simple_fit_fn, averaged_assay = "BOGUS"),
    regexp = "assay is required"
  )
  expect_error(
    apply_fit_to_se(se_small, simple_fit_fn, merge = "bad_value"),
    regexp = "merge"
  )
  expect_error(
    apply_fit_to_se(se_small, simple_fit_fn, on_error = "bad_value"),
    regexp = "on_error"
  )
})


test_that("apply_fit_to_se iterates over all (drug x cell line x normalization_type) triplets", {
  n_drugs <- nrow(se_small)
  n_cls <- ncol(se_small)
  norm_types <- c("GR", "RV")

  call_log <- list()
  log_fn <- function(dt) {
    call_log[[length(call_log) + 1L]] <<- unique(dt[["normalization_type"]])
    list(x_mean = 0.5, ec50 = 1e-6)
  }

  suppressWarnings(result_se <- apply_fit_to_se(se_small, log_fn,
    normalization_types = norm_types, fit_source = "test"))

  # Called once per triplet (skipping empties)
  expect_lte(length(call_log), n_drugs * n_cls * length(norm_types))
  # Each call received only one normalization_type value
  lapply(call_log, function(nt) expect_length(nt, 1L))

  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row", column.field = "column"
  )
  new_rows <- metrics_df[metrics_df[["fit_source"]] == "test", ]
  # One row per norm_type per non-empty triplet
  expect_equal(
    sort(unique(new_rows[["normalization_type"]])),
    sort(norm_types)
  )
})


test_that("apply_fit_to_se filters by normalization_type before calling fit_fn", {
  # fit_fn receives only rows matching the requested normalization_type
  seen_types <- character(0)
  check_fn <- function(dt) {
    seen_types <<- c(seen_types, as.character(unique(dt[["normalization_type"]])))
    list(x_mean = 0.5)
  }

  suppressWarnings(apply_fit_to_se(se_small, check_fn,
    normalization_types = "GR", fit_source = "chk"))

  expect_true(all(seen_types == "GR"))
})


test_that("apply_fit_to_se calls fit_fn exactly once per triplet on a fully-populated SE", {
  drug_ids <- c("DRUG_A", "DRUG_B")
  cl_ids <- c("CL_1", "CL_2")
  norm_types <- c("GR", "RV")
  built <- .build_test_se(drug_ids, cl_ids, norm_types, c(0.1, 1.0, 10.0), seed = 7L)
  ctrl_se <- built$se

  call_count <- 0L
  count_fn <- function(dt) {
    call_count <<- call_count + 1L
    list(metric1 = 1.0)
  }

  n_expected <- length(drug_ids) * length(cl_ids) * length(norm_types)
  suppressWarnings(apply_fit_to_se(ctrl_se, count_fn,
    normalization_types = norm_types, fit_source = "exact_count_test"))

  expect_equal(call_count, n_expected)
})


test_that("apply_fit_to_se passes Averaged assay columns unmodified to fit_fn", {
  # Build a minimal SE with a custom data column to verify fit_fn receives
  # exactly those columns (no extras, no omissions).
  norm_types <- c("GR", "RV")
  built <- .build_test_se(
    drug_ids = "DRUG_X", cl_ids = "CL_X",
    norm_types = norm_types, conc_vals = c(0.01, 0.1, 1.0),
    x_col = "viability"
  )
  col_se <- built$se
  data_cols <- setdiff(names(built$avg_rows), c("row", "column"))

  captured_cols <- list()
  col_fn <- function(dt) {
    captured_cols[[length(captured_cols) + 1L]] <<- sort(names(dt))
    list(metric_a = 1.0)
  }

  suppressWarnings(apply_fit_to_se(col_se, col_fn,
    normalization_types = norm_types, fit_source = "col_check_test"))

  expected_cols <- sort(data_cols)
  expect_true(length(captured_cols) > 0L)
  lapply(captured_cols, function(cols) {
    expect_equal(cols, expected_cols)
  })
})


test_that("apply_fit_to_se records fit_source and normalization_type in output", {
  suppressWarnings(result_se <- apply_fit_to_se(se_small, simple_fit_fn,
    normalization_types = c("GR", "RV"), fit_source = "my_fit"))

  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row", column.field = "column"
  )
  new_rows <- metrics_df[
    !is.na(metrics_df[["fit_source"]]) & metrics_df[["fit_source"]] == "my_fit",
  ]
  expect_true(nrow(new_rows) > 0L)
  expect_true(all(new_rows[["fit_source"]] == "my_fit"))
  expect_true(all(new_rows[["normalization_type"]] %in% c("GR", "RV")))
})


test_that("apply_fit_to_se on_error='warn' emits warning and skips failed triplet", {
  err_fn <- function(dt) stop("deliberate test error")

  warnings_caught <- character(0)
  result_se <- withCallingHandlers(
    apply_fit_to_se(se_small, err_fn, on_error = "warn"),
    warning = function(w) {
      warnings_caught <<- c(warnings_caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Warnings emitted for each failing triplet
  expect_true(length(warnings_caught) > 0L)
  expect_true(any(grepl("deliberate test error", warnings_caught)))

  # No new custom rows written (all triplets failed)
  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_false("custom" %in% metrics_df[["fit_source"]])
})


test_that("apply_fit_to_se on_error='stop' propagates error immediately", {
  err_fn <- function(dt) stop("deliberate stop error")

  expect_error(
    apply_fit_to_se(se_small, err_fn, on_error = "stop"),
    regexp = "deliberate stop error"
  )
})


test_that("apply_fit_to_se merge='merge' is idempotent for same fit_source+normalization_type", {
  suppressWarnings(result_se1 <- apply_fit_to_se(se_small, simple_fit_fn,
    fit_source = "idempotent_test"))
  suppressWarnings(result_se2 <- apply_fit_to_se(result_se1, simple_fit_fn,
    fit_source = "idempotent_test"))

  metrics1 <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se1, "Metrics"),
    row.field = "row", column.field = "column"
  )
  metrics2 <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se2, "Metrics"),
    row.field = "row", column.field = "column"
  )

  expect_equal(nrow(metrics1), nrow(metrics2))
  test_rows1 <- metrics1[
    !is.na(metrics1[["fit_source"]]) &
      metrics1[["fit_source"]] == "idempotent_test",
  ]
  test_rows2 <- metrics2[
    !is.na(metrics2[["fit_source"]]) &
      metrics2[["fit_source"]] == "idempotent_test",
  ]
  expect_equal(nrow(test_rows1), nrow(test_rows2))
})


test_that("apply_fit_to_se merge='merge' preserves rows with a different fit_source", {
  existing_fit_source <- "gDR"

  suppressWarnings(result_se <- apply_fit_to_se(se_small, simple_fit_fn,
    fit_source = "new_fit"))

  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_true(existing_fit_source %in% metrics_df[["fit_source"]])
  expect_true("new_fit" %in% metrics_df[["fit_source"]])
})


test_that("apply_fit_to_se merge='replace' overwrites the entire Metrics assay", {
  suppressWarnings(result_se <- apply_fit_to_se(se_small, simple_fit_fn,
    merge = "replace", fit_source = "replace_fit"))

  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row", column.field = "column"
  )
  # Only the new fit_source should appear — old "gDR" rows gone
  fit_sources <- unique(metrics_df[["fit_source"]])
  expect_false("gDR" %in% fit_sources)
  expect_true("replace_fit" %in% fit_sources)
})


test_that("apply_fit_to_se warns for missing response_metrics columns in fit_fn output", {
  # fit_fn returns a list with none of the standard response_metrics columns
  incomplete_fn <- function(dt) list(my_custom_metric = 42)

  expect_warning(
    apply_fit_to_se(se_small, incomplete_fn, fit_source = "incomplete"),
    regexp = "missing"
  )
})


test_that("apply_fit_to_se returns unchanged se when no triplets produce results", {
  # Use a normalization_type that does not exist in the data so all sub_dt
  # are empty and the function returns the SE unchanged.
  result_se <- apply_fit_to_se(se_small, simple_fit_fn,
    normalization_types = "NONEXISTENT", fit_source = "empty_test")

  # Metrics assay should be identical to the original (no new rows written)
  orig_metrics <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_small, "Metrics"),
    row.field = "row", column.field = "column"
  )
  result_metrics <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_equal(nrow(orig_metrics), nrow(result_metrics))
  expect_false("empty_test" %in% result_metrics[["fit_source"]])
})


test_that("apply_fit_to_se creates Metrics assay when not present", {
  no_metrics_se <- se_small
  SummarizedExperiment::assay(no_metrics_se, "Metrics") <- NULL

  suppressWarnings(result_se <- apply_fit_to_se(no_metrics_se, simple_fit_fn,
    fit_source = "new_metrics"))

  expect_true("Metrics" %in% SummarizedExperiment::assayNames(result_se))
  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_true("new_metrics" %in% metrics_df[["fit_source"]])
})


test_that("apply_fit_to_se end-to-end integration: minimal SE, trivial fit_fn, Metrics assay verified", {
  drug_ids <- c("DRUG_A", "DRUG_B")
  cl_ids <- c("CL_1", "CL_2")
  norm_types <- c("GR", "RV")
  conc_vals <- c(0.1, 0.3, 1.0, 3.0, 10.0)

  built <- .build_test_se(drug_ids, cl_ids, norm_types, conc_vals,
    seed = 42L, x_min = 0.3, x_max = 1.0)
  minimal_se <- built$se
  avg_rows <- built$avg_rows

  trivial_fit_fn <- function(dt) {
    list(x_mean = mean(dt[["x"]], na.rm = TRUE), ec50 = 1e-6)
  }

  result_se <- suppressWarnings(
    apply_fit_to_se(
      minimal_se,
      trivial_fit_fn,
      normalization_types = norm_types,
      fit_source = "integration_test"
    )
  )

  expect_true(is(result_se, "SummarizedExperiment"))
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(result_se))

  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field = "row",
    column.field = "column"
  )

  our_rows <- metrics_df[
    !is.na(metrics_df[["fit_source"]]) &
      metrics_df[["fit_source"]] == "integration_test",
  ]

  # One row per (drug × cell line × normalization_type) triplet: 2×2×2 = 8
  n_expected <- length(drug_ids) * length(cl_ids) * length(norm_types)
  expect_equal(nrow(our_rows), n_expected)

  # fit_source correctly stamped on every row.
  expect_true(all(our_rows[["fit_source"]] == "integration_test"))

  # Both normalization_types present in output.
  expect_equal(
    sort(unique(as.character(our_rows[["normalization_type"]]))),
    sort(norm_types)
  )

  # fit_fn output columns are present.
  expect_true("x_mean" %in% names(our_rows))
  expect_true("ec50"   %in% names(our_rows))

  # ec50 is always 1e-6; x_mean is non-NA.
  expect_true(all(!is.na(our_rows[["x_mean"]])))
  expect_true(all(as.numeric(our_rows[["ec50"]]) == 1e-6))

  # Spot-check: x_mean for DRUG_A × CL_1 × GR matches the hand-computed mean.
  sub_dt <- avg_rows[
    row == "DRUG_A" & column == "CL_1" & normalization_type == "GR"
  ]
  expected_mean_val <- mean(sub_dt[["x"]])

  obs_mean_val <- as.numeric(
    our_rows[
      our_rows[["row"]] == "DRUG_A" &
        our_rows[["column"]] == "CL_1" &
        our_rows[["normalization_type"]] == "GR",
    ][["x_mean"]]
  )
  expect_equal(obs_mean_val, expected_mean_val, tolerance = 1e-10)
})


####
# fit_drug_response_metrics tests
####

test_that("fit_drug_response_metrics returns named list with expected fields", {
  dt <- data.table::data.table(
    Concentration = c(0.001, 0.01, 0.1, 1, 10),
    x = c(0.95, 0.8, 0.5, 0.2, 0.1),
    normalization_type = "RV"
  )
  result <- fit_drug_response_metrics(dt)

  expect_true(is.list(result))
  expected_names <- c("fit_source", "normalization_type", "x_mean", "x_AOC",
                      "N_conc", "maxlog10Concentration", "xc50", "h", "r2",
                      "x_0", "x_inf", "fit_type")
  expect_true(all(expected_names %in% names(result)))
  expect_equal(result$normalization_type, "RV")
  expect_equal(result$fit_source, "custom")
  expect_equal(result$N_conc, 5L)
})


test_that("fit_drug_response_metrics returns NA metrics for all-NA input", {
  dt <- data.table::data.table(
    Concentration = c(0.1, 1, 10),
    x = c(NA_real_, NA_real_, NA_real_),
    normalization_type = "GR"
  )
  result <- fit_drug_response_metrics(dt)

  expect_equal(result$fit_type, "DRCInvalidFitResult")
  expect_true(is.na(result$xc50))
  expect_true(is.na(result$h))
  expect_equal(result$N_conc, 0L)
})


test_that("fit_drug_response_metrics uses GR priors for GR normalization", {
  dt <- data.table::data.table(
    Concentration = c(0.001, 0.01, 0.1, 1, 10),
    x = c(0.95, 0.8, 0.5, 0.2, 0.1),
    normalization_type = "GR"
  )
  result <- fit_drug_response_metrics(dt)

  expect_equal(result$normalization_type, "GR")
  expect_true(result$fit_type %in% c("DRC4pHillFitModel", "DRCInvalidFitResult"))
})


test_that("fit_drug_response_metrics computes correct x_mean and x_AOC", {
  x_vals <- c(0.9, 0.7, 0.5, 0.3, 0.1)
  dt <- data.table::data.table(
    Concentration = c(0.001, 0.01, 0.1, 1, 10),
    x = x_vals,
    normalization_type = "RV"
  )
  result <- fit_drug_response_metrics(dt)

  expect_equal(result$x_mean, mean(x_vals))
  expect_equal(result$x_AOC, 1 - mean(x_vals))
})


test_that("fit_drug_response_metrics estimates xc50 fallback when fit fails", {
  # All x > 0.5 — xc50 should be Inf (drug has no effect)
  dt_high <- data.table::data.table(
    Concentration = c(0.1, 1),
    x = c(0.9, 0.8),
    normalization_type = "RV"
  )
  result_high <- fit_drug_response_metrics(dt_high)
  if (result_high$fit_type == "DRCInvalidFitResult") {
    expect_equal(result_high$xc50, Inf)
  }

  # All x <= 0.5 — xc50 should be -Inf (drug very effective)
  dt_low <- data.table::data.table(
    Concentration = c(0.1, 1),
    x = c(0.3, 0.1),
    normalization_type = "RV"
  )
  result_low <- fit_drug_response_metrics(dt_low)
  if (result_low$fit_type == "DRCInvalidFitResult") {
    expect_equal(result_low$xc50, -Inf)
  }
})


####
# .persist_metrics tests
####

test_that(".persist_metrics creates Metrics assay from scratch", {
  built <- .build_test_se(seed = 1L)
  se <- built$se

  new_metrics <- data.table::data.table(
    row = "DRUG_A", column = "CL_1",
    fit_source = "test", normalization_type = "GR",
    x_mean = 0.5, ec50 = 1e-6
  )

  result <- gDRcore:::.persist_metrics(
    se, new_metrics, "replace", "Metrics", "row", "column"
  )
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(result))
  df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_equal(nrow(df), 1L)
  expect_equal(df[["fit_source"]], "test")
})


test_that(".persist_metrics merge mode upserts by fit_source + normalization_type", {
  built <- .build_test_se(seed = 2L)
  se <- built$se

  first <- data.table::data.table(
    row = "DRUG_A", column = "CL_1",
    fit_source = "gDR", normalization_type = "GR",
    x_mean = 0.5
  )
  se <- gDRcore:::.persist_metrics(se, first, "replace", "Metrics", "row", "column")

  second <- data.table::data.table(
    row = "DRUG_A", column = "CL_1",
    fit_source = "custom", normalization_type = "GR",
    x_mean = 0.8
  )
  result <- gDRcore:::.persist_metrics(se, second, "merge", "Metrics", "row", "column")

  df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_equal(nrow(df), 2L)
  expect_true("gDR" %in% df[["fit_source"]])
  expect_true("custom" %in% df[["fit_source"]])
})


test_that(".persist_metrics merge mode replaces rows with matching fit_source", {
  built <- .build_test_se(seed = 3L)
  se <- built$se

  first <- data.table::data.table(
    row = "DRUG_A", column = "CL_1",
    fit_source = "custom", normalization_type = "GR",
    x_mean = 0.5
  )
  se <- gDRcore:::.persist_metrics(se, first, "replace", "Metrics", "row", "column")

  updated <- data.table::data.table(
    row = "DRUG_A", column = "CL_1",
    fit_source = "custom", normalization_type = "GR",
    x_mean = 0.9
  )
  result <- gDRcore:::.persist_metrics(se, updated, "merge", "Metrics", "row", "column")

  df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_equal(nrow(df), 1L)
  expect_equal(as.numeric(df[["x_mean"]]), 0.9)
})

####
# apply_custom_fit
####

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


test_that("apply_custom_fit validates output_assay is provided", {
  built <- .build_test_se()
  se <- built$se
  expect_error(
    apply_custom_fit(se, simple_fit_fn, "single-agent", fit_source = "test"),
    regexp = "output_assay"
  )
})

test_that("apply_custom_fit requires summary_assay when summary_fn is provided", {
  built <- .build_test_se()
  se <- built$se
  expect_error(
    apply_custom_fit(se, simple_fit_fn, "single-agent",
                     output_assay = "out",
                     summary_fn = function(dt) list(n = nrow(dt)),
                     fit_source = "test"),
    regexp = "summary_assay"
  )
})

test_that("apply_custom_fit single-agent produces same result as apply_fit_to_se", {
  built <- .build_test_se(seed = 11L)
  se <- built$se

  via_old <- suppressWarnings(
    apply_fit_to_se(se, simple_fit_fn, metrics_assay = "Metrics",
                    fit_source = "test")
  )
  via_new <- apply_custom_fit(
    se, simple_fit_fn, "single-agent",
    output_assay = "Metrics",
    fit_source = "test"
  )

  df_old <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(via_old, "Metrics"),
    row.field = "row", column.field = "column"
  )
  df_new <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(via_new, "Metrics"),
    row.field = "row", column.field = "column"
  )
  expect_equal(sort(names(df_old)), sort(names(df_new)))
  expect_equal(nrow(df_old), nrow(df_new))
})

test_that("apply_custom_fit writes to a custom-named assay", {
  built <- .build_test_se(seed = 12L)
  se <- built$se

  se_out <- apply_custom_fit(
    se, simple_fit_fn, "single-agent",
    output_assay = "my_custom_metrics",
    fit_source = "test"
  )
  expect_true("my_custom_metrics" %in% SummarizedExperiment::assayNames(se_out))
  expect_false("Metrics" %in% SummarizedExperiment::assayNames(se_out))
})

test_that("apply_custom_fit chaining adds two independent assays", {
  built <- .build_test_se(drug_ids = c("D1", "D2"), cl_ids = c("CL1", "CL2"),
                          seed = 13L)
  se <- built$se

  fn_a <- function(dt) list(metric_a = mean(dt$x, na.rm = TRUE))
  fn_b <- function(dt) list(metric_b = sd(dt$x, na.rm = TRUE))

  se_out <- se |>
    apply_custom_fit(fn_a, "single-agent", output_assay = "out_a",
                     fit_source = "src_a") |>
    apply_custom_fit(fn_b, "single-agent", output_assay = "out_b",
                     fit_source = "src_b")

  expect_true("out_a" %in% SummarizedExperiment::assayNames(se_out))
  expect_true("out_b" %in% SummarizedExperiment::assayNames(se_out))
})

test_that("apply_custom_fit summary_fn writes to summary_assay with one row per cell", {
  built <- .build_test_se(drug_ids = c("D1", "D2"),
                          cl_ids = c("CL1", "CL2"),
                          seed = 14L)
  se <- built$se

  sum_fn <- function(fit_dt) {
    list(mean_x_mean = mean(fit_dt[["x_mean"]], na.rm = TRUE),
         n_slices = nrow(fit_dt))
  }

  se_out <- apply_custom_fit(
    se, simple_fit_fn, "single-agent",
    output_assay = "out",
    summary_fn = sum_fn,
    summary_assay = "out_summary",
    fit_source = "test"
  )

  expect_true("out_summary" %in% SummarizedExperiment::assayNames(se_out))
  sumdf <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "out_summary"),
    row.field = "row", column.field = "column"
  )
  # One summary row per (drug x cell line) cell — 2 x 2 = 4 cells
  expect_equal(nrow(sumdf), 4L)
  expect_true("mean_x_mean" %in% names(sumdf))
  expect_true("n_slices" %in% names(sumdf))
})

test_that("apply_custom_fit explicit slicing_cols overrides data_type default", {
  # Build SE where slicing column is "custom_group" instead of normalization_type
  dt <- data.table::CJ(
    row = "D1", column = "CL1",
    custom_group = c("A", "B"),
    Concentration = c(0.1, 1.0)
  )
  set.seed(77L)
  dt[, x := runif(.N)]
  bumpy <- BumpyMatrix::splitAsBumpyMatrix(
    dt[, c("custom_group", "Concentration", "x")],
    row = dt$row, col = dt$column
  )
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(Averaged = bumpy))

  fn <- function(d) list(n = nrow(d), grp = d$custom_group[1])
  se_out <- apply_custom_fit(
    se, fn, "single-agent",
    slicing_cols = "custom_group",
    slicing_values = c("A", "B"),
    output_assay = "grouped_out",
    fit_source = "grp_test"
  )

  df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "grouped_out"),
    row.field = "row", column.field = "column"
  )
  expect_equal(nrow(df), 2L)          # one row per group
  expect_setequal(df[["grp"]], c("A", "B"))
})

test_that("apply_custom_fit merge upsert replaces rows with matching fit_source + slice", {
  built <- .build_test_se(seed = 15L)
  se <- built$se

  fn_first  <- function(dt) list(val = 1.0)
  fn_second <- function(dt) list(val = 9.9)

  se1 <- apply_custom_fit(se, fn_first,  "single-agent",
                          output_assay = "out", fit_source = "src")
  se2 <- apply_custom_fit(se1, fn_second, "single-agent",
                          output_assay = "out", fit_source = "src", merge = "merge")

  df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se2, "out"),
    row.field = "row", column.field = "column"
  )
  expect_true(all(as.numeric(df[["val"]]) == 9.9))
})


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

test_that("bliss_fit_fn integrates with apply_custom_fit on combination data", {
  se <- .build_combo_se()
  se_out <- apply_custom_fit(
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
  expect_equal(nrow(df), 2L)
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

test_that("hss_fit_fn integrates with apply_custom_fit on combination data", {
  se <- .build_combo_se()
  se_out <- apply_custom_fit(
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
  se_out <- se |>
    apply_custom_fit(bliss_fit_fn, "combination",
                     output_assay = "custom_bliss", fit_source = "bliss") |>
    apply_custom_fit(hss_fit_fn,   "combination",
                     output_assay = "custom_hss",   fit_source = "hss")

  expect_true("custom_bliss" %in% SummarizedExperiment::assayNames(se_out))
  expect_true("custom_hss"   %in% SummarizedExperiment::assayNames(se_out))
})


####
# apply_custom_fits (single-pass multi-fit)
####

test_that("apply_custom_fits requires named fit_fns list", {
  built <- .build_test_se(seed = 30L)
  se <- built$se
  expect_error(
    apply_custom_fits(se, list(simple_fit_fn), "single-agent", fit_source = "t"),
    regexp = "names"
  )
})

test_that("apply_custom_fits produces one assay per fit_fn name", {
  built <- .build_test_se(seed = 31L)
  se <- built$se
  fn_a <- function(dt) list(metric_a = mean(dt$x, na.rm = TRUE))
  fn_b <- function(dt) list(metric_b = sd(dt$x, na.rm = TRUE))

  se_out <- apply_custom_fits(
    se,
    fit_fns = list(assay_a = fn_a, assay_b = fn_b),
    data_type = "single-agent",
    fit_source = "test"
  )
  expect_true("assay_a" %in% SummarizedExperiment::assayNames(se_out))
  expect_true("assay_b" %in% SummarizedExperiment::assayNames(se_out))
})

test_that("apply_custom_fits results match chained apply_custom_fit", {
  built <- .build_test_se(drug_ids = c("D1", "D2"), cl_ids = c("CL1", "CL2"),
                          seed = 32L)
  se <- built$se
  fn_a <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
  fn_b <- function(dt) list(x_sd = sd(dt$x, na.rm = TRUE))

  # single-pass
  se_multi <- apply_custom_fits(
    se, fit_fns = list(out_a = fn_a, out_b = fn_b),
    data_type = "single-agent", fit_source = "src"
  )

  # chained
  se_chain <- se |>
    apply_custom_fit(fn_a, "single-agent", output_assay = "out_a", fit_source = "src") |>
    apply_custom_fit(fn_b, "single-agent", output_assay = "out_b", fit_source = "src")

  extract <- function(s, nm) {
    data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(s, nm),
      row.field = "row", column.field = "column"
    ))
  }

  df_a_multi <- extract(se_multi, "out_a")
  df_a_chain <- extract(se_chain, "out_a")
  expect_equal(sort(as.numeric(df_a_multi$x_mean)),
               sort(as.numeric(df_a_chain$x_mean)))

  df_b_multi <- extract(se_multi, "out_b")
  df_b_chain <- extract(se_chain, "out_b")
  expect_equal(sort(as.numeric(df_b_multi$x_sd)),
               sort(as.numeric(df_b_chain$x_sd)))
})

test_that("apply_custom_fits multi-output pattern: one fn writes two assays", {
  built <- .build_test_se(seed = 33L)
  se <- built$se

  # fit_fn returns named list of named lists → writes two assays from one pass
  combo_fn <- function(dt) {
    list(
      assay_mean = list(x_mean = mean(dt$x, na.rm = TRUE)),
      assay_sd = list(x_sd = sd(dt$x, na.rm = TRUE))
    )
  }

  se_out <- apply_custom_fits(
    se,
    fit_fns = list(assay_mean = combo_fn, assay_sd = combo_fn),
    data_type = "single-agent",
    fit_source = "test"
  )

  expect_true("assay_mean" %in% SummarizedExperiment::assayNames(se_out))
  expect_true("assay_sd" %in% SummarizedExperiment::assayNames(se_out))

  df_mean <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "assay_mean"),
    row.field = "row", column.field = "column"
  )
  df_sd <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se_out, "assay_sd"),
    row.field = "row", column.field = "column"
  )
  expect_true("x_mean" %in% names(df_mean))
  expect_true("x_sd"   %in% names(df_sd))
})

test_that("apply_custom_fits on combination data: bliss + hss in one pass", {
  se <- .build_combo_se(seed = 34L)

  se_out <- apply_custom_fits(
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
  expect_equal(nrow(df_bliss), 2L)
  expect_equal(nrow(df_hss), 2L)
})

test_that("apply_custom_fits on_error='warn' skips failing fn without stopping others", {
  built <- .build_test_se(seed = 35L)
  se <- built$se

  good_fn <- function(dt) list(x_good = mean(dt$x, na.rm = TRUE))
  error_fn <- function(dt) stop("deliberate failure")

  expect_warning(
    se_out <- apply_custom_fits(
      se,
      fit_fns = list(good_assay = good_fn, bad_assay = error_fn),
      data_type = "single-agent",
      fit_source = "test",
      on_error = "warn"
    ),
    regexp = "deliberate failure"
  )
  expect_true("good_assay" %in% SummarizedExperiment::assayNames(se_out))
  expect_false("bad_assay" %in% SummarizedExperiment::assayNames(se_out))
})
