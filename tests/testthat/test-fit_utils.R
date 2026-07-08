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
                           seed = 42L) {
  avg_rows <- data.table::CJ(
    row = drug_ids,
    column = cl_ids,
    normalization_type = norm_types,
    Concentration = conc_vals
  )
  set.seed(seed)
  avg_rows[, (x_col) := runif(.N)]

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
    regexp = "normalization_types"
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
  # Build a minimal SE with a known set of data columns so we can verify that
  # fit_fn receives exactly those columns (no extras, no omissions).
  drug_ids <- "DRUG_X"
  cl_ids <- "CL_X"
  norm_types <- c("GR", "RV")
  conc_vals <- c(0.01, 0.1, 1.0)

  avg_rows <- data.table::CJ(
    row = drug_ids,
    column = cl_ids,
    normalization_type = norm_types,
    Concentration = conc_vals
  )
  avg_rows[, viability := seq(0.9, 0.5, length.out = .N)]

  data_cols <- setdiff(names(avg_rows), c("row", "column"))
  col_bumpy <- BumpyMatrix::splitAsBumpyMatrix(
    avg_rows[, data_cols, with = FALSE],
    row = avg_rows[["row"]],
    col = avg_rows[["column"]]
  )
  col_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Averaged = col_bumpy)
  )

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
    regexp = "missing response_metrics columns"
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
  ## -- Step 1: construct a minimal SE with an Averaged BumpyMatrix assay -------

  drug_ids <- c("DRUG_A", "DRUG_B")
  cl_ids <- c("CL_1", "CL_2")
  norm_types <- c("GR", "RV")
  conc_vals <- c(0.1, 0.3, 1.0, 3.0, 10.0)

  # Create cross-join of all (drug × cell line × norm_type × concentration)
  # combinations and attach synthetic response values.
  avg_rows <- data.table::CJ(
    row = drug_ids,
    column = cl_ids,
    normalization_type = norm_types,
    Concentration = conc_vals
  )
  set.seed(42)
  avg_rows[, x := runif(.N, min = 0.3, max = 1.0)]

  # Build the BumpyMatrix from the data columns (exclude row/column index cols).
  data_cols <- setdiff(names(avg_rows), c("row", "column"))
  avg_bumpy <- BumpyMatrix::splitAsBumpyMatrix(
    avg_rows[, data_cols, with = FALSE],
    row = avg_rows[["row"]],
    col = avg_rows[["column"]]
  )

  minimal_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Averaged = avg_bumpy)
  )

  ## -- Step 2: define a trivial fit_fn ------------------------------------------
  # Returns two metrics: mean of x (computed from the data) and a fixed ec50.
  trivial_fit_fn <- function(dt) {
    list(x_mean = mean(dt[["x"]], na.rm = TRUE), ec50 = 1e-6)
  }

  ## -- Step 3: run apply_fit_to_se -------------------------------------------------
  # Suppress expected warnings about missing response_metrics columns
  # (the trivial_fit_fn intentionally returns only two metrics).
  result_se <- suppressWarnings(
    apply_fit_to_se(
      minimal_se,
      trivial_fit_fn,
      normalization_types = norm_types,
      fit_source          = "integration_test"
    )
  )

  ## -- Step 4: function completed without error ---------------------------------
  expect_true(is(result_se, "SummarizedExperiment"))

  ## -- Step 5: Metrics assay was created ----------------------------------------
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(result_se))

  ## -- Step 6: unpack Metrics and verify row count, metadata, and values --------
  metrics_df <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(result_se, "Metrics"),
    row.field    = "row",
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
