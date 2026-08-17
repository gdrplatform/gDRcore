####
# Tests for fit_SE.timecourse, .compute_growth_rates, .growth_dt_to_se
####

# ---------------------------------------------------------------------------
# Helper: build a minimal LogFoldChange SE without going through create_SE
# (create_SE for time-course requires a full gDRimport-style input with ODBC
# identifiers — this helper constructs the assay directly).
# ---------------------------------------------------------------------------
.make_tc_se <- function(drugs = c("DrugA", "DrugB"),
                        concs = c(0, 0.01, 0.1, 1, 10),
                        durs = c(0, 12, 24, 36, 48, 60, 72, 84, 96),
                        cell_lines = "MCF7",
                        wells = "A",
                        noise_sd = 0.05,
                        seed = 1L) {
  set.seed(seed)
  rows <- data.table::CJ(
    DrugName = c("DMSO", drugs),
    Concentration = concs,
    Duration = durs,
    CellLineName = cell_lines,
    Barcode = "PL01",
    WellRow = wells,
    WellColumn = "1"
  )
  # Remove non-zero DMSO concentrations and zero-drug concentrations for
  # treated rows (simplify: DMSO is always Concentration=0)
  rows <- rows[!(DrugName == "DMSO" & Concentration > 0)]
  rows <- rows[!(DrugName != "DMSO" & Concentration == 0)]

  rows[, LogFoldChange := {
    base <- 0.03 * Duration                       # baseline growth
    drug_effect <- if (all(DrugName == "DMSO")) {
      0
    } else {
      -log10(1 + Concentration) * Duration * 0.01
    }
    base + drug_effect + rnorm(.N, 0, noise_sd)
  }, by = DrugName]

  # row_id = DrugName|CellLineName, col_id = Barcode|plate placeholder
  rows[, row_id := paste0(DrugName, "|", CellLineName)]
  rows[, col_id := Barcode]

  assay_cols <- c("DrugName", "Concentration", "Duration",
                  "CellLineName", "Barcode", "WellRow", "WellColumn",
                  "LogFoldChange")
  bm <- BumpyMatrix::splitAsBumpyMatrix(
    rows[, assay_cols, with = FALSE],
    row = rows$row_id,
    col = rows$col_id
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(LogFoldChange = bm)
  )
}


# ---------------------------------------------------------------------------
# Helper: combination time-course SE — single agent, two doublet arms
# (same partner drug at two concentrations) and a triplet arm, plus two
# distinct control tags ("vehicle" and "untreated") so that several control
# rows exist per (CellLine, period).
# ---------------------------------------------------------------------------
.make_tc_combo_se <- function(concs = c(0.01, 0.1, 1, 10),
                              durs = c(0, 24, 48, 72),
                              wells = c("A", "B"),
                              noise_sd = 0,
                              seed = 2L) {
  set.seed(seed)
  ctrl_tags <- gDRutils::get_env_identifiers("untreated_tag")
  vehicle <- ctrl_tags[1L]

  arms <- data.table::rbindlist(list(
    data.table::data.table(arm = "sa", DrugName_2 = vehicle, Concentration_2 = 0,
                           DrugName_3 = vehicle, Concentration_3 = 0),
    data.table::data.table(arm = "doublet_lo", DrugName_2 = "DrugB", Concentration_2 = 0.4,
                           DrugName_3 = vehicle, Concentration_3 = 0),
    data.table::data.table(arm = "doublet_hi", DrugName_2 = "DrugB", Concentration_2 = 0.8,
                           DrugName_3 = vehicle, Concentration_3 = 0),
    data.table::data.table(arm = "triplet", DrugName_2 = "DrugB", Concentration_2 = 0.4,
                           DrugName_3 = "DrugC", Concentration_3 = 0.2)
  ))

  grid <- data.table::CJ(arm = arms$arm, Concentration = concs, Duration = durs,
                         WellRow = wells, sorted = FALSE)
  treated <- merge(grid, arms, by = "arm")
  treated[, DrugName := "DrugA"]

  # One control block per control tag — both are vehicle-like rows that must
  # collapse into a single baseline per (CellLine, period).
  controls <- data.table::CJ(DrugName = ctrl_tags, Duration = durs,
                             WellRow = wells, sorted = FALSE)
  controls[, c("arm", "Concentration", "DrugName_2", "Concentration_2",
               "DrugName_3", "Concentration_3") :=
             .("control", 0, vehicle, 0, vehicle, 0)]

  rows <- data.table::rbindlist(list(treated, controls), use.names = TRUE)
  rows[, CellLineName := "MCF7"]
  rows[, Barcode := "PL01"]
  rows[, WellColumn := "1"]

  # Growth slows with the dose of every drug present in the well.
  rows[, LogFoldChange :=
         0.03 * Duration -
         0.01 * Duration * (log10(1 + Concentration) +
                            0.5 * log10(1 + Concentration_2) +
                            0.5 * log10(1 + Concentration_3)) +
         rnorm(.N, 0, noise_sd)]

  rows[, row_id := paste(DrugName, DrugName_2, DrugName_3, CellLineName, sep = "|")]
  rows[, col_id := Barcode]

  assay_cols <- c("DrugName", "Concentration", "DrugName_2", "Concentration_2",
                  "DrugName_3", "Concentration_3", "Duration", "CellLineName",
                  "Barcode", "WellRow", "WellColumn", "LogFoldChange")
  bm <- BumpyMatrix::splitAsBumpyMatrix(
    rows[, assay_cols, with = FALSE],
    row = rows$row_id,
    col = rows$col_id
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(LogFoldChange = bm)
  )
}


# ---------------------------------------------------------------------------
# Helper: noiseless SE whose LogFoldChange is quadratic in Duration, so the
# lm() slope depends on exactly which timepoints fall into a window.
# ---------------------------------------------------------------------------
.make_tc_curved_se <- function(durs, curvature = 5e-4) {
  rows <- data.table::CJ(
    DrugName = c("vehicle", "DrugA"),
    Duration = durs,
    CellLineName = "MCF7",
    Barcode = "PL01",
    WellRow = "A",
    WellColumn = "1"
  )
  rows[, Concentration := ifelse(DrugName == "vehicle", 0, 1)]
  rows[, LogFoldChange := curvature * Duration^2]
  rows[, row_id := paste0(DrugName, "|", CellLineName)]
  rows[, col_id := Barcode]

  assay_cols <- c("DrugName", "Concentration", "Duration", "CellLineName",
                  "Barcode", "WellRow", "WellColumn", "LogFoldChange")
  bm <- BumpyMatrix::splitAsBumpyMatrix(
    rows[, assay_cols, with = FALSE],
    row = rows$row_id,
    col = rows$col_id
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(LogFoldChange = bm)
  )
}

# Slope (doublings/day) of the quadratic curve over the given timepoints.
.expected_curved_rate <- function(durs, curvature = 5e-4) {
  y <- curvature * durs^2
  unname(stats::coef(stats::lm(y ~ durs))[2L]) * 24
}


se_tc_small <- .make_tc_se()
# Periods are half-open [start, end).
# With durs = c(0,12,24,36,48,60,72,84,96):
#   early: [0, 48)  → 0, 12, 24, 36    ≥ 2 pts per group ✓
#   late:  [48, 96) → 48, 60, 72, 84   ≥ 2 pts per group ✓
periods_std <- list(early = c(0, 48), late = c(48, 96))
norm_map_std <- c(early = "None", late = "early")


# ---------------------------------------------------------------------------
# .compute_growth_rates
# ---------------------------------------------------------------------------

test_that(".compute_growth_rates returns expected columns", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  expected_cols <- c("CellLineName", "DrugName", "Concentration", "period",
                     "GrowthRate", "sd_GrowthRate",
                     "NormalizedGrowthRate", "normalization_type")
  expect_true(all(expected_cols %in% names(gr)),
              info = paste("Missing:", toString(setdiff(expected_cols, names(gr)))))
})


test_that(".compute_growth_rates returns one row per (Drug, Conc, CellLine, period)", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  key_cols <- c("DrugName", "Concentration", "CellLineName", "period")
  expect_equal(NROW(gr), data.table::uniqueN(gr, by = key_cols))
})


test_that(".compute_growth_rates covers all requested periods", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  expect_setequal(unique(gr$period), names(periods_std))
})


test_that(".compute_growth_rates: NormalizedGrowthRate is numeric, finite for treated rows", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  treated <- gr[Concentration > 0 & DrugName != "DMSO"]
  expect_true(is.numeric(treated$NormalizedGrowthRate))
  expect_true(all(is.finite(treated$NormalizedGrowthRate)))
})


test_that(".compute_growth_rates 'None' keeps raw GrowthRate, not DMSO-normalised", {
  # early = "None" → NormalizedGrowthRate == GrowthRate (no division by DMSO)
  # late  = "early" → NormalizedGrowthRate == GrowthRate / DMSO_rate[early]
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  early_dmso <- gr[period == "early" & DrugName == "DMSO"]
  early_drug <- gr[period == "early" & DrugName == "DrugA" & Concentration == 1]
  # For "None" period: NormalizedGrowthRate must equal GrowthRate exactly
  expect_equal(early_dmso$NormalizedGrowthRate, early_dmso$GrowthRate)
  expect_equal(early_drug$NormalizedGrowthRate, early_drug$GrowthRate)

  # For "early"-referenced late period: NGR = GrowthRate / DMSO_early_rate
  late_drug <- gr[period == "late" & DrugName == "DrugA" & Concentration == 1]
  dmso_early_rate <- gr[period == "early" & DrugName == "DMSO"]$GrowthRate
  expect_equal(late_drug$NormalizedGrowthRate,
               late_drug$GrowthRate / dmso_early_rate,
               tolerance = 1e-10)
})


test_that(".compute_growth_rates returns empty data.table when no rows in window", {
  # Period range that has no data
  empty_periods <- list(impossible = c(1000, 2000))
  empty_map <- c(impossible = "None")
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, empty_periods, empty_map, "LogFoldChange"
  )
  expect_equal(NROW(gr), 0L)
})


test_that(".compute_growth_rates windows are half-open [start, end)", {
  durs <- c(0, 12, 24, 36)
  se_curved <- .make_tc_curved_se(durs)
  gr <- gDRcore:::.compute_growth_rates(
    se_curved,
    list(early = c(0, 24), late = c(24, 48)),
    c(early = "None", late = "None"),
    "LogFoldChange"
  )
  treated <- gr[DrugName == "DrugA"]
  # early takes {0, 12} (24 belongs to the next window), late takes {24, 36}
  expect_equal(treated[period == "early"]$GrowthRate,
               .expected_curved_rate(c(0, 12)), tolerance = 1e-10)
  expect_equal(treated[period == "late"]$GrowthRate,
               .expected_curved_rate(c(24, 36)), tolerance = 1e-10)
  # Contiguous windows partition the timepoints: no timepoint is used twice
  # and none is dropped, so the two slopes differ on a curved profile.
  expect_false(isTRUE(all.equal(treated[period == "early"]$GrowthRate,
                                treated[period == "late"]$GrowthRate)))
})


test_that(".compute_growth_rates aggregates replicate wells", {
  se_rep <- .make_tc_se(wells = c("A", "B"), noise_sd = 0)
  gr <- gDRcore:::.compute_growth_rates(
    se_rep, periods_std, norm_map_std, "LogFoldChange"
  )
  # Two noiseless replicate wells → sd is defined and equal to zero
  expect_true(all(is.finite(gr$sd_GrowthRate)))
  expect_equal(unique(gr$sd_GrowthRate), 0, tolerance = 1e-12)

  # A single well cannot yield a standard deviation
  gr_single <- gDRcore:::.compute_growth_rates(
    .make_tc_se(wells = "A", noise_sd = 0), periods_std, norm_map_std, "LogFoldChange"
  )
  expect_true(all(is.na(gr_single$sd_GrowthRate)))
})


# ---------------------------------------------------------------------------
# .compute_growth_rates — combination layouts
# ---------------------------------------------------------------------------

test_that(".compute_growth_rates keeps single-agent, doublet and triplet arms separate", {
  se_combo <- .make_tc_combo_se()
  gr <- gDRcore:::.compute_growth_rates(
    se_combo, periods_std, norm_map_std, "LogFoldChange"
  )
  key_cols <- c("DrugName", "Concentration", "DrugName_2", "Concentration_2",
                "DrugName_3", "Concentration_3", "CellLineName", "period")
  expect_true(all(key_cols %in% names(gr)))
  expect_equal(NROW(gr), data.table::uniqueN(gr, by = key_cols))

  # DrugA at 1 uM in the early period: SA + two doublet arms + triplet
  arms <- gr[DrugName == "DrugA" & Concentration == 1 & period == "early"]
  expect_equal(NROW(arms), 4L)
  expect_setequal(arms$Concentration_2, c(0, 0.4, 0.8, 0.4))
  # Doublet arms differ only by partner concentration and must not be pooled
  expect_equal(NROW(arms[DrugName_2 == "DrugB" & DrugName_3 == "vehicle"]), 2L)
})


test_that(".compute_growth_rates collapses several control rows into one baseline", {
  # "vehicle" and "untreated" rows both qualify as controls; without
  # aggregation the baseline merge would be many-to-many and duplicate rows.
  se_combo <- .make_tc_combo_se()
  gr <- gDRcore:::.compute_growth_rates(
    se_combo, periods_std, norm_map_std, "LogFoldChange"
  )
  key_cols <- c("DrugName", "Concentration", "DrugName_2", "Concentration_2",
                "DrugName_3", "Concentration_3", "CellLineName", "period")
  late_treated <- gr[period == "late" & DrugName == "DrugA"]
  expect_equal(NROW(late_treated), data.table::uniqueN(late_treated, by = key_cols))
  expect_true(all(is.finite(late_treated$NormalizedGrowthRate)))

  # NGR = GrowthRate / mean(control rates of the reference period)
  ctrl_rate <- mean(gr[period == "early" &
                         DrugName %in% gDRutils::get_env_identifiers("untreated_tag")]$GrowthRate)
  expect_equal(late_treated$NormalizedGrowthRate,
               late_treated$GrowthRate / ctrl_rate,
               tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# .growth_dt_to_se
# ---------------------------------------------------------------------------

test_that(".growth_dt_to_se returns a SummarizedExperiment with GrowthRates assay", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  se_gr <- gDRcore:::.growth_dt_to_se(gr)
  expect_s4_class(se_gr, "SummarizedExperiment")
  expect_true("GrowthRates" %in% SummarizedExperiment::assayNames(se_gr))
})


test_that(".growth_dt_to_se columns = period names", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  se_gr <- gDRcore:::.growth_dt_to_se(gr)
  expect_setequal(colnames(se_gr), names(periods_std))
})


test_that(".growth_dt_to_se GrowthRates cells contain NormalizedGrowthRate column", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  se_gr <- gDRcore:::.growth_dt_to_se(gr)
  bm <- SummarizedExperiment::assay(se_gr, "GrowthRates")
  pk <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(bm))
  expect_true("NormalizedGrowthRate" %in% names(pk))
  expect_true("Concentration" %in% names(pk))
})


test_that(".growth_dt_to_se gives every combo arm its own row", {
  se_combo <- .make_tc_combo_se()
  gr <- gDRcore:::.compute_growth_rates(
    se_combo, periods_std, norm_map_std, "LogFoldChange"
  )
  se_gr <- gDRcore:::.growth_dt_to_se(gr)
  # SA + doublet@0.4 + doublet@0.8 + triplet
  expect_equal(NROW(se_gr), 4L)
  expect_true(all(grepl("^DrugA\\|", rownames(se_gr))))
  # Partner concentration is part of the row key and of the cell content
  expect_true(any(grepl("0\\.8", rownames(se_gr))))
  bm <- SummarizedExperiment::assay(se_gr, "GrowthRates")
  pk <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(bm))
  expect_true("Concentration_2" %in% names(pk))
})


# ---------------------------------------------------------------------------
# fit_SE.timecourse — full pipeline
# ---------------------------------------------------------------------------

test_that("fit_SE.timecourse returns a SummarizedExperiment with Metrics assay", {
  se_fit <- suppressWarnings(fit_SE.timecourse(
    se_tc_small,
    periods = periods_std,
    normalization_map = norm_map_std
  ))
  expect_s4_class(se_fit, "SummarizedExperiment")
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(se_fit))
})


test_that("fit_SE.timecourse Metrics contains expected fit columns", {
  se_fit <- suppressWarnings(fit_SE.timecourse(
    se_tc_small,
    periods = periods_std,
    normalization_map = norm_map_std
  ))
  metrics_dt <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se_fit, "Metrics"),
      row.field = "row", column.field = "column"
    )
  )
  expected_fit_cols <- c("ec50", "xc50", "h", "r2", "x_mean",
                         "x_AOC", "N_conc", "fit_type")
  expect_true(all(expected_fit_cols %in% names(metrics_dt)),
              info = paste("Missing:", toString(setdiff(expected_fit_cols, names(metrics_dt)))))
})


test_that("fit_SE.timecourse columns = period names", {
  se_fit <- suppressWarnings(fit_SE.timecourse(
    se_tc_small,
    periods = periods_std,
    normalization_map = norm_map_std
  ))
  expect_setequal(colnames(se_fit), names(periods_std))
})


test_that("fit_SE.timecourse: custom metrics_assay name is honoured", {
  se_fit <- suppressWarnings(fit_SE.timecourse(
    se_tc_small,
    periods = periods_std,
    normalization_map = norm_map_std,
    metrics_assay = "TC_Metrics"
  ))
  expect_true("TC_Metrics" %in% SummarizedExperiment::assayNames(se_fit))
  expect_false("Metrics" %in% SummarizedExperiment::assayNames(se_fit))
})


test_that("fit_SE.timecourse is deterministic", {
  args <- list(se_tc_small, periods = periods_std, normalization_map = norm_map_std)
  first <- suppressWarnings(do.call(fit_SE.timecourse, args))
  second <- suppressWarnings(do.call(fit_SE.timecourse, args))
  as_dt <- function(se) {
    data.table::as.data.table(
      BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se, "Metrics"))
    )
  }
  expect_equal(as_dt(first), as_dt(second))
})


test_that("fit_SE.timecourse fits every combo arm", {
  se_combo <- .make_tc_combo_se()
  se_fit <- suppressWarnings(fit_SE.timecourse(
    se_combo,
    periods = periods_std,
    normalization_map = norm_map_std
  ))
  expect_equal(NROW(se_fit), 4L)
  expect_setequal(colnames(se_fit), names(periods_std))
  metrics_dt <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se_fit, "Metrics"),
      row.field = "row", column.field = "column"
    )
  )
  expect_equal(data.table::uniqueN(metrics_dt, by = c("row", "column")), NROW(metrics_dt))
})


test_that("fit_SE.timecourse errors when 'early' period is missing", {
  bad_periods <- list(mid = c(24, 72))
  bad_norm_map <- c(mid = "None")
  expect_error(
    fit_SE.timecourse(se_tc_small, periods = bad_periods,
                      normalization_map = bad_norm_map),
    "'periods' must contain at least an 'early' entry"
  )
})


test_that("fit_SE.timecourse errors on missing LogFoldChange assay", {
  se_no_lfc <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Averaged = SummarizedExperiment::assay(se_tc_small, "LogFoldChange"))
  )
  expect_error(
    fit_SE.timecourse(se_no_lfc, periods = periods_std,
                      normalization_map = norm_map_std),
    regexp = "LogFoldChange|assay"
  )
})


test_that("fit_SE.timecourse three periods produces three-column SE", {
  se_tc3 <- .make_tc_se(durs = c(0, 24, 48, 72, 96, 120, 144))
  periods3 <- list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))
  norm_map3 <- c(early = "None", mid = "early", late = "early")
  se_fit <- fit_SE.timecourse(se_tc3, periods = periods3,
                              normalization_map = norm_map3)
  expect_setequal(colnames(se_fit), c("early", "mid", "late"))
})


# ---------------------------------------------------------------------------
# fit_SE.timecourse — input validation
# ---------------------------------------------------------------------------

test_that("fit_SE.timecourse requires normalization_map to cover every period", {
  # 'late' has no entry — its rows would otherwise vanish from the output
  expect_error(
    fit_SE.timecourse(se_tc_small, periods = periods_std,
                      normalization_map = c(early = "None")),
    regexp = "equal"
  )
})


test_that("fit_SE.timecourse rejects a normalization target that is not a period", {
  expect_error(
    fit_SE.timecourse(se_tc_small, periods = periods_std,
                      normalization_map = c(early = "None", late = "middle")),
    regexp = "subset"
  )
})


test_that("fit_SE.timecourse rejects an unnamed normalization_map", {
  expect_error(
    fit_SE.timecourse(se_tc_small, periods = periods_std,
                      normalization_map = c("None", "early")),
    regexp = "names|named"
  )
})


test_that("fit_SE.timecourse rejects a non-increasing period window", {
  expect_error(
    fit_SE.timecourse(se_tc_small, periods = list(early = c(48, 48)),
                      normalization_map = c(early = "None")),
    regexp = "TRUE"
  )
})


test_that("fit_SE.timecourse rejects a non-increasing range_conc", {
  expect_error(
    fit_SE.timecourse(se_tc_small, periods = periods_std,
                      normalization_map = norm_map_std,
                      range_conc = c(5, 5e-3)),
    regexp = "TRUE"
  )
})


test_that("fit_SE.timecourse errors when no timepoint falls into any window", {
  expect_error(
    fit_SE.timecourse(se_tc_small, periods = list(early = c(1000, 2000)),
                      normalization_map = c(early = "None")),
    regexp = "no growth rates could be computed"
  )
})
