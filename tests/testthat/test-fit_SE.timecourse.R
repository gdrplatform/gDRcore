####
# Tests for fit_SE.timecourse, .compute_growth_rates, .growth_dt_to_se
####

# ---------------------------------------------------------------------------
# Helper: build a minimal LogFoldChange SE without going through create_SE
# (create_SE for time-course requires a full gDRimport-style input with ODBC
# identifiers — this helper constructs the assay directly).
# ---------------------------------------------------------------------------
.make_tc_se <- function(drugs      = c("DrugA", "DrugB"),
                        concs      = c(0, 0.01, 0.1, 1, 10),
                        durs       = c(0, 12, 24, 36, 48, 60, 72, 84, 96),
                        cell_lines = "MCF7",
                        seed       = 1L) {
  set.seed(seed)
  rows <- data.table::CJ(
    DrugName     = c("DMSO", drugs),
    Concentration = concs,
    Duration     = durs,
    CellLineName = cell_lines,
    Barcode      = "PL01",
    WellRow      = "A",
    WellColumn   = "1"
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
    base + drug_effect + rnorm(.N, 0, 0.05)
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

se_tc_small <- .make_tc_se()
# Periods use open intervals (Duration > lo & < hi).
# With durs = c(0,12,24,36,48,60,72,84,96):
#   early: (0, 48)  → 12, 24, 36     ≥ 2 pts per group ✓
#   late:  (48, 96) → 60, 72, 84     ≥ 2 pts per group ✓
periods_std  <- list(early = c(0, 48), late = c(48, 96))
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
  early_drug  <- gr[period == "early" & DrugName == "DrugA" & Concentration == 1]
  # For "None" period: NormalizedGrowthRate must equal GrowthRate exactly
  expect_equal(early_dmso$NormalizedGrowthRate, early_dmso$GrowthRate)
  expect_equal(early_drug$NormalizedGrowthRate,  early_drug$GrowthRate)

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
  empty_map     <- c(impossible = "None")
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, empty_periods, empty_map, "LogFoldChange"
  )
  expect_equal(NROW(gr), 0L)
})


# ---------------------------------------------------------------------------
# .growth_dt_to_se
# ---------------------------------------------------------------------------

test_that(".growth_dt_to_se returns a SummarizedExperiment with GrowthRates assay", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  se_gr <- gDRcore:::.growth_dt_to_se(gr, se_tc_small)
  expect_s4_class(se_gr, "SummarizedExperiment")
  expect_true("GrowthRates" %in% SummarizedExperiment::assayNames(se_gr))
})


test_that(".growth_dt_to_se columns = period names", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  se_gr <- gDRcore:::.growth_dt_to_se(gr, se_tc_small)
  expect_setequal(colnames(se_gr), names(periods_std))
})


test_that(".growth_dt_to_se GrowthRates cells contain NormalizedGrowthRate column", {
  gr <- gDRcore:::.compute_growth_rates(
    se_tc_small, periods_std, norm_map_std, "LogFoldChange"
  )
  se_gr <- gDRcore:::.growth_dt_to_se(gr, se_tc_small)
  bm    <- SummarizedExperiment::assay(se_gr, "GrowthRates")
  pk    <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(bm))
  expect_true("NormalizedGrowthRate" %in% names(pk))
  expect_true("Concentration" %in% names(pk))
})


# ---------------------------------------------------------------------------
# fit_SE.timecourse — full pipeline
# ---------------------------------------------------------------------------

test_that("fit_SE.timecourse returns a SummarizedExperiment with Metrics assay", {
  se_fit <- fit_SE.timecourse(
    se_tc_small,
    periods          = periods_std,
    normalization_map = norm_map_std
  )
  expect_s4_class(se_fit, "SummarizedExperiment")
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(se_fit))
})


test_that("fit_SE.timecourse Metrics contains expected fit columns", {
  se_fit <- fit_SE.timecourse(
    se_tc_small,
    periods          = periods_std,
    normalization_map = norm_map_std
  )
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
  se_fit <- fit_SE.timecourse(
    se_tc_small,
    periods          = periods_std,
    normalization_map = norm_map_std
  )
  expect_setequal(colnames(se_fit), names(periods_std))
})


test_that("fit_SE.timecourse: custom metrics_assay name is honoured", {
  se_fit <- fit_SE.timecourse(
    se_tc_small,
    periods          = periods_std,
    normalization_map = norm_map_std,
    metrics_assay    = "TC_Metrics"
  )
  expect_true("TC_Metrics" %in% SummarizedExperiment::assayNames(se_fit))
  expect_false("Metrics" %in% SummarizedExperiment::assayNames(se_fit))
})


test_that("fit_SE.timecourse errors when 'early' period is missing", {
  bad_periods  <- list(mid = c(24, 72))
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
  periods3  <- list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))
  norm_map3 <- c(early = "None", mid = "early", late = "early")
  se_fit <- fit_SE.timecourse(se_tc3, periods = periods3,
                              normalization_map = norm_map3)
  expect_setequal(colnames(se_fit), c("early", "mid", "late"))
})
