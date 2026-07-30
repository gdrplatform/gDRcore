####
# fit_SE.timecourse — time-course growth-rate fitting
####
#
# Two-stage pipeline:
#   Stage 1 (.compute_growth_rates): per time-window lm() slope +
#            DMSO normalization → NormalizedGrowthRate per
#            (Drug, Concentration, CellLine, period).
#   Stage 2 (apply_fit): logistic fit on NormalizedGrowthRate ~ log10(Conc)
#            per (Drug, CellLine, period) → GR50/GRmax/EC50/Hill/r2.


#' fit_SE for time-course (Incucyte) screens
#'
#' Two-stage pipeline for time-course dose-response data:
#' \enumerate{
#'   \item Compute log2 growth rates per user-defined time window via a linear
#'     model (\code{lm(LogFoldChange ~ Duration)}) and normalise to the
#'     DMSO control of the reference period.
#'   \item Fit a 3-parameter Hill (logistic) curve on the
#'     \code{NormalizedGrowthRate ~ log10(Concentration)} dose-response for
#'     each \code{(Drug, CellLine, period)} combination.
#' }
#'
#' The \code{LogFoldChange} assay consumed by this function is produced by
#' \code{\link[gDRcore]{normalize_SE}} when \code{data_type = "time-course"}.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with a
#'   \code{LogFoldChange} assay (output of \code{normalize_SE} for
#'   \code{data_type = "time-course"}).
#' @param periods named list of numeric vectors of length 2; each element
#'   defines a time window \code{c(start_hour, end_hour)} used for growth-rate
#'   calculation.  At minimum \code{periods$early} must be provided.
#'   Example: \code{list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))}.
#' @param normalization_map named character vector mapping each period name to
#'   the reference period for DMSO normalization, or \code{"None"} to keep
#'   the raw growth rate.
#'   Example: \code{c(early = "None", mid = "early", late = "early")}.
#' @param lfc_assay string; name of the log-fold-change assay in \code{se}.
#'   Default \code{"LogFoldChange"}.
#' @param metrics_assay string; name of the output metrics assay.
#'   Default \code{"Metrics"}.
#' @param n_point_cutoff integer; minimum number of unique concentrations
#'   required to attempt a curve fit. Default \code{4L}.
#' @param range_conc numeric vector of length 2; concentration range for
#'   \code{x_AOC_range}. Default \code{c(5e-3, 5)}.
#' @param force_fit logical; if \code{TRUE} fit is kept even if not
#'   statistically significant. Default \code{FALSE}.
#' @param pcutoff numeric; p-value threshold for the F-test. Default \code{0.05}.
#' @param cap numeric; upper-capping fold above \code{x_0 = 1}. Default \code{0.1}.
#'
#' @return A new \code{SummarizedExperiment} where rows = (Drug, CellLine) and
#'   columns = period names.  The \code{Metrics} assay (or \code{metrics_assay})
#'   holds the fitted parameters: \code{ec50}, \code{xc50}, \code{h},
#'   \code{x_inf}, \code{x_0}, \code{r2}, \code{x_mean}, \code{x_AOC},
#'   \code{x_AOC_range}, \code{x_max}, \code{fit_type}, \code{N_conc},
#'   \code{maxlog10Concentration}, and \code{period}.
#'
#' @examples
#' \dontrun{
#' # After normalize_SE(se_tc, data_type = "time-course"):
#' periods  <- list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))
#' norm_map <- c(early = "None", mid = "early", late = "early")
#' se_fit   <- fit_SE.timecourse(se_tc, periods = periods,
#'                               normalization_map = norm_map)
#' gDRutils::convert_se_assay_to_dt(se_fit, "Metrics")
#' }
#'
#' @seealso \code{\link{fit_SE}} for single-agent fitting,
#'   \code{\link{fit_SE.combinations}} for combination screens,
#'   \code{\link{apply_fit}} for the generic fitting engine.
#'
#' @keywords runDrugResponseProcessingPipeline
#' @export
fit_SE.timecourse <- function(se,
                              periods,
                              normalization_map,
                              lfc_assay    = "LogFoldChange",
                              metrics_assay = "Metrics",
                              n_point_cutoff = 4L,
                              range_conc   = c(5e-3, 5),
                              force_fit    = FALSE,
                              pcutoff      = 0.05,
                              cap          = 0.1) {

  # --- Input validation ---
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_list(periods, min.len = 1L, types = "numeric")
  checkmate::assert_character(normalization_map, min.len = 1L)
  checkmate::assert_true(all(names(normalization_map) %in% names(periods)))
  checkmate::assert_string(lfc_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_int(n_point_cutoff, lower = 1L)
  checkmate::assert_numeric(range_conc, len = 2L)
  checkmate::assert_logical(force_fit, len = 1L)
  checkmate::assert_number(pcutoff, lower = 0, upper = 1)
  checkmate::assert_number(cap, lower = 0)
  gDRutils::validate_se_assay_name(se, lfc_assay)

  if (!"early" %in% names(periods)) {
    stop("'periods' must contain at least an 'early' entry.")
  }

  # --- Stage 1: growth rates ---
  growth_dt <- .compute_growth_rates(se, periods, normalization_map, lfc_assay)

  if (NROW(growth_dt) == 0L) {
    warning("fit_SE.timecourse: .compute_growth_rates returned 0 rows; returning original SE.")
    return(se)
  }

  # --- Build a period-level SE for apply_fit ---
  se_gr <- .growth_dt_to_se(growth_dt, se)

  # --- Stage 2: logistic fit via apply_fit ---
  fit_fn <- function(dt) {
    fit_drug_response_metrics(
      dt,
      x_col          = "NormalizedGrowthRate",
      n_point_cutoff = n_point_cutoff,
      range_conc     = range_conc,
      force_fit      = force_fit,
      pcutoff        = pcutoff,
      cap            = cap
    )
  }

  # Each column of se_gr IS one period (the BumpyMatrix col dimension).
  # Stage 2 uses the "time-course-metrics" profile (slicing by normalization_type="GR",
  # input_assay="GrowthRates", nested_cols=["concentration"]).
  se_gr <- apply_fit(
    se             = se_gr,
    fit_fn         = fit_fn,
    data_type      = "time-course-metrics",
    input_assay    = "GrowthRates",
    output_assay   = metrics_assay,
    merge          = "replace",
    on_error       = "warn",
    fit_source     = "gDR"
  )

  se_gr
}


####
# Internal helpers
####

#' Compute per-period growth rates from a LogFoldChange SE
#'
#' For each \code{(Barcode, WellRow, WellColumn)} within each time window
#' defined in \code{periods}, fits \code{lm(LogFoldChange ~ Duration)} and
#' extracts the slope (doublings / hour).  Slopes are then:
#' \enumerate{
#'   \item Aggregated (mean ± sd) per
#'     \code{(CellLineName, DrugName, Concentration, period)}.
#'   \item DMSO-normalised according to \code{normalization_map}.
#' }
#'
#' @param se \code{SummarizedExperiment} with \code{lfc_assay}.
#' @param periods named list; see \code{\link{fit_SE.timecourse}}.
#' @param normalization_map named character; see \code{\link{fit_SE.timecourse}}.
#' @param lfc_assay string; assay name.
#' @return \code{data.table} with columns:
#'   \code{CellLineName}, \code{DrugName}, \code{Concentration},
#'   \code{period}, \code{GrowthRate}, \code{sd_GrowthRate},
#'   \code{NormalizedGrowthRate}, \code{normalization_type}.
#' @keywords internal
.compute_growth_rates <- function(se, periods, normalization_map, lfc_assay) {

  all_data <- gDRutils::convert_se_assay_to_dt(se, lfc_assay)
  data.table::setDT(all_data)

  dur_col  <- gDRutils::get_env_identifiers("duration")
  conc_col <- gDRutils::get_env_identifiers("concentration")
  cl_col   <- gDRutils::get_env_identifiers("cellline_name")
  drug_col <- gDRutils::get_env_identifiers("drug_name")

  # Partner drug/conc identifiers (present in combo experiments, NA in single-agent)
  conc_col2 <- tryCatch(gDRutils::get_env_identifiers("concentration2"),
                        error = function(e) NULL)
  drug_col2 <- tryCatch(gDRutils::get_env_identifiers("drug_name2"),
                        error = function(e) NULL)

  # Well-level keys for lm grouping.
  # Include partner drug/conc so SA and combo wells are kept separate —
  # the same DrugA at the same Concentration appears in both SA rows
  # (DrugName2=NA) and combo rows (DrugName2=DrugB), which must not be pooled.
  barcode_col <- gDRutils::get_env_identifiers("barcode")[1L]
  well_cols   <- gDRutils::get_env_identifiers("well_position")
  lm_by <- intersect(
    c(barcode_col, well_cols, cl_col, drug_col, conc_col, drug_col2, conc_col2),
    names(all_data)
  )

  # --- Stage 1a: per-well lm slope per period ---
  all_rates <- data.table::rbindlist(lapply(names(periods), function(period_name) {
    cutoff <- as.numeric(periods[[period_name]])
    window <- all_data[get(dur_col) > cutoff[1L] & get(dur_col) < cutoff[2L]]
    if (NROW(window) == 0L) {
      return(data.table::data.table())
    }
    rates <- window[,
      .(rate = if (.N >= 2L) {
          stats::coef(stats::lm(LogFoldChange ~ Duration))[2L]
        } else {
          NA_real_
        }),
      by = lm_by
    ]
    rates[, period := period_name]
    rates
  }))

  if (NROW(all_rates) == 0L) return(data.table::data.table())

  # convert slope: doublings/hour → doublings/day
  all_rates[, rate := rate * 24]

  # --- Stage 1b: aggregate replicates ---
  # Include partner drug/conc in aggregation key so SA vs combo are distinct rows
  agg_by <- intersect(
    c(cl_col, drug_col, conc_col, drug_col2, conc_col2, "period"),
    names(all_rates)
  )
  av_rates <- all_rates[,
    .(GrowthRate = mean(rate, na.rm = TRUE),
      sd_GrowthRate = stats::sd(rate, na.rm = TRUE)),
    by = agg_by
  ]

  # --- Stage 1c: DMSO normalization ---
  # A row is DMSO/vehicle only when DrugName is a control tag AND either
  # Concentration == 0 or both partner slots are also absent/zero.
  dmso_tag <- gDRutils::get_env_identifiers("untreated_tag")
  is_primary_ctrl <- av_rates[[drug_col]] %in% c(dmso_tag, "DMSO", "vehicle") |
                     av_rates[[conc_col]] == 0
  # Exclude rows that have a real partner drug (combo with one arm being vehicle)
  has_partner <- if (!is.null(drug_col2) && drug_col2 %in% names(av_rates)) {
    !is.na(av_rates[[drug_col2]]) &
    !av_rates[[drug_col2]] %in% c(dmso_tag, "DMSO", "vehicle", NA_character_)
  } else {
    rep(FALSE, NROW(av_rates))
  }
  is_dmso <- is_primary_ctrl & !has_partner

  dmso_baselines <- av_rates[is_dmso, .(
    CellLineName = get(cl_col),
    period,
    rate_0 = GrowthRate
  )]

  if (NROW(dmso_baselines) == 0L) {
    warning(".compute_growth_rates: no DMSO/vehicle rows found; NormalizedGrowthRate will be NA.")
    av_rates[, NormalizedGrowthRate := NA_real_]
    av_rates[, normalization_type   := "NGR"]
    return(av_rates)
  }

  # Apply normalization_map.
  # "None" → NormalizedGrowthRate = GrowthRate (raw doublings/day, no DMSO division).
  # Any other value (e.g. "early") → NormalizedGrowthRate = GrowthRate / DMSO_rate[ref_period].
  norm_dt <- data.table::data.table(
    period      = names(normalization_map),
    norm_target = unname(normalization_map)
  )

  av_rates <- merge(av_rates, norm_dt, by = "period", all.x = TRUE)

  # Rows where normalization is requested (norm_target != "None")
  needs_norm <- av_rates$norm_target != "None"

  if (any(needs_norm, na.rm = TRUE)) {
    av_norm <- av_rates[needs_norm == TRUE]
    av_norm <- merge(
      av_norm,
      data.table::setnames(
        data.table::copy(dmso_baselines),
        c("CellLineName", "norm_period", "rate_0")
      ),
      by.x = c(cl_col, "norm_target"),
      by.y = c("CellLineName", "norm_period"),
      all.x = TRUE
    )
    av_norm[, NormalizedGrowthRate := GrowthRate / rate_0]
    av_rates <- data.table::rbindlist(
      list(
        av_norm,
        av_rates[needs_norm == FALSE][, NormalizedGrowthRate := GrowthRate][, rate_0 := NA_real_]
      ),
      fill = TRUE
    )
  } else {
    # All periods are "None" — keep raw growth rate
    av_rates[, NormalizedGrowthRate := GrowthRate]
    av_rates[, rate_0 := NA_real_]
  }

  # normalization_type = "NGR" (NormalizedGrowthRate) — not to be confused with
  # single-agent "GR" (GR value from gDRutils). Time-course has no RV/GR split.
  av_rates[, normalization_type := "NGR"]
  av_rates
}


#' Convert growth rate data.table to a period-level SummarizedExperiment
#'
#' Builds a thin SE where rows = (Drug, CellLine) and columns = period,
#' with a \code{"GrowthRates"} BumpyMatrix assay containing
#' \code{Concentration}, \code{NormalizedGrowthRate}, and
#' \code{normalization_type} per cell.  This SE is the input to
#' \code{\link{apply_fit}} in stage 2 of \code{\link{fit_SE.timecourse}}.
#'
#' @param growth_dt \code{data.table}; output of \code{.compute_growth_rates}.
#' @param se original time-course SE (used only for column name resolution).
#' @return \code{SummarizedExperiment}.
#' @keywords internal
.growth_dt_to_se <- function(growth_dt, se) {

  cl_col   <- gDRutils::get_env_identifiers("cellline_name")
  drug_col <- gDRutils::get_env_identifiers("drug_name")
  conc_col <- gDRutils::get_env_identifiers("concentration")

  # Exclude DMSO and zero-concentration rows from fitting
  conc_gt0 <- growth_dt[[conc_col]] > 0
  drug_not_ctrl <- !growth_dt[[drug_col]] %in%
    c(gDRutils::get_env_identifiers("untreated_tag"), "DMSO", "vehicle")
  fit_dt <- growth_dt[conc_gt0 & drug_not_ctrl]

  if (NROW(fit_dt) == 0L) {
    stop(".growth_dt_to_se: no treated rows with Concentration > 0 remain after filtering.")
  }

  # row = DrugName|CellLineName key; col = period
  fit_dt[, row    := paste0(get(drug_col), "|", get(cl_col))]
  fit_dt[, column := period]

  assay_cols <- c(conc_col, "NormalizedGrowthRate", "normalization_type")
  assay_dt   <- fit_dt[, assay_cols, with = FALSE]

  bm <- BumpyMatrix::splitAsBumpyMatrix(
    assay_dt,
    row = fit_dt$row,
    col = fit_dt$column
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(GrowthRates = bm)
  )
}
