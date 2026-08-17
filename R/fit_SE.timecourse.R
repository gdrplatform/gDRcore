####
# fit_SE.timecourse — time-course growth-rate fitting
####
#
# Two-stage pipeline:
#   Stage 1 (compute_growth_rates): per time-window growth rate +
#            DMSO normalization → NormalizedGrowthRate per
#            (Drug, Concentration, CellLine, period).
#   Stage 2 (apply_fit): logistic fit on NormalizedGrowthRate ~ log10(Conc)
#            per (Drug, CellLine, period) → GR50/GRmax/EC50/Hill/r2.
#
# Both stages are pluggable: `rate_fn` replaces the growth-rate estimator and
# `fit_fn` replaces the dose-response fit, mirroring the custom fit interface
# available for single-agent and combination data.


#' fit_SE for time-course (Incucyte) screens
#'
#' Two-stage pipeline for time-course dose-response data:
#' \enumerate{
#'   \item Compute growth rates per user-defined time window and normalise them
#'     to the control of the reference period (\code{\link{compute_growth_rates}}).
#'   \item Fit a 3-parameter Hill (logistic) curve on the
#'     \code{NormalizedGrowthRate ~ log10(Concentration)} dose-response for
#'     each \code{(Drug, CellLine, period)} combination.
#' }
#'
#' The \code{LogFoldChange} assay consumed by this function is produced by
#' \code{\link[gDRcore]{normalize_SE}} when \code{data_type = "time-course"}.
#'
#' Both stages accept a replacement function, so custom fitting works for
#' time-course data the same way it does for single-agent and combination
#' screens.  The two stages can also be run separately:
#' \code{\link{compute_growth_rates}} returns the intermediate
#' \code{GrowthRates} SummarizedExperiment, which can be passed to
#' \code{\link{apply_fit}} with \code{data_type = "time-course-metrics"}.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with a
#'   \code{LogFoldChange} assay (output of \code{normalize_SE} for
#'   \code{data_type = "time-course"}).
#' @param periods named list of numeric vectors of length 2; each element
#'   defines a half-open time window \code{c(start_hour, end_hour)} used for
#'   growth-rate calculation: a timepoint is used when
#'   \code{start_hour <= Duration < end_hour}, so contiguous windows partition
#'   the time-course without dropping the shared boundary.
#'   At minimum \code{periods$early} must be provided.
#'   Example: \code{list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))}.
#' @param normalization_map named character vector mapping \emph{every} period
#'   name to the reference period for DMSO normalization, or \code{"None"} to
#'   keep the raw growth rate.  Its names must match \code{names(periods)}
#'   exactly and its values must be \code{"None"} or an existing period name.
#'   Example: \code{c(early = "None", mid = "early", late = "early")}.
#' @param lfc_assay string; name of the log-fold-change assay in \code{se}.
#'   Defaults to the \code{input_assay} of the \code{"time-course"} fit profile.
#' @param metrics_assay string; name of the output metrics assay.
#'   Default \code{"Metrics"}.
#' @param rate_fn function or \code{NULL}; growth-rate estimator passed to
#'   \code{\link{compute_growth_rates}}.  \code{NULL} uses the default
#'   \code{lm(LogFoldChange ~ Duration)} slope.
#' @param fit_fn function or \code{NULL}; dose-response fit applied to each
#'   \code{normalization_type} slice of the \code{GrowthRates} assay.
#'   \code{NULL} uses \code{\link{fit_drug_response_metrics}} on
#'   \code{NormalizedGrowthRate}.  When a custom function is supplied,
#'   \code{n_point_cutoff}, \code{range_conc}, \code{force_fit},
#'   \code{pcutoff} and \code{cap} are not used.
#' @param fit_source string or \code{NULL}; value recorded in the
#'   \code{fit_source} column.  \code{NULL} records \code{"gDR"} for the default
#'   fit and \code{"custom"} for a user-supplied \code{fit_fn}.
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
#'   columns = period names, holding the \code{GrowthRates} assay from stage 1
#'   and the \code{Metrics} assay (or \code{metrics_assay}) from stage 2.  With
#'   the default \code{fit_fn} the metrics are: \code{ec50}, \code{xc50},
#'   \code{h}, \code{x_inf}, \code{x_0}, \code{r2}, \code{x_mean}, \code{x_AOC},
#'   \code{x_AOC_range}, \code{x_max}, \code{fit_type}, \code{N_conc},
#'   \code{maxlog10Concentration}, and \code{period}.
#'
#' @examples
#' \dontrun{
#' # After normalize_SE(se_tc, data_type = "time-course"):
#' periods <- list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))
#' norm_map <- c(early = "None", mid = "early", late = "early")
#' se_fit <- fit_SE.timecourse(se_tc, periods = periods,
#'                             normalization_map = norm_map)
#' gDRutils::convert_se_assay_to_dt(se_fit, "Metrics")
#'
#' # Custom growth-rate estimator and custom dose-response fit:
#' se_custom <- fit_SE.timecourse(
#'   se_tc,
#'   periods = periods,
#'   normalization_map = norm_map,
#'   rate_fn = function(dt) 24 * stats::median(diff(dt$LogFoldChange) /
#'                                             diff(dt$Duration)),
#'   fit_fn = function(dt) data.table::data.table(x_mean = mean(dt$NormalizedGrowthRate))
#' )
#' }
#'
#' @seealso \code{\link{compute_growth_rates}} for stage 1 on its own,
#'   \code{\link{fit_SE}} for single-agent fitting,
#'   \code{\link{fit_SE.combinations}} for combination screens,
#'   \code{\link{apply_fit}} for the generic fitting engine.
#'
#' @keywords runDrugResponseProcessingPipeline
#' @export
fit_SE.timecourse <- function(se,
                              periods,
                              normalization_map,
                              lfc_assay = NULL,
                              metrics_assay = "Metrics",
                              rate_fn = NULL,
                              fit_fn = NULL,
                              fit_source = NULL,
                              n_point_cutoff = 4L,
                              range_conc = c(5e-3, 5),
                              force_fit = FALSE,
                              pcutoff = 0.05,
                              cap = 0.1) {

  # --- Input validation (stage 1 arguments are validated by compute_growth_rates) ---
  checkmate::assert_string(metrics_assay)
  checkmate::assert_function(fit_fn, null.ok = TRUE)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_int(n_point_cutoff, lower = 1L)
  checkmate::assert_numeric(range_conc, len = 2L, lower = 0, any.missing = FALSE,
                            finite = TRUE)
  checkmate::assert_true(range_conc[1L] < range_conc[2L])
  checkmate::assert_logical(force_fit, len = 1L)
  checkmate::assert_number(pcutoff, lower = 0, upper = 1)
  checkmate::assert_number(cap, lower = 0)

  # --- Stage 1: growth rates ---
  se_gr <- compute_growth_rates(
    se,
    periods = periods,
    normalization_map = normalization_map,
    lfc_assay = lfc_assay,
    rate_fn = rate_fn
  )

  # --- Stage 2: dose-response fit via apply_fit ---
  if (is.null(fit_fn)) {
    fit_fn <- function(dt) {
      fit_drug_response_metrics(
        dt,
        x_col = "NormalizedGrowthRate",
        n_point_cutoff = n_point_cutoff,
        range_conc = range_conc,
        force_fit = force_fit,
        pcutoff = pcutoff,
        cap = cap
      )
    }
    default_source <- "gDR"
  } else {
    default_source <- "custom"
  }
  if (is.null(fit_source)) {
    fit_source <- default_source
  }

  # Each column of se_gr IS one period (the BumpyMatrix col dimension).
  # Stage 2 uses the "time-course-metrics" profile (slicing by normalization_type="NGR",
  # input_assay="GrowthRates", nested_cols=["concentration"]).
  apply_fit(
    se = se_gr,
    fit_fn = fit_fn,
    data_type = "time-course-metrics",
    output_assay = metrics_assay,
    merge = "replace",
    on_error = "warn",
    fit_source = fit_source
  )
}


#' Compute per-period growth rates for time-course screens
#'
#' Stage 1 of the time-course pipeline: fits a growth rate per well within each
#' time window, aggregates replicates, and normalises the rates to the control
#' wells of a reference period.  The result is the input of stage 2 — pass it to
#' \code{\link{apply_fit}} with \code{data_type = "time-course-metrics"} to fit
#' dose-response curves with an arbitrary fit function, or use
#' \code{\link{fit_SE.timecourse}} to run both stages at once.
#'
#' Growth rates are computed per \code{(Barcode, WellRow, WellColumn, CellLine,
#' Drug, Concentration, partner drug/concentration slots)} group, so
#' single-agent, doublet and triplet arms are never pooled.  Control wells are
#' detected via \code{gDRutils::get_env_identifiers("untreated_tag")} and a
#' single baseline per \code{(CellLine, period)} is used for normalization.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with a
#'   log-fold-change assay.
#' @param periods named list of numeric vectors of length 2; half-open time
#'   windows \code{c(start_hour, end_hour)}, see \code{\link{fit_SE.timecourse}}.
#' @param normalization_map named character vector mapping every period to its
#'   reference period or to \code{"None"}, see \code{\link{fit_SE.timecourse}}.
#' @param lfc_assay string or \code{NULL}; name of the log-fold-change assay.
#'   \code{NULL} uses the \code{input_assay} of the \code{"time-course"} fit
#'   profile (see \code{\link{get_fit_profile}}).
#' @param rate_fn function or \code{NULL}; growth-rate estimator.  It receives a
#'   \code{data.table} holding the timepoints of one well within one window —
#'   including the duration column and the log-fold-change column — and must
#'   return a single numeric value: the growth rate in doublings per day.
#'   \code{NULL} uses the slope of \code{lm(LogFoldChange ~ Duration)}.
#'
#' @return A \code{SummarizedExperiment} where rows = (Drug, CellLine, combo arm)
#'   and columns = period names, with a \code{GrowthRates} assay holding
#'   \code{Concentration}, the partner concentrations, \code{NormalizedGrowthRate}
#'   and \code{normalization_type} per cell.
#'
#' @examples
#' \dontrun{
#' periods <- list(early = c(0, 48), late = c(48, 96))
#' norm_map <- c(early = "None", late = "early")
#' se_gr <- compute_growth_rates(se_tc, periods, norm_map)
#'
#' # Stage 2 with any custom fit function:
#' se_fit <- apply_fit(
#'   se_gr,
#'   fit_fn = my_fit_fn,
#'   data_type = "time-course-metrics",
#'   output_assay = "Metrics",
#'   merge = "replace",
#'   fit_source = "custom"
#' )
#' }
#'
#' @seealso \code{\link{fit_SE.timecourse}}, \code{\link{apply_fit}}
#'
#' @keywords runDrugResponseProcessingPipeline
#' @export
compute_growth_rates <- function(se,
                                 periods,
                                 normalization_map,
                                 lfc_assay = NULL,
                                 rate_fn = NULL) {

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_list(periods, min.len = 1L, types = "numeric", names = "named")
  checkmate::assert_character(normalization_map, min.len = 1L, names = "named",
                              any.missing = FALSE)
  # Both directions: a period absent from the map would get no normalization
  # target and its rows would silently drop out of the output.
  checkmate::assert_set_equal(names(normalization_map), names(periods))
  # Every target must be either "None" or an existing period.
  checkmate::assert_subset(unname(normalization_map), c("None", names(periods)))
  checkmate::assert_string(lfc_assay, null.ok = TRUE)
  checkmate::assert_function(rate_fn, null.ok = TRUE)

  invisible(lapply(names(periods), function(period_name) {
    window <- periods[[period_name]]
    checkmate::assert_numeric(window, len = 2L, any.missing = FALSE, finite = TRUE,
                              .var.name = sprintf("periods[['%s']]", period_name))
    checkmate::assert_true(window[1L] < window[2L])
  }))

  if (is.null(lfc_assay)) {
    lfc_assay <- get_fit_profile("time-course")$input_assay
  }
  gDRutils::validate_se_assay_name(se, lfc_assay)

  if (!"early" %in% names(periods)) {
    stop("'periods' must contain at least an 'early' entry.")
  }

  growth_dt <- .compute_growth_rates(se, periods, normalization_map, lfc_assay, rate_fn)

  if (NROW(growth_dt) == 0L) {
    stop("compute_growth_rates: no growth rates could be computed - check that 'periods' ",
         "overlap the durations present in the '", lfc_assay, "' assay.")
  }

  .growth_dt_to_se(growth_dt)
}


####
# Internal helpers
####

#' Default growth-rate estimator: slope of lm(LogFoldChange ~ Duration)
#'
#' @param dur_col string; duration column name.
#' @param lfc_col string; log-fold-change column name.
#' @return function taking a \code{data.table} and returning doublings per day.
#' @keywords internal
#' @noRd
.default_rate_fn <- function(dur_col, lfc_col) {
  lm_formula <- stats::as.formula(paste(lfc_col, "~", dur_col))
  function(dt) {
    if (NROW(dt) < 2L) {
      NA_real_
    } else {
      # the slope is in doublings per hour; downstream works in doublings per day
      unname(stats::coef(stats::lm(lm_formula, data = dt))[2L]) * 24
    }
  }
}


#' Resolve partner drug/concentration column names present in a data object
#'
#' Queries \code{gDRutils::get_env_identifiers()} for slots 2 and 3 and
#' returns only those that actually appear in \code{names(x)}.
#'
#' @param x named object (\code{data.table}, \code{data.frame}, etc.) whose
#'   column names are tested for presence.
#' @return list with elements \code{drug} and \code{conc}, each a character
#'   vector of resolved column names present in \code{x}.
#' @keywords internal
#' @noRd
.resolve_partner_cols <- function(x) {
  ids_drug <- c("drug_name2", "drug_name3")
  ids_conc <- c("concentration2", "concentration3")
  drug_cols <- unlist(lapply(ids_drug,
    function(k) tryCatch(gDRutils::get_env_identifiers(k), error = function(e) NULL)))
  conc_cols <- unlist(lapply(ids_conc,
    function(k) tryCatch(gDRutils::get_env_identifiers(k), error = function(e) NULL)))
  list(
    drug = intersect(drug_cols, names(x)),
    conc = intersect(conc_cols, names(x))
  )
}


#' Compute per-period growth rates from a LogFoldChange SE
#'
#' For each \code{(Barcode, WellRow, WellColumn)} within each half-open time
#' window defined in \code{periods}, applies \code{rate_fn} and collects the
#' resulting growth rate (doublings / day).  Rates are then:
#' \enumerate{
#'   \item Aggregated (mean ± sd) per
#'     \code{(CellLineName, DrugName, Concentration, partner slots, period)}.
#'   \item DMSO-normalised according to \code{normalization_map}.
#' }
#'
#' @param se \code{SummarizedExperiment} with \code{lfc_assay}.
#' @param periods named list; see \code{\link{fit_SE.timecourse}}.
#' @param normalization_map named character; see \code{\link{fit_SE.timecourse}}.
#' @param lfc_assay string; assay name.
#' @param rate_fn function or \code{NULL}; see \code{\link{compute_growth_rates}}.
#' @return \code{data.table} with columns:
#'   \code{CellLineName}, \code{DrugName}, \code{Concentration},
#'   \code{period}, \code{GrowthRate}, \code{sd_GrowthRate},
#'   \code{NormalizedGrowthRate}, \code{normalization_type}.
#' @keywords internal
#' @noRd
.compute_growth_rates <- function(se, periods, normalization_map, lfc_assay,
                                  rate_fn = NULL) {

  all_data <- gDRutils::convert_se_assay_to_dt(se, lfc_assay)
  data.table::setDT(all_data)

  dur_col <- gDRutils::get_env_identifiers("duration")
  conc_col <- gDRutils::get_env_identifiers("concentration")
  cl_col <- gDRutils::get_env_identifiers("cellline_name")
  drug_col <- gDRutils::get_env_identifiers("drug_name")

  if (is.null(rate_fn)) {
    rate_fn <- .default_rate_fn(dur_col, lfc_assay)
  }

  # Partner drug/conc identifiers for all slots (2, 3, ...).
  # Include every slot present in the data so that doublet (DrugName3=vehicle)
  # and triplet (DrugName3=RealDrug) combos are treated as distinct treatments
  # and not pooled together in the per-well fit or replicate aggregation.
  partner_cols <- .resolve_partner_cols(all_data)
  partner_drug_cols <- partner_cols$drug
  partner_conc_cols <- partner_cols$conc

  # Well-level keys for rate_fn grouping.
  # Include all partner drug/conc slots so SA and every combo arm
  # (doublet vs triplet, different partner concentrations) stay separate.
  barcode_col <- gDRutils::get_env_identifiers("barcode")[1L]
  well_cols <- gDRutils::get_env_identifiers("well_position")
  rate_by <- intersect(
    c(barcode_col, well_cols, cl_col, drug_col, conc_col,
      partner_drug_cols, partner_conc_cols),
    names(all_data)
  )

  # --- Stage 1a: per-well growth rate per period ---
  # Windows are half-open [start, end): a timepoint on the boundary between two
  # contiguous windows belongs to the later one, and t = 0 is not dropped.
  all_rates <- data.table::rbindlist(lapply(names(periods), function(period_name) {
    cutoff <- as.numeric(periods[[period_name]])
    window <- all_data[get(dur_col) >= cutoff[1L] & get(dur_col) < cutoff[2L]]
    if (NROW(window) == 0L) {
      return(data.table::data.table())
    }
    rates <- window[, .(rate = rate_fn(.SD)), by = rate_by]
    rates[, period := period_name]
    rates
  }))

  if (NROW(all_rates) == 0L) return(data.table::data.table())

  checkmate::assert_numeric(all_rates$rate, .var.name = "value returned by 'rate_fn'")
  if (anyDuplicated(all_rates, by = c(rate_by, "period")) > 0L) {
    stop("'rate_fn' must return a single growth rate per well and time window.")
  }

  # --- Stage 1b: aggregate replicates ---
  # Include all partner drug/conc slots so each combo arm aggregates separately
  agg_by <- intersect(
    c(cl_col, drug_col, conc_col, partner_drug_cols, partner_conc_cols, "period"),
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
  # untreated_tag = c("vehicle", "untreated") by default — the canonical list of
  # control drug name values in gDRutils; no need to hardcode "DMSO"/"vehicle".
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  is_primary_ctrl <- av_rates[[drug_col]] %in% untreated_tag |
                     av_rates[[conc_col]] == 0
  # Exclude rows that have any real partner drug in any slot (doublet or triplet)
  has_partner <- if (length(partner_drug_cols) == 0L) {
    rep(FALSE, NROW(av_rates))
  } else {
    Reduce("|", lapply(partner_drug_cols, function(dc) {
      !is.na(av_rates[[dc]]) & !av_rates[[dc]] %in% untreated_tag
    }))
  }
  is_dmso <- is_primary_ctrl & !has_partner

  # One baseline per (CellLine, period): combo layouts hold several control rows
  # per cell line (different partner slots/concentrations), which would turn the
  # merge below into a many-to-many join and duplicate NormalizedGrowthRate.
  dmso_baselines <- av_rates[is_dmso,
    .(rate_0 = mean(GrowthRate, na.rm = TRUE)),
    by = c(cl_col, "period")
  ]

  if (NROW(dmso_baselines) == 0L) {
    warning(".compute_growth_rates: no DMSO/vehicle rows found; NormalizedGrowthRate will be NA.")
    av_rates[, NormalizedGrowthRate := NA_real_]
    av_rates[, normalization_type := "NGR"]
    return(av_rates)
  }

  # Apply normalization_map.
  # "None" → NormalizedGrowthRate = GrowthRate (raw doublings/day, no DMSO division).
  # Any other value (e.g. "early") → NormalizedGrowthRate = GrowthRate / DMSO_rate[ref_period].
  norm_dt <- data.table::data.table(
    period = names(normalization_map),
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
        "period", "norm_period"
      ),
      by.x = c(cl_col, "norm_target"),
      by.y = c(cl_col, "norm_period"),
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
#' @return \code{SummarizedExperiment}.
#' @keywords internal
#' @noRd
.growth_dt_to_se <- function(growth_dt) {

  cl_col <- gDRutils::get_env_identifiers("cellline_name")
  drug_col <- gDRutils::get_env_identifiers("drug_name")
  conc_col <- gDRutils::get_env_identifiers("concentration")

  # Resolve all partner drug/conc slots (2, 3, ...) that exist in the data
  partner_cols <- .resolve_partner_cols(growth_dt)
  partner_drug_cols <- partner_cols$drug
  partner_conc_cols <- partner_cols$conc

  # Exclude control (vehicle/untreated) and zero-concentration rows from fitting
  conc_gt0 <- growth_dt[[conc_col]] > 0
  drug_not_ctrl <- !growth_dt[[drug_col]] %in%
    gDRutils::get_env_identifiers("untreated_tag")
  fit_dt <- growth_dt[conc_gt0 & drug_not_ctrl]

  if (NROW(fit_dt) == 0L) {
    stop(".growth_dt_to_se: no treated rows with Concentration > 0 remain after filtering.")
  }

  # row key: DrugName | all_partner_drugs | all_partner_concs | CellLineName
  # Includes all slots (2, 3, ...) so that:
  #   SA            (DrugName2=vehicle, DrugName3=vehicle)      → own row
  #   doublet        (DrugName2=DrugB,   DrugName3=vehicle)      → own row
  #   triplet        (DrugName2=DrugB,   DrugName3=DrugC)        → own row
  #   doublet @ 0.4µM vs 0.8µM partner conc                     → two separate rows
  key_parts <- intersect(
    c(drug_col, partner_drug_cols, partner_conc_cols, cl_col),
    names(fit_dt)
  )

  fit_dt[, row := do.call(paste, c(
    lapply(key_parts, function(k) {
      v <- fit_dt[[k]]
      if (is.numeric(v)) sprintf("%.4g", v) else as.character(v)
    }),
    sep = "|"
  ))]
  fit_dt[, column := period]

  # Assay columns: primary + all partner concentrations + response + label
  assay_cols <- intersect(
    c(conc_col, partner_conc_cols, "NormalizedGrowthRate", "normalization_type"),
    names(fit_dt)
  )
  assay_dt <- fit_dt[, assay_cols, with = FALSE]

  bm <- BumpyMatrix::splitAsBumpyMatrix(
    assay_dt,
    row = fit_dt$row,
    col = fit_dt$column
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(GrowthRates = bm)
  )
}
