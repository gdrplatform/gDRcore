# Compute per-period growth rates for time-course screens

Stage 1 of the time-course pipeline: fits a growth rate per well within
each time window, aggregates replicates, and normalises the rates to the
control wells of a reference period. The returned table keeps every
treatment, controls included, together with the raw rate and its
standard deviation, so it can be reported on directly. Pass it through
[`growth_rates_to_se`](https://gdrplatform.github.io/gDRcore/reference/growth_rates_to_se.md)
to obtain the stage 2 input, or use
[`fit_SE.timecourse`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.timecourse.md)
to run both stages at once.

## Usage

``` r
compute_growth_rates(
  se,
  periods,
  normalization_map,
  lfc_assay = NULL,
  rate_fn = NULL
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with a log-fold-change assay.

- periods:

  named list of numeric vectors of length 2; half-open time windows
  `c(start_hour, end_hour)`, see
  [`fit_SE.timecourse`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.timecourse.md).

- normalization_map:

  named character vector mapping every period to its reference period or
  to `"None"`, see
  [`fit_SE.timecourse`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.timecourse.md).

- lfc_assay:

  string or `NULL`; name of the log-fold-change assay. `NULL` uses the
  `input_assay` of the `"time-course"` fit profile (see
  [`get_fit_profile`](https://gdrplatform.github.io/gDRcore/reference/get_fit_profile.md)).

- rate_fn:

  function or `NULL`; growth-rate estimator. It receives a `data.table`
  holding the timepoints of one well within one window — including the
  duration column and the log-fold-change column — and must return a
  single numeric value: the growth rate in doublings per day. `NULL`
  uses the slope of `lm(LogFoldChange ~ Duration)`.

## Value

A `data.table` with one row per
`(CellLine, Drug, Concentration, partner slots, period)` — controls
included — holding `GrowthRate` (doublings/day), `sd_GrowthRate` across
replicate wells, `rate_0` (the control baseline used, `NA` for periods
mapped to `"None"` and when no control rows are present),
`NormalizedGrowthRate` and `normalization_type`. All of these columns
are always present, whatever the normalization map or the availability
of controls.

## Details

Growth rates are computed per
`(Barcode, WellRow, WellColumn, CellLine, Drug, Concentration, partner drug/concentration slots)`
group, so single-agent, doublet and triplet arms are never pooled.
Control wells are detected via
`gDRutils::get_env_identifiers("untreated_tag")` and a single baseline
per `(CellLine, period)` is used for normalization.

## See also

[`growth_rates_to_se`](https://gdrplatform.github.io/gDRcore/reference/growth_rates_to_se.md),
[`fit_SE.timecourse`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.timecourse.md),
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
periods <- list(early = c(0, 48), late = c(48, 96))
norm_map <- c(early = "None", late = "early")
growth_dt <- compute_growth_rates(se_tc, periods, norm_map)

# Stage 2 with any custom fit function:
se_fit <- apply_fit(
  growth_rates_to_se(growth_dt),
  fit_fn = my_fit_fn,
  data_type = "time-course-metrics",
  output_assay = "Metrics",
  merge = "replace",
  fit_source = "custom"
)
} # }
```
