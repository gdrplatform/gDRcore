# fit_SE for time-course (Incucyte) screens

Two-stage pipeline for time-course dose-response data:

1.  Compute growth rates per user-defined time window and normalise them
    to the control of the reference period
    ([`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md)).

2.  Fit a 3-parameter Hill (logistic) curve on the
    `NormalizedGrowthRate ~ log10(Concentration)` dose-response for each
    `(Drug, CellLine, period)` combination.

## Usage

``` r
fit_SE.timecourse(
  se,
  periods,
  normalization_map,
  lfc_assay = NULL,
  metrics_assay = "Metrics",
  rate_fn = NULL,
  fit_fn = NULL,
  fit_source = NULL,
  n_point_cutoff = 4L,
  range_conc = c(0.005, 5),
  force_fit = FALSE,
  pcutoff = 0.05,
  cap = 0.1
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with a `LogFoldChange` assay (output of `normalize_SE` for
  `data_type = "time-course"`).

- periods:

  named list of numeric vectors of length 2; each element defines a
  half-open time window `c(start_hour, end_hour)` used for growth-rate
  calculation: a timepoint is used when
  `start_hour <= Duration < end_hour`, so contiguous windows partition
  the time-course without dropping the shared boundary. At minimum
  `periods$early` must be provided. Example:
  `list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))`.

- normalization_map:

  named character vector mapping *every* period name to the reference
  period for DMSO normalization, or `"None"` to keep the raw growth
  rate. Its names must match `names(periods)` exactly and its values
  must be `"None"` or an existing period name. Example:
  `c(early = "None", mid = "early", late = "early")`.

- lfc_assay:

  string; name of the log-fold-change assay in `se`. Defaults to the
  `input_assay` of the `"time-course"` fit profile.

- metrics_assay:

  string; name of the output metrics assay. Default `"Metrics"`.

- rate_fn:

  function or `NULL`; growth-rate estimator passed to
  [`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md).
  `NULL` uses the default `lm(LogFoldChange ~ Duration)` slope.

- fit_fn:

  function or `NULL`; dose-response fit applied to each
  `normalization_type` slice of the `GrowthRates` assay. `NULL` uses
  [`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics.md)
  on `NormalizedGrowthRate`. When a custom function is supplied,
  `n_point_cutoff`, `range_conc`, `force_fit`, `pcutoff` and `cap` are
  not used.

- fit_source:

  string or `NULL`; value recorded in the `fit_source` column. `NULL`
  records `"gDR"` for the default fit and `"custom"` for a user-supplied
  `fit_fn`.

- n_point_cutoff:

  integer; minimum number of unique concentrations required to attempt a
  curve fit. Default `4L`.

- range_conc:

  numeric vector of length 2; concentration range for `x_AOC_range`.
  Default `c(5e-3, 5)`.

- force_fit:

  logical; if `TRUE` fit is kept even if not statistically significant.
  Default `FALSE`.

- pcutoff:

  numeric; p-value threshold for the F-test. Default `0.05`.

- cap:

  numeric; upper-capping fold above `x_0 = 1`. Default `0.1`.

## Value

A new `SummarizedExperiment` where rows = (Drug, CellLine) and columns =
period names, holding the `GrowthRates` assay from stage 1 and the
`Metrics` assay (or `metrics_assay`) from stage 2. With the default
`fit_fn` the metrics are: `ec50`, `xc50`, `h`, `x_inf`, `x_0`, `r2`,
`x_mean`, `x_AOC`, `x_AOC_range`, `x_max`, `fit_type`, `N_conc`,
`maxlog10Concentration`, and `period`.

## Details

The `LogFoldChange` assay consumed by this function is produced by
[`normalize_SE`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)
when `data_type = "time-course"`.

Both stages accept a replacement function, so custom fitting works for
time-course data the same way it does for single-agent and combination
screens. The two stages can also be run separately:
[`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md)
returns the growth rate table (controls included) and
[`growth_rates_to_se`](https://gdrplatform.github.io/gDRcore/reference/growth_rates_to_se.md)
turns it into the input of
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
with `data_type = "time-course-metrics"`.

## See also

[`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md)
for stage 1 on its own,
[`fit_SE`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)
for single-agent fitting,
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
for combination screens,
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
for the generic fitting engine.

## Examples

``` r
if (FALSE) { # \dontrun{
# After normalize_SE(se_tc, data_type = "time-course"):
periods <- list(early = c(0, 48), mid = c(48, 96), late = c(96, 144))
norm_map <- c(early = "None", mid = "early", late = "early")
se_fit <- fit_SE.timecourse(se_tc, periods = periods,
                            normalization_map = norm_map)
gDRutils::convert_se_assay_to_dt(se_fit, "Metrics")

# Custom growth-rate estimator and custom dose-response fit:
se_custom <- fit_SE.timecourse(
  se_tc,
  periods = periods,
  normalization_map = norm_map,
  rate_fn = function(dt) 24 * stats::median(diff(dt$LogFoldChange) /
                                            diff(dt$Duration)),
  fit_fn = function(dt) data.table::data.table(x_mean = mean(dt$NormalizedGrowthRate))
)
} # }
```
