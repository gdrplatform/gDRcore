# Convert a growth rate table into the stage 2 SummarizedExperiment

Builds the `GrowthRates` SummarizedExperiment consumed by
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
with `data_type = "time-course-metrics"`: rows are
`(Drug, partner drug/concentration slots, CellLine)` — so every
combination arm is fitted as its own dose-response curve — and columns
are the period names. Control and zero-concentration rows are dropped,
since they carry no dose-response information.

## Usage

``` r
growth_rates_to_se(growth_dt)
```

## Arguments

- growth_dt:

  `data.table`; output of
  [`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md).

## Value

A `SummarizedExperiment` with a `GrowthRates` assay holding
`Concentration`, the partner concentrations, `NormalizedGrowthRate` and
`normalization_type` per cell.

## See also

[`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md),
[`fit_SE.timecourse`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.timecourse.md)

## Examples

``` r
if (FALSE) { # \dontrun{
growth_dt <- compute_growth_rates(se_tc, periods, norm_map)
se_gr <- growth_rates_to_se(growth_dt)
} # }
```
