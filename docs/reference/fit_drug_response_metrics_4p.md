# fit_drug_response_metrics_4p

Variant of
[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics.md)
using a 4-parameter Hill model — `x_0` is fitted freely rather than
fixed at 1. Useful when the upper asymptote is expected to deviate from
1 (e.g. partial agonists, non-standard assay readouts).

## Usage

``` r
fit_drug_response_metrics_4p(
  avg_dt,
  x_col = "x",
  capping_fold = 5,
  range_conc = c(0.005, 5),
  pcutoff = 0.05,
  n_point_cutoff = 4L,
  force_fit = FALSE,
  cap = 0.1
)
```

## Arguments

- avg_dt:

  `data.table` of averaged data for one (drug x cell line x
  normalization_type) triplet. Should contain columns `Concentration`,
  the response column named by `x_col`, `normalization_type`, and
  optionally `x_std` (used for `x_sd_avg`).

- x_col:

  string; name of the response column in `avg_dt`. Defaults to `"x"`
  (the standard `Averaged` assay column name). Set to e.g.
  `"NormalizedGrowthRate"` when fitting time-course growth-rate data.

- capping_fold:

  numeric capping fold passed to
  [`cap_xc50`](https://gdrplatform.github.io/gDRstyle/reference/cap_xc50.html)

- range_conc:

  numeric vector of length 2 specifying the concentration range used to
  compute `x_AOC_range`. Default `c(5e-3, 5)`, matching the
  [`fit_SE()`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)
  default.

- pcutoff:

  numeric p-value cutoff for the fit. Default `0.05`.

- n_point_cutoff:

  integer minimum number of concentration points required to attempt
  fitting. Default `4L`.

- force_fit:

  logical; if `TRUE`, force fitting even when the number of points is
  below `n_point_cutoff`. Default `FALSE`.

- cap:

  numeric; upper cap on `xc50` relative to the concentration range.
  Default `0.1`.

## Value

Named list of fit metrics. `fit_type` will be `"DRC4pHillFitModel"` on
success.

## See also

[`fit_drug_response_metrics`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics.md)
for the default 3-parameter variant that matches
[`fit_SE()`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)
output.

## Examples

``` r
dt <- data.table::data.table(
  Concentration = c(0.001, 0.01, 0.1, 1, 10),
  x = c(0.95, 0.8, 0.5, 0.2, 0.1),
  normalization_type = "RV"
)
fit_drug_response_metrics_4p(dt)
#> $normalization_type
#> [1] "RV"
#> 
#> $x_mean
#> [1] 0.5069374
#> 
#> $x_AOC
#> [1] 0.4930626
#> 
#> $x_AOC_range
#> [1] 0.5470636
#> 
#> $x_max
#> [1] 0.1
#> 
#> $x_sd_avg
#> [1] NaN
#> 
#> $N_conc
#> [1] 5
#> 
#> $maxlog10Concentration
#> [1] 1
#> 
#> $ec50
#> [1] 0.07987809
#> 
#> $xc50
#> [1] 0.09549121
#> 
#> $h
#> h:(Intercept) 
#>      0.659329 
#> 
#> $r2
#> [1] 0.9996739
#> 
#> $rss
#> [1] 0.0001767251
#> 
#> $p_value
#> [1] 0.0004890517
#> 
#> $x_0
#> [1] 0.996396
#> 
#> $x_inf
#> [1] 0.05872765
#> 
#> $fit_type
#> [1] "DRC4pHillFitModel"
#> 
```
