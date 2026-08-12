# fit_drug_response_metrics

Reference fit function replicating standard gDR Hill curve fitting
(3-parameter model, `x_0` fixed at 1). Produces output column-compatible
with
[`fit_SE`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)
and
[`logisticFit`](https://gdrplatform.github.io/gDRstyle/reference/logisticFit.html),
including `p_value`, `rss`, `x_AOC_range`, `x_max`, and `x_sd_avg`. For
use with
[`apply_fit_to_se`](https://gdrplatform.github.io/gDRcore/reference/apply_fit_to_se.md)
or
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
on single-agent or time-course data.

## Usage

``` r
fit_drug_response_metrics(
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

Named list of fit metrics fully compatible with the standard gDR
`Metrics` assay schema: `ec50`, `xc50`, `h`, `x_inf`, `x_0`, `r2`,
`rss`, `p_value`, `x_mean`, `x_AOC`, `x_AOC_range`, `x_max`, `x_sd_avg`,
`N_conc`, `maxlog10Concentration`, `fit_type`, `normalization_type`.

## See also

[`fit_drug_response_metrics_4p`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics_4p.md)
for a 4-parameter variant that fits `x_0` freely.

## Examples

``` r
dt <- data.table::data.table(
  Concentration = c(0.001, 0.01, 0.1, 1, 10),
  x = c(0.95, 0.8, 0.5, 0.2, 0.1),
  normalization_type = "RV"
)
fit_drug_response_metrics(dt)
#> $normalization_type
#> [1] "RV"
#> 
#> $x_mean
#> [1] 0.5070152
#> 
#> $x_AOC
#> [1] 0.4929848
#> 
#> $x_AOC_range
#> [1] 0.5471268
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
#> [1] 0.07913255
#> 
#> $xc50
#> [1] 0.09530136
#> 
#> $h
#> h:(Intercept) 
#>     0.6523032 
#> 
#> $r2
#> [1] 0.9996676
#> 
#> $rss
#> [1] 0.00018015
#> 
#> $p_value
#> [1] 0.000006059718
#> 
#> $x_0
#> [1] 1
#> 
#> $x_inf
#> [1] 0.05710528
#> 
#> $fit_type
#> [1] "DRC3pHillFitModelFixS0"
#> 

# Custom response column (e.g. time-course growth rates):
dt2 <- data.table::data.table(
  Concentration = c(0.001, 0.01, 0.1, 1, 10),
  NormalizedGrowthRate = c(0.98, 0.82, 0.55, 0.21, 0.09),
  normalization_type = "GR"
)
fit_drug_response_metrics(dt2, x_col = "NormalizedGrowthRate")
#> $normalization_type
#> [1] "GR"
#> 
#> $x_mean
#> [1] 0.5299304
#> 
#> $x_AOC
#> [1] 0.4700696
#> 
#> $x_AOC_range
#> [1] 0.5198366
#> 
#> $x_max
#> [1] 0.09
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
#> [1] 0.111215
#> 
#> $xc50
#> [1] 0.1259556
#> 
#> $h
#> h:(Intercept) 
#>     0.6665989 
#> 
#> $r2
#> [1] 0.9981679
#> 
#> $rss
#> [1] 0.001068109
#> 
#> $p_value
#> [1] 0.00007841885
#> 
#> $x_0
#> [1] 1
#> 
#> $x_inf
#> [1] 0.03980943
#> 
#> $fit_type
#> [1] "DRC3pHillFitModelFixS0"
#> 
```
