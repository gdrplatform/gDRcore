# apply_fit_to_se

Apply a user-supplied fit function per (drug x cell line x
normalization_type) triplet and persist results into the Metrics assay.

## Usage

``` r
apply_fit_to_se(
  se,
  fit_fn,
  normalization_types = c("GR", "RV"),
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  merge = "merge",
  on_error = "warn",
  fit_source = "custom"
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with an Averaged assay

- fit_fn:

  function taking a `data.table` and returning a named list

- normalization_types:

  character vector of types to iterate

- averaged_assay:

  name of the input assay

- metrics_assay:

  name of the output assay

- merge:

  `"merge"` (idempotent upsert) or `"replace"`

- on_error:

  `"warn"` (skip + warning) or `"stop"`

- fit_source:

  string recorded in the `fit_source` column

## Value

updated `SummarizedExperiment`

## Details

This is a convenience wrapper around
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
for the single-agent use case with the standard Metrics assay.

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
se <- mae[["single-agent"]]
simple_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
se_out <- apply_fit_to_se(se, simple_fn, fit_source = "demo")
#> Warning: fit_fn output missing columns: x_mean_sd, x_AOC_sd, x_AOC_range_sd, xc50_sd, x_max_sd, ec50_sd, x_inf_sd, x_0_sd, h_sd, r2_sd, x_sd_avg_sd
```
