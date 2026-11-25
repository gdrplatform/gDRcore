# apply_combo_sa_fits

Fit single-agent dose-response curves for each co-treatment
concentration in a combination SE and store the results in the `Metrics`
assay. This is **step 1+2** of the
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
pipeline:

## Usage

``` r
apply_combo_sa_fits(
  se,
  series_identifiers = NULL,
  normalization_types = c("GR", "RV"),
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  fit_source = "gDR",
  fit_args = list()
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with an `Averaged` assay (combination experiment).

- series_identifiers:

  character vector of length 2: column names for the two concentration
  axes (e.g. `c("Concentration", "Concentration_2")`). Defaults to
  `get_default_nested_identifiers(se, data_model("combination"))`.

- normalization_types:

  character vector of normalization types to process. Default
  `c("GR", "RV")`.

- averaged_assay:

  string; name of the input assay. Default `"Averaged"`.

- metrics_assay:

  string; name of the output assay. Default `"Metrics"`.

- fit_source:

  string recorded in the `fit_source` column. Default `"gDR"`.

- fit_args:

  named list of arguments for the curve fits, passed to
  [`gDRutils::fit_drug_response_metrics()`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.html)
  (e.g. `solver`). Each series keeps its own `x_0`, `cap` and
  `force_fit`. Default [`list()`](https://rdrr.io/r/base/list.html).

## Value

Updated `SummarizedExperiment` with a `metrics_assay` assay containing
SA fit parameters per co-treatment concentration.

## Details

1.  Standardise concentrations and build the complete dose-response
    matrix.

2.  Fit SA curves per co-treatment concentration via
    `fit_combo_cotreatments()` and `fit_combo_codilutions()`.

3.  Compute **smooth** SA predictions at every combo point via
    [`map_ids_to_fits()`](https://gdrplatform.github.io/gDRcore/reference/map_ids_to_fits.md)
    and store them together with the SA fit parameters.

The resulting `Metrics` assay contains `dilution_drug`, `cotrt_value`,
`ec50`, `h`, `x_inf`, `x_0`, and `normalization_type` — the same schema
as produced by
[`fit_SE.combinations()`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md).

## See also

[`apply_combo_excess`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_excess.md),
[`apply_combo_isobolograms`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_isobolograms.md),
[`apply_combo_scores`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_scores.md),
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
SummarizedExperiment::assays(combo_se) <-
  SummarizedExperiment::assays(combo_se)["Averaged"]
combo_se_fitted <- apply_combo_sa_fits(combo_se[1, 1])
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
"Metrics" %in% SummarizedExperiment::assayNames(combo_se_fitted)
#> [1] TRUE
```
