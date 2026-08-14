# apply_combo_excess

Compute per-point HSA and Bliss excess values for a combination SE,
writing the result to the `excess` assay. This is **step 3** of the
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
pipeline.

## Usage

``` r
apply_combo_excess(
  se,
  series_identifiers = NULL,
  normalization_types = c("GR", "RV"),
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  excess_assay = "excess"
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with `Averaged` and `Metrics` assays.

- series_identifiers:

  character vector of length 2: concentration column names. Defaults to
  `get_default_nested_identifiers(se, ...)`.

- normalization_types:

  character vector. Default `c("GR", "RV")`.

- averaged_assay:

  string. Default `"Averaged"`.

- metrics_assay:

  string. Default `"Metrics"`.

- excess_assay:

  string; name of the output assay. Default `"excess"`.

## Value

Updated SE with an `excess_assay` assay containing per-point `smooth`,
`hsa_excess`, and `bliss_excess` columns.

## Details

Requires the `Averaged` and `Metrics` assays to be present (the latter
as produced by
[`apply_combo_sa_fits`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_sa_fits.md)
or `fit_SE.combinations`).

## See also

[`apply_combo_sa_fits`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_sa_fits.md),
[`apply_combo_isobolograms`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_isobolograms.md),
[`apply_combo_scores`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_scores.md),
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
combo_se_excess <- apply_combo_excess(combo_se[1, 1])
#> Loading required namespace: BumpyMatrix
"excess" %in% SummarizedExperiment::assayNames(combo_se_excess)
#> [1] TRUE
```
