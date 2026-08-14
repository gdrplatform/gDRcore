# apply_combo_isobolograms

Compute Loewe combination index isobolograms for a combination SE,
writing results to the `isobolograms` and `all_iso_points` assays. This
is **step 4** of the
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
pipeline.

## Usage

``` r
apply_combo_isobolograms(
  se,
  series_identifiers = NULL,
  normalization_types = c("GR", "RV"),
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  isobolograms_assay = "isobolograms",
  iso_points_assay = "all_iso_points"
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with `Averaged` and `Metrics` assays.

- series_identifiers:

  character vector of length 2. Defaults to
  `get_default_nested_identifiers(se, ...)`.

- normalization_types:

  character vector. Default `c("GR", "RV")`.

- averaged_assay:

  string. Default `"Averaged"`.

- metrics_assay:

  string. Default `"Metrics"`.

- isobolograms_assay:

  string; name of the Loewe curves output assay. Default
  `"isobolograms"`.

- iso_points_assay:

  string; name of the iso points output assay. Default
  `"all_iso_points"`.

## Value

Updated SE with `isobolograms_assay` and `iso_points_assay` assays.

## Details

Requires the `Averaged` and `Metrics` assays (as produced by
[`apply_combo_sa_fits`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_sa_fits.md)
or `fit_SE.combinations`). Isobologram analysis is only performed when
at least 9 combo concentration points are available — otherwise `NA` is
written, matching `fit_SE.combinations` behaviour.

## See also

[`apply_combo_sa_fits`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_sa_fits.md),
[`apply_combo_excess`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_excess.md),
[`apply_combo_scores`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_scores.md),
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
combo_se_iso <- apply_combo_isobolograms(combo_se[1, 1])
"isobolograms" %in% SummarizedExperiment::assayNames(combo_se_iso)
#> [1] TRUE
```
