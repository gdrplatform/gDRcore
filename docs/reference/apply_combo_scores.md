# apply_combo_scores

Compute Bliss and HSA synergy scores for a combination SE, replicating
the exact scoring logic of
[`fit_SE.combinations`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md).

## Usage

``` r
apply_combo_scores(
  se,
  scores_assay = "scores",
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  normalization_types = c("GR", "RV"),
  fit_source = "gDR"
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with `"Averaged"` and `"Metrics"` assays (combination experiment).

- scores_assay:

  string; name of the output assay to write scores into. Default
  `"scores"`.

- averaged_assay:

  string; name of the input assay. Default `"Averaged"`.

- metrics_assay:

  string; name of the assay containing SA fit parameters. Default
  `"Metrics"`.

- normalization_types:

  character vector of normalization types to process. Default
  `c("GR", "RV")`.

- fit_source:

  string recorded in the `fit_source` column. Default `"gDR"`.

## Value

Updated `SummarizedExperiment` with a `scores_assay` assay containing
`bliss_score` and `hsa_score` per triplet.

## Details

Unlike the low-level
[`bliss_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/bliss_fit_fn.md)
and
[`hss_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/hss_fit_fn.md)
(which work on raw Averaged data per triplet), this function uses the SA
fit parameters from the `Metrics` assay to generate curve-smoothed
single-agent responses via
[`predict_efficacy_from_conc`](https://gdrplatform.github.io/gDRstyle/reference/predict_efficacy_from_conc.html)
before computing excess. This produces results numerically identical to
`fit_SE.combinations`.

The function requires a `Metrics` assay containing columns
`dilution_drug`, `cotrt_value`, `ec50`, `h`, `x_inf`, `x_0`, and
`normalization_type` — as produced by `fit_SE.combinations` or
[`apply_fit_to_se`](https://gdrplatform.github.io/gDRcore/reference/apply_fit_to_se.md)
with `fit_drug_response_metrics`.

Scoring steps (per drug-combo × cell-line × normalization_type):

1.  Predict smooth SA responses at every combo concentration using
    `predict_efficacy_from_conc` on `drug_1` and `drug_2` parameter sets
    from `Metrics`.

2.  Average `col_values` (drug-1-along-conc1) and `row_values`
    (drug-2-along-conc2) → `smooth`.

3.  Compute Bliss-expected using
    [`calculate_Bliss`](https://gdrplatform.github.io/gDRcore/reference/calculate_matrix_metric.md)
    and HSA-expected using
    [`calculate_HSA`](https://gdrplatform.github.io/gDRcore/reference/calculate_matrix_metric.md)
    on the smooth SA edges.

4.  Compute excess via
    [`calculate_excess`](https://gdrplatform.github.io/gDRcore/reference/calculate_excess.md).

5.  Score = mean of top 10-percentile excess via
    [`calculate_score`](https://gdrplatform.github.io/gDRcore/reference/calculate_score.md).

## See also

[`bliss_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/bliss_fit_fn.md),
[`hss_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/hss_fit_fn.md)
for the simplified raw-data variants.

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
combo_se_out <- apply_combo_scores(combo_se)
#> Loading required namespace: BumpyMatrix
"scores" %in% SummarizedExperiment::assayNames(combo_se_out)
#> [1] TRUE
```
