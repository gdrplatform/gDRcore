# hss_fit_fn

Reference fit function for Highest Single Agent (HSA) synergy scoring on
combination data. Uses raw averaged single-agent edges.

## Usage

``` r
hss_fit_fn(avg_dt)
```

## Arguments

- avg_dt:

  `data.table` for one (drug1 x drug2 x cell line x normalization_type)
  combination. Must contain columns `Concentration`, `Concentration_2`,
  `x`, and `normalization_type`.

## Value

Named list with `hss_score`, `hss_excess_mean`, `n_combo_points`, and
`normalization_type`.

## Details

Intended for use with
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
on combination SEs:

      apply_fit(
        combo_se, hss_fit_fn, "combination",
        output_assay = "custom_hss", fit_source = "hss"
      )

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
combo_se_out <- apply_fit(
  combo_se, hss_fit_fn, "combination",
  output_assay = "custom_hss", fit_source = "hss"
)
```
