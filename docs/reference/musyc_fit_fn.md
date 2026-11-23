# musyc_fit_fn

Fit function for the MuSyC (multi-dimensional synergy of combinations)
model on combination dose-response data, for use with
[`apply_fit`](https://gdrplatform.github.io/gDRstyle/reference/apply_fit.html).

## Usage

``` r
musyc_fit_fn(avg_dt, n_restarts = 10L, seed = 1L)
```

## Arguments

- avg_dt:

  `data.table` for one (drug1 x drug2 x cell line x normalization_type)
  combination. Must contain columns `Concentration`, `Concentration_2`,
  `x`, and `normalization_type`.

- n_restarts:

  integer; number of restarts from random initial values.

- seed:

  integer; seed for the random initial values. The global random number
  generator is left as it was.

## Value

Named list with `normalization_type`, the parameters `E0`, `E1`, `E2`,
`E3`, `h1`, `h2`, `C1`, `C2`, `alpha12`, `alpha21`, the efficacy synergy
`beta`, `r2` and `converged`. Numeric fields are `NA_real_` when no fit
succeeds.

## Details

The surface (with `gamma12 = gamma21 = 1`) is fitted in R by bounded
least squares, [`nlminb`](https://rdrr.io/r/stats/nlminb.html) with the
analytic Jacobian; the settings follow the Python `synergy` package,
which an earlier version of this function called.

The whole grid is fitted, single-agent edges included, with three
additions to the `synergy` defaults: the untreated corner is set to 1
(added when absent), `E` is bounded from -1 for GR and from 0 for RV up
to 1.5, and `C1`/`C2` stay within a factor of 1000 of the tested doses.
`h`, `C` and `alpha` are fitted on a log scale. The fit starts from the
`synergy` initial guess (Hill fits of the two single-agent edges) and
from `n_restarts` random initial values; the fit with the best `r2`,
computed on the measured points, is kept, preferring fits that converged
within the iteration limit.

`alpha > 1` means synergistic potency, `beta > 0` synergistic efficacy.
`alpha` is weakly constrained when a drug has little effect, so extreme
values of it are not a measure of synergy.

## See also

[`musyc_summary_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_summary_fn.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
combo_se_out <- gDRutils::apply_fit(
  combo_se, musyc_fit_fn, "combination",
  output_assay = "musyc_params",
  summary_fn = musyc_summary_fn, summary_assay = "musyc_summary",
  fit_source = "musyc"
)
```
