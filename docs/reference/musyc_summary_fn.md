# musyc_summary_fn

Summary function for
[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md):
called once per (drug1 x drug2 x cell line) cell on the fits of all
normalization types.

## Usage

``` r
musyc_summary_fn(fit_dt)
```

## Arguments

- fit_dt:

  `data.table` of the `musyc_fit_fn` results for one cell.

## Value

Named list with `mean_beta`, `mean_log2_alpha` (mean of `log2(alpha12)`
and `log2(alpha21)`, each clipped to `[-10, 10]`), `mean_r2`,
`any_converged` and `synergy_call` (`"synergistic"`, `"antagonistic"`,
`"additive"` or `"unknown"`, from `mean_beta` with a dead zone of 0.1).
`synergy_call` reflects efficacy only: a pair with strong potency
synergy (`alpha > 1`) and `beta` near 0 is called `"additive"`, so read
it together with `mean_log2_alpha`.

## See also

[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md)

## Examples

``` r
fit_dt <- data.table::data.table(
  normalization_type = c("GR", "RV"), alpha12 = c(2, 1.5),
  alpha21 = c(1.2, 1), beta = c(0.3, 0.2), r2 = c(0.9, 0.95),
  converged = TRUE
)
musyc_summary_fn(fit_dt)
#> $mean_beta
#> [1] 0.25
#> 
#> $mean_log2_alpha
#> [1] 0.4619992
#> 
#> $mean_r2
#> [1] 0.925
#> 
#> $any_converged
#> [1] TRUE
#> 
#> $synergy_call
#> [1] "synergistic"
#> 
```
