# predict_musyc_surface

Response of a fitted MuSyC surface (`gamma12 = gamma21 = 1`) at the
given concentrations, e.g. to compare a
[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md)
fit with the measured grid.

## Usage

``` r
predict_musyc_surface(conc, conc_2, params)
```

## Arguments

- conc, conc_2:

  numeric vectors of the same length; concentrations of drug 1 and drug
  2.

- params:

  named list or one-row `data.table` with `E0`, `E1`, `E2`, `E3`, `h1`,
  `h2`, `C1`, `C2`, `alpha12` and `alpha21`, as returned by
  `musyc_fit_fn`.

## Value

Numeric vector of responses, one per concentration pair.

## See also

[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md)

## Examples

``` r
params <- list(E0 = 1, E1 = 0.4, E2 = 0.3, E3 = -0.2, h1 = 1.2, h2 = 0.9,
               C1 = 0.5, C2 = 0.2, alpha12 = 3, alpha21 = 0.5)
predict_musyc_surface(c(0, 0.1, 1), c(0, 1, 1), params)
#> [1] 1.0000000 0.3891412 0.1064444
```
