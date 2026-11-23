# musyc_available

Deprecated. Whether
[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md)
can run, which is always the case now that the fit no longer needs
Python. Will be removed in the next release.

## Usage

``` r
musyc_available()
```

## Value

Logical scalar, `TRUE`.

## See also

[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md)

## Examples

``` r
suppressWarnings(musyc_available())
#> [1] TRUE
```
