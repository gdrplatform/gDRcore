# musyc_available

Whether
[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md)
can run: `reticulate` is installed and the configured Python has the
`synergy` package. The check never provisions a Python of its own, which
recent `reticulate` versions do when no Python is configured.

## Usage

``` r
musyc_available()
```

## Value

Logical scalar.

## See also

[`musyc_fit_fn`](https://gdrplatform.github.io/gDRcore/reference/musyc_fit_fn.md)

## Examples

``` r
musyc_available()
#> [1] FALSE
```
