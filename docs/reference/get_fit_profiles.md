# Get all registered fit profiles

Returns the list of all experiment-type profiles used by
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
and
[`apply_fits`](https://gdrplatform.github.io/gDRcore/reference/apply_fits.md)
to resolve default slicing configuration.

## Usage

``` r
get_fit_profiles()
```

## Value

Named list of profile definitions. Each element contains `slicing_cols`,
`slicing_values`, `input_assay`, and `description`.

## Examples

``` r
get_fit_profiles()
#> $`single-agent`
#> $`single-agent`$slicing_cols
#> [1] "normalization_type"
#> 
#> $`single-agent`$slicing_values
#> [1] "GR" "RV"
#> 
#> $`single-agent`$input_assay
#> [1] "Averaged"
#> 
#> $`single-agent`$nested_cols
#> [1] "concentration"
#> 
#> $`single-agent`$nested_cols_note
#> [1] "identifier keys resolved at runtime via gDRutils::get_env_identifiers() — actual column names may differ in custom environments"
#> 
#> $`single-agent`$description
#> [1] "1D dose-response curve; fit_fn receives one data.table per normalization_type slice, rows indexed by concentration"
#> 
#> 
#> $`time-course`
#> $`time-course`$slicing_cols
#> [1] "normalization_type"
#> 
#> $`time-course`$slicing_values
#> [1] "GR" "RV"
#> 
#> $`time-course`$input_assay
#> [1] "Averaged"
#> 
#> $`time-course`$nested_cols
#> [1] "concentration" "duration"     
#> 
#> $`time-course`$nested_cols_note
#> [1] "identifier keys resolved at runtime via gDRutils::get_env_identifiers() — actual column names may differ in custom environments"
#> 
#> $`time-course`$description
#> [1] "Time-series per dose; fit_fn receives one data.table per normalization_type slice, rows indexed by concentration x duration (time axis)"
#> 
#> 
#> $combination
#> $combination$slicing_cols
#> [1] "normalization_type"
#> 
#> $combination$slicing_values
#> [1] "GR" "RV"
#> 
#> $combination$input_assay
#> [1] "Averaged"
#> 
#> $combination$nested_cols
#> [1] "concentration"  "concentration2"
#> 
#> $combination$nested_cols_note
#> [1] "identifier keys resolved at runtime via gDRutils::get_env_identifiers() — actual column names may differ in custom environments"
#> 
#> $combination$description
#> [1] "2D concentration grid; fit_fn receives one data.table per normalization_type slice, rows indexed by concentration x concentration2"
#> 
#> 
```
