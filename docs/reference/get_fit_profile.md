# Get a single fit profile by name

Get a single fit profile by name

## Usage

``` r
get_fit_profile(name)
```

## Arguments

- name:

  string; profile name (e.g. `"single-agent"`).

## Value

Named list with the following fields:

- `slicing_cols`:

  Column(s) used to split each BumpyMatrix cell into sub-experiments
  before calling `fit_fn` (e.g. `"normalization_type"`).

- `slicing_values`:

  Values of `slicing_cols` to iterate; `NULL` means all unique values
  found in each cell.

- `input_assay`:

  Default source assay name.

- `nested_cols`:

  Informational: the column(s) that index rows *inside* each BumpyMatrix
  cell (i.e., what `fit_fn` receives as row structure). Values are the
  canonical default identifiers resolved by
  [`gDRutils::get_env_identifiers()`](https://gdrplatform.github.io/gDRstyle/reference/identifiers.html)
  at package load — they may differ if identifiers are customised at
  runtime.

- `nested_cols_note`:

  Human-readable note about how `nested_cols` are resolved.

- `description`:

  Human-readable description of the profile.

Stops if the profile does not exist.

## Examples

``` r
get_fit_profile("single-agent")
#> $slicing_cols
#> [1] "normalization_type"
#> 
#> $slicing_values
#> [1] "GR" "RV"
#> 
#> $input_assay
#> [1] "Averaged"
#> 
#> $nested_cols
#> [1] "concentration"
#> 
#> $nested_cols_note
#> [1] "identifier keys resolved at runtime via gDRutils::get_env_identifiers() — actual column names may differ in custom environments"
#> 
#> $description
#> [1] "1D dose-response curve; fit_fn receives one data.table per normalization_type slice, rows indexed by concentration"
#> 
```
