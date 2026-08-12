# Register or update a fit profile

Adds a new experiment-type profile or overwrites an existing one. The
profile is available immediately in the current R session and is picked
up by
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
and
[`apply_fits`](https://gdrplatform.github.io/gDRcore/reference/apply_fits.md).

## Usage

``` r
register_fit_profile(
  name,
  slicing_cols,
  slicing_values = NULL,
  input_assay,
  description = ""
)
```

## Arguments

- name:

  string; unique profile identifier (e.g. `"biochemical"`).

- slicing_cols:

  character vector of column names used to split each BumpyMatrix cell
  into sub-experiments.

- slicing_values:

  character vector of values to iterate, or `NULL` to iterate all unique
  values found in each cell.

- input_assay:

  string; name of the default source assay.

- description:

  string; human-readable description (optional).

## Value

Invisibly returns the registered profile list.

## Examples

``` r
register_fit_profile(
  "biochemical",
  slicing_cols   = "assay_type",
  slicing_values = c("Ki", "IC50"),
  input_assay    = "Biochemical",
  description    = "Biochemical activity assay (Ki / IC50)"
)
get_fit_profile("biochemical")
#> $slicing_cols
#> [1] "assay_type"
#> 
#> $slicing_values
#> [1] "Ki"   "IC50"
#> 
#> $input_assay
#> [1] "Biochemical"
#> 
#> $description
#> [1] "Biochemical activity assay (Ki / IC50)"
#> 
```
