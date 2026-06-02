# Standardize concentrations.

Utilize a map to standardize concentrations.

## Usage

``` r
replace_conc_with_standardized_conc(
  original_concs,
  conc_map,
  original_conc_col,
  standardized_conc_col
)
```

## Arguments

- original_concs:

  numeric vector of concentrations to replace using `conc_map`.

- conc_map:

  data.table of two columns named `original_conc_col` and
  `standardized_conc_col`.

- original_conc_col:

  string of the name of the column in `conc_map` containing the original
  concentrations to replace.

- standardized_conc_col:

  string of the name of the column in `conc_map` containing the
  standardized concentrations to use for replacement.

## Value

numeric vector of standardized concentrations.

## See also

map_conc_to_standardized_conc

## Examples

``` r
conc_map <- data.table::data.table(
  orig = c(0.99, 0.6, 0.456, 0.4),
  std = c(1, 0.6, 0.46, 0.4)
)
original_concs <- c(0.456, 0.456, 0.4, 0.99)
exp <- c(0.46, 0.46, 0.4, 1)
obs <- replace_conc_with_standardized_conc(
  original_concs,
  conc_map,
  original_conc_col = "orig",
  standardized_conc_col = "std"
)
```
