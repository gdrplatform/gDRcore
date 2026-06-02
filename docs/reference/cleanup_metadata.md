# cleanup_metadata

Cleanup a data.table with metadata

## Usage

``` r
cleanup_metadata(df_metadata)
```

## Arguments

- df_metadata:

  a data.table with metadata

## Value

a data.table with cleaned metadata

## Details

Adds annotations and check whether user provided correct input data.

## Examples

``` r
df <- data.table::data.table(
  clid = "CELL_LINE",
  Gnumber = "DRUG_1",
  Concentration = c(0, 1),
  Duration = 72
)
cleanup_df <- cleanup_metadata(df)
```
