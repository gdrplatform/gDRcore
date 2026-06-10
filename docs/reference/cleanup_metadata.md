# cleanup_metadata

Cleanup a data.table with metadata

## Usage

``` r
cleanup_metadata(
  df_metadata,
  cell_line_annotation = NULL,
  drug_annotation = NULL
)
```

## Arguments

- df_metadata:

  a data.table with metadata

- cell_line_annotation:

  optional data.table with cell line annotations; if NULL (default),
  annotations are looked up from gDRinternal or gDRtestData

- drug_annotation:

  optional data.table with drug annotations; if NULL (default),
  annotations are looked up from gDRinternal or gDRtestData

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
