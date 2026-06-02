# Map references

Map references

## Usage

``` r
.map_references(
  mat_elem,
  rowData_colnames = c(gDRutils::get_env_identifiers("duration"), paste0(c("drug",
    "drug_name", "drug_moa"), "3"))
)
```

## Arguments

- mat_elem:

  data.table input containing experimental metadata and row identifiers.

- rowData_colnames:

  character vector of variables (column names) used to identify and map
  reference treatments.

## Value

list
