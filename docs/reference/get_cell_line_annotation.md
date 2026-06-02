# get_cell_line_annotation

Get cell line annotation data table

## Usage

``` r
get_cell_line_annotation(
  data,
  fname = "cell_lines.csv",
  fill = "unknown",
  annotation_package = if ("gDRinternal" %in% .packages(all.available = TRUE)) {
    
    "gDRinternal"
 } else {
     "gDRtestData"
 }
)
```

## Arguments

- data:

  data.table with cell line identifiers to be matched

- fname:

  string with file name containing the annotation

- fill:

  string indicating how unknown cell lines should be filled in the DB

- annotation_package:

  string indicating name of the package containing cell line annotation

## Value

data.table with cell line annotations

## Examples

``` r
data <- data.table::data.table(clid = c("CL1", "CL2", "CL3"))
cell_line_annotation <- get_cell_line_annotation(data)
```
