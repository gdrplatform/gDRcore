# annotate_dt_with_cell_line

Annotate cell line data with the provided annotation table

## Usage

``` r
annotate_dt_with_cell_line(data, cell_line_annotation, fill = "unknown")
```

## Arguments

- data:

  data.table with dose-response data

- cell_line_annotation:

  data.table with cell line annotations

- fill:

  string indicating how unknown cell lines should be filled in the DB

## Value

data.table with annotated cell lines

## Examples

``` r
data <- data.table::data.table(
   clid = c("CL1", "CL2", "CL3"),
   Gnumber = c("D1", "D2", "D3")
)
cell_line_annotation <- get_cell_line_annotation(data)
annotated_metadata <- annotate_dt_with_cell_line(data, cell_line_annotation)
```
