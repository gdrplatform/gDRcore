# annotate_se_with_cell_line

Annotate SummarizedExperiment object with cell line annotations

## Usage

``` r
annotate_se_with_cell_line(se, cell_line_annotation, fill = "unknown")
```

## Arguments

- se:

  SummarizedExperiment object containing dose-response data

- cell_line_annotation:

  data.table with cell line annotations

- fill:

  string indicating how unknown cell lines should be filled in the DB

## Value

SummarizedExperiment object with annotated cell lines

## Examples

``` r
se <- SummarizedExperiment::SummarizedExperiment(
   rowData = data.table::data.table(clid = c("CL1", "CL2", "CL3"))
)
cell_line_annotation <- get_cell_line_annotation(data.table::as.data.table(SummarizedExperiment::rowData(se)))
annotated_se <- annotate_se_with_cell_line(se, cell_line_annotation)
```
