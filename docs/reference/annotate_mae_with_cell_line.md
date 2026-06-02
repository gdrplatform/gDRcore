# annotate_mae_with_cell_line

Annotate MultiAssayExperiment object with cell line annotations

## Usage

``` r
annotate_mae_with_cell_line(mae, cell_line_annotation, fill = "unknown")
```

## Arguments

- mae:

  MultiAssayExperiment object containing dose-response data

- cell_line_annotation:

  data.table with cell line annotations

- fill:

  string indicating how unknown cell lines should be filled in the DB

## Value

MultiAssayExperiment object with annotated cell lines

## Examples

``` r
mae <- MultiAssayExperiment::MultiAssayExperiment(
   experiments = list(exp1 = SummarizedExperiment::SummarizedExperiment(
     rowData = data.table::data.table(clid = c("CL1", "CL2", "CL3"))
   ))
)
cell_line_annotation <- get_cell_line_annotation(data.table::as.data.table(
   SummarizedExperiment::rowData(
     MultiAssayExperiment::experiments(mae)[[1]])))
annotated_mae <- annotate_mae_with_cell_line(mae, cell_line_annotation)
```
