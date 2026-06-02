# annotate_mae_with_drug

Annotate MultiAssayExperiment object with drug annotations

## Usage

``` r
annotate_mae_with_drug(mae, drug_annotation, fill = "unknown")
```

## Arguments

- mae:

  MultiAssayExperiment object containing dose-response data

- drug_annotation:

  data.table with drug annotations

- fill:

  string indicating how unknown drugs should be filled in the DB

## Value

MultiAssayExperiment object with annotated drugs

## Examples

``` r
mae <- MultiAssayExperiment::MultiAssayExperiment(
   experiments = list(exp1 = SummarizedExperiment::SummarizedExperiment(
     rowData = data.table::data.table(Gnumber = c("D1", "D2", "D3"))
   ))
)
drug_annotation <- get_drug_annotation(data.table::as.data.table(
   SummarizedExperiment::rowData(
     MultiAssayExperiment::experiments(mae)[[1]])))
annotated_mae <- annotate_mae_with_drug(mae, drug_annotation)
```
