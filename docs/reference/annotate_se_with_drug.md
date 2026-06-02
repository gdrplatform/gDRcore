# annotate_se_with_drug

Annotate SummarizedExperiment object with drug annotations

## Usage

``` r
annotate_se_with_drug(se, drug_annotation, fill = "unknown")
```

## Arguments

- se:

  SummarizedExperiment object containing dose-response data

- drug_annotation:

  data.table with drug annotations

- fill:

  string indicating how unknown drugs should be filled in the DB

## Value

SummarizedExperiment object with annotated drugs

## Examples

``` r
se <- SummarizedExperiment::SummarizedExperiment(
   rowData = data.table::data.table(Gnumber = c("D1", "D2", "D3"))
)
drug_annotation <- get_drug_annotation(data.table::as.data.table(SummarizedExperiment::rowData(se)))
annotated_se <- annotate_se_with_drug(se, drug_annotation)
```
