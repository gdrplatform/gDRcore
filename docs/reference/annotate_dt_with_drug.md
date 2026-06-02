# annotate_dt_with_drug

Annotate drug data with the provided annotation table

## Usage

``` r
annotate_dt_with_drug(data, drug_annotation, fill = "unknown")
```

## Arguments

- data:

  data.table with dose-response data

- drug_annotation:

  data.table with drug annotations

- fill:

  string indicating how unknown drugs should be filled in the DB

## Value

data.table with annotated drugs

## Examples

``` r
data <- data.table::data.table(
   clid = c("CL1", "CL2", "CL3"),
   Gnumber = c("D1", "D2", "D3")
)
drug_annotation <- get_drug_annotation(data)
annotated_metadata <- annotate_dt_with_drug(data, drug_annotation)
```
