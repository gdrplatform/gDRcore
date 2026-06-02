# Retrieve the drug annotation from the annotated dt input

Retrieve the drug annotation from the annotated dt input

## Usage

``` r
get_drug_annotation_from_dt(dt)
```

## Arguments

- dt:

  annotated data.table

## Value

data.table with drug annotation

## Examples

``` r
dt <- data.table::data.table(Gnumber = "A",
DrugName = "drugA",
drug_moa = "drug_moa_A")
get_drug_annotation_from_dt(dt)
#>    Gnumber DrugName   drug_moa
#>     <char>   <char>     <char>
#> 1:       A    drugA drug_moa_A
```
