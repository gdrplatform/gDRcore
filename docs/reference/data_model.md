# Detect model of data

Detect model of data

## Usage

``` r
data_model(x)
```

## Arguments

- x:

  data.table with raw data or SummarizedExperiment object with gDR
  assays

## Value

string with the information of the raw data follows single-agent or
combination data model

## Examples

``` r
data_model("single-agent")
#> [1] "single-agent"
```
