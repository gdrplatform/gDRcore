# Get default nested identifiers

Get default nested identifiers

## Usage

``` r
get_default_nested_identifiers(x, data_model = NULL)

# S3 method for class 'data.table'
get_default_nested_identifiers(x, data_model = NULL)

# S3 method for class 'SummarizedExperiment'
get_default_nested_identifiers(x, data_model = NULL)
```

## Arguments

- x:

  data.table with raw data or `SummarizedExperiment` object with gDR
  assays

- data_model:

  single-agent vs combination

## Value

vector of nested identifiers

## Examples

``` r
get_default_nested_identifiers(data.table::data.table())
#> $`single-agent`
#> [1] "Concentration"
#> 
#> $combination
#> [1] "Concentration"   "Concentration_2"
#> 
#> $`time-course`
#> [1] "Concentration"   "Concentration_2"
#> 
```
