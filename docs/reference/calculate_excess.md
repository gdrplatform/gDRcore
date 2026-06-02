# Calculate the difference between values in two data.tables

Calculate the difference between values, likely representing the same
metric, from two data.tables.

## Usage

``` r
calculate_excess(
  metric,
  measured,
  series_identifiers,
  metric_col,
  measured_col
)
```

## Arguments

- metric:

  data.table often representing readouts derived by calculating some
  metric. Examples of this could include hsa or bliss calculations from
  single-agent data.

- measured:

  data.table often representing measured data from an experiment.

- series_identifiers:

  character vector of identifiers in `measured` or `metric` which define
  a unique data point.

- metric_col:

  string of the column in `metric` to use in excess calculation.

- measured_col:

  string of the column in `measured` to use in excess calculation.

## Value

data.table of `measured`, now with an additional column named `excess`
(positive values for synergy/benefit).

## Examples

``` r
metric <- data.table::data.table(
  Concentration = c(1, 2, 3, 1, 2, 3),
  Concentration_2 = c(1, 1, 1, 2, 2, 2),
  GRvalue = c(100, 200, 300, 400, 500, 600)
)
measured <- data.table::data.table(
  Concentration = c(3, 1, 2, 2, 1, 3),
  Concentration_2 = c(1, 1, 1, 2, 2, 2),
  testvalue = c(200, 0, 100, 400, 300, 500)
)
series_identifiers <- c("Concentration", "Concentration_2")
metric_col <- "GRvalue"
measured_col <- "testvalue"
calculate_excess(
  metric,
  measured,
  series_identifiers,
  metric_col,
  measured_col
)
#>    Concentration Concentration_2     x
#>            <num>           <num> <num>
#> 1:             3               1   100
#> 2:             1               1   100
#> 3:             2               1   100
#> 4:             2               2   100
#> 5:             1               2   100
#> 6:             3               2   100
```
