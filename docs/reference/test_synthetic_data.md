# Testing synthetic data form gDRtestData package

Testing synthetic data form gDRtestData package

## Usage

``` r
test_synthetic_data(
  original,
  data,
  dataName,
  override_untrt_controls = NULL,
  assays = c("Normalized", "Averaged", "Metrics"),
  tolerance = 0.001
)
```

## Arguments

- original:

  original MAE assay

- data:

  datase MAE or data.table

- dataName:

  dataset name

- override_untrt_controls:

  named list containing defining factors in the treatments

- assays:

  assays to test

- tolerance:

  tolerance factor

## Value

`NULL`

## Examples

``` r
set.seed(2)
cell_lines <- gDRtestData::create_synthetic_cell_lines()
drugs <- gDRtestData::create_synthetic_drugs()
data <- "finalMAE_small"
original <- gDRutils::get_synthetic_data(data)
test_synthetic_data(original, original, "test")
#> Test passed with 1 success 🥇.
#> Test passed with 1 success 🎊.
#> Test passed with 1 success 🎉.
```
