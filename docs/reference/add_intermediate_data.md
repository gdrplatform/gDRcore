# add intermediate data (qs2 files) for given MAE

add intermediate data (qs2 files) for given MAE

## Usage

``` r
add_intermediate_data(mae, data_dir, steps = get_pipeline_steps())
```

## Arguments

- mae:

  `MultiAssayExperiment` with dose-response data

- data_dir:

  output directory

- steps:

  character vector with pipeline steps for which intermediate data
  should be saved

## Value

`NULL`
