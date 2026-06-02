# get info about created/present assays in SE at the given pipeline step

get info about created/present assays in SE at the given pipeline step

## Usage

``` r
get_assays_per_pipeline_step(
  step,
  data_model,
  status = c("created", "present")
)
```

## Arguments

- step:

  string with pipeline step

- data_model:

  single-agent vs combination

- status:

  string return vector of assays created or present at the given step?

## Value

assay
