# check if the given step is preceding the step chosen in the partial run

check if the given step is preceding the step chosen in the partial run

## Usage

``` r
is_preceding_step(current_step, start_from, steps = get_pipeline_steps())
```

## Arguments

- current_step, :

  string with the step to be evaluated

- start_from:

  string indicating the pipeline step from which partial run should be
  launched

- steps:

  charvect with all available steps

## Value

logical
