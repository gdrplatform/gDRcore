# check if the given step can be skipped if partial run is chosen

check if the given step can be skipped if partial run is chosen

## Usage

``` r
do_skip_step(current_step, start_from, steps = get_pipeline_steps())
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
