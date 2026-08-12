# Apply summary_fn across cells and persist into summary_assay

Apply summary_fn across cells and persist into summary_assay

## Usage

``` r
.apply_summary_fn(
  se,
  fit_dt,
  summary_fn,
  summary_assay,
  merge,
  fit_source,
  on_error
)
```

## Arguments

- se:

  SummarizedExperiment

- fit_dt:

  flat data.table with "row", "column", and all fit columns

- summary_fn:

  function(data.table) -\> named list

- summary_assay:

  assay name for summary output

- merge:

  merge strategy

- fit_source:

  fit_source label

- on_error:

  error handling strategy

## Value

updated SummarizedExperiment
