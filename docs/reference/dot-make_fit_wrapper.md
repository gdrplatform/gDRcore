# Build a per-cell wrapper for apply_bumpy_function

Iterates over each unique combination of slicing_col values, filters the
cell data.table, calls fit_fn, stamps fit_source and slicing column
values, and rbinds all results. Respects on_error semantics.

## Usage

``` r
.make_fit_wrapper(fit_fn, slicing_cols, slicing_values, on_error, fit_source)
```

## Arguments

- fit_fn:

  user-supplied fit function

- slicing_cols:

  character vector of column names to slice by

- slicing_values:

  character vector of values to iterate, or NULL for all

- on_error:

  `"warn"` or `"stop"`

- fit_source:

  fit_source label

## Value

function(DFrame) -\> data.table
