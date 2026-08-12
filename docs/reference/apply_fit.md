# apply_fit

Generic layer for applying a user-supplied fit function across all (row
x column x slicing_col) triplets in a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
and persisting the results into a named output assay.

## Usage

``` r
apply_fit(
  se,
  fit_fn,
  data_type,
  slicing_cols = NULL,
  slicing_values = NULL,
  input_assay = NULL,
  output_assay,
  summary_fn = NULL,
  summary_assay = NULL,
  merge = "merge",
  on_error = "warn",
  fit_source
)

apply_custom_fit(...)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- fit_fn:

  function(`data.table`) returning a named list or single-row
  data.frame. Receives data already filtered to one `slicing_cols`
  combination.

- data_type:

  one of `"single-agent"`, `"combination"`, `"time-course"`. Resolves
  profile defaults.

- slicing_cols:

  character vector of column name(s) used to split each BumpyMatrix cell
  into sub-experiments. `NULL` uses the profile default.

- slicing_values:

  character vector of values to iterate over in `slicing_cols`. `NULL`
  iterates all unique values found in each cell.

- input_assay:

  name of the source BumpyMatrix assay. `NULL` uses the profile default.

- output_assay:

  name of the assay to write fit results into. Required. Any assay name
  is accepted, including the native `"Metrics"` assay. With
  `merge = "merge"` (default), existing rows keyed by `fit_source` are
  replaced while all other rows are preserved — safe for co-existing
  alongside native gDR metrics. With `merge = "replace"`, the entire
  assay is overwritten.

- summary_fn:

  optional function(`data.table`) → named list called once per (row x
  column) cell on all rows returned by `fit_fn` for that cell. Requires
  `summary_assay`.

- summary_assay:

  name of the assay to write summary results into. Required when
  `summary_fn` is provided.

- merge:

  `"merge"` (idempotent upsert keyed by `fit_source`

  - `slicing_cols`) or `"replace"` (overwrite the whole assay).

- on_error:

  `"warn"` (skip failed cells with a warning) or `"stop"` (propagate the
  error).

- fit_source:

  character string recorded in the `fit_source` column of every output
  row. Forms part of the upsert key.

## Value

updated `SummarizedExperiment` with `output_assay` (and optionally
`summary_assay`) added or updated.

## Details

The experiment type is declared via `data_type`, which resolves a
built-in slicing profile (default `slicing_cols`, `slicing_values`, and
`input_assay`). Any profile field can be overridden with the
corresponding argument.

An optional `summary_fn` may be provided. It is called once per (row x
column) cell and receives all rows that `fit_fn` produced for that cell,
returning a single summary row written to `summary_assay`. This is the
right place for aggregated metrics (e.g. mean synergy across
normalization types) that span multiple slices.

Use the pipe to apply several custom fits to the same SE:

      se |>
        apply_fit(bliss_fn, "combination",
                         output_assay = "custom_bliss", ...) |>
        apply_fit(musyc_fn,  "combination",
                         output_assay = "musyc_params",
                         summary_fn   = musyc_summary_fn,
                         summary_assay = "musyc_summary", ...)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
se <- mae[["single-agent"]]
mean_fn <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
se_out <- apply_fit(
  se, mean_fn, "single-agent",
  output_assay = "custom_mean", fit_source = "demo"
)
```
