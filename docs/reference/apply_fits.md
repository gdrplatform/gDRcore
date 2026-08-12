# apply_fits

Apply multiple fit functions in a **single pass** over each BumpyMatrix
cell, writing each function's results into its own named output assay.

## Usage

``` r
apply_fits(
  se,
  fit_fns,
  data_type,
  slicing_cols = NULL,
  slicing_values = NULL,
  input_assay = NULL,
  merge = "merge",
  on_error = "warn",
  fit_source
)

apply_custom_fits(...)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- fit_fns:

  named list of functions. Each name becomes an output assay name; each
  value is a function(`data.table`) → named list (or a named list of
  named lists for the shared pre-computation pattern).

- data_type:

  one of `"single-agent"`, `"combination"`, `"time-course"`.

- slicing_cols:

  character vector; `NULL` uses profile default.

- slicing_values:

  character vector; `NULL` uses profile default.

- input_assay:

  string; `NULL` uses profile default.

- merge:

  `"merge"` or `"replace"`.

- on_error:

  `"warn"` or `"stop"`.

- fit_source:

  character string stamped as `fit_source` in every output row.

## Value

updated `SummarizedExperiment` with one new (or updated) assay per entry
in `fit_fns`.

## Details

`apply_fits()` is the performance-efficient alternative to chaining
multiple
[`apply_fit`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
calls when two or more fit functions operate on the **same input
assay**. Instead of unsplitting the BumpyMatrix K times (once per
function), it traverses each cell once and applies all functions in that
single pass.

Use this when:

- You have two or more independent fit functions on the same data (e.g.
  Bliss + HSS on combination data).

- A single fit function produces results for multiple output assays
  (shared pre-computation pattern — see below).

### Independent fit functions (most common)

Names of `fit_fns` become the output assay names:

      apply_fits(
        combo_se,
        fit_fns = list(
          custom_bliss = bliss_fit_fn,
          custom_hss   = hss_fit_fn
        ),
        data_type  = "combination",
        fit_source = "synergy"
      )

### Shared pre-computation (advanced)

A single function can return a **named list of named lists** to write
multiple assays while computing expensive intermediates only once:

      bliss_and_hss <- function(dt) {
        sa_curves <- fit_hill_curves(dt)  # expensive — done once
        list(
          custom_bliss = list(bliss_score = compute_bliss(sa_curves, dt)),
          custom_hss   = list(hss_score   = compute_hss(sa_curves, dt))
        )
      }
      apply_fits(
        combo_se,
        fit_fns    = list(bliss_and_hss = bliss_and_hss),
        output_assay_map = c(bliss_and_hss = NA),  # ignored; keys from return value
        ...
      )

The multi-output pattern is detected automatically when a fit function
returns a named list whose values are themselves named lists. Each
top-level name maps to an assay; the inner named list provides the row
columns.

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
se <- mae[["single-agent"]]
fn_a <- function(dt) list(x_mean = mean(dt$x, na.rm = TRUE))
fn_b <- function(dt) list(x_sd = sd(dt$x, na.rm = TRUE))
se_out <- apply_fits(
  se,
  fit_fns = list(mean_metrics = fn_a, sd_metrics = fn_b),
  data_type = "single-agent",
  fit_source = "demo"
)
```
