# List the measurement timepoints entering each growth-rate window

Reports which durations actually go into the fit for every period, so
that a report or a log can state the measurements it used rather than
only the window bounds. The bounds alone are ambiguous: whether a
timepoint sitting exactly on a boundary is included depends on the
window convention, and a window may silently contain far fewer points
than intended.

## Usage

``` r
get_period_timepoints(se, periods, lfc_assay = NULL)
```

## Arguments

- se:

  `SummarizedExperiment` with `lfc_assay`.

- periods:

  named list of two-element numeric windows, as passed to
  [`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md).

- lfc_assay:

  string; assay to read durations from. `NULL` uses the default input
  assay of the `"time-course"` fit profile.

## Value

`data.table` with one row per period: `period`, `window` (the bounds as
`"start - end"`), `n` (number of distinct timepoints) and `timepoints`
(those timepoints, sorted, comma-separated). Periods matching no
measurement yield `n = 0`.

`timepoints` is the union of the durations present anywhere in `se`, not
a per-plate or per-well listing. On a uniform sampling grid the two
coincide; if plates were read on different schedules, a duration
contributed by only one barcode still appears here, and `n` is therefore
an upper bound on what any single well contributed to the fit.

## Details

Uses the same half-open `[start, end)` rule as
[`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md),
so the answer always matches what was fitted.

## See also

[`compute_growth_rates`](https://gdrplatform.github.io/gDRcore/reference/compute_growth_rates.md),
[`fit_SE.timecourse`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.timecourse.md)

## Examples

``` r
if (FALSE) { # \dontrun{
get_period_timepoints(se_tc, list(early = c(44, 92), late = c(92, 140)))
} # }
```
