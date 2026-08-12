# Custom fitting with apply_fit()

``` r
library(gDRcore)
library(gDRtestData)
library(gDRutils)
library(SummarizedExperiment)
library(BumpyMatrix)
library(data.table)
```

## Overview

The standard gDR pipeline fits dose-response curves using a fixed
4-parameter log-logistic model
([`fit_SE()`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)).
Sometimes you need something different: an alternative curve model,
custom synergy metrics, Bayesian estimates, or bespoke pharmacology
metrics for your assay.

[`apply_fit()`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
and
[`apply_fits()`](https://gdrplatform.github.io/gDRcore/reference/apply_fits.md)
let you plug **any R function** into the gDR SE/MAE pipeline without
touching the pipeline internals.

Key properties:

- Works with **single-agent**, **combination**, and **time-course**
  data.
- Each fit function writes to its **own named assay** — no risk of
  overwriting native gDR assays (`"Metrics"`, `"scores"`, `"excess"`,
  …).
- Results are **idempotent**: calling twice with the same `fit_source`
  overwrites rather than duplicates.
- [`apply_fits()`](https://gdrplatform.github.io/gDRcore/reference/apply_fits.md)
  applies N fit functions in a **single** BumpyMatrix traversal —
  efficient when you have several metrics to compute on the same data.

------------------------------------------------------------------------

## Data setup

We use a small synthetic single-agent SE from `gDRtestData`.

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
sa_se <- mae[["single-agent"]]
sa_se
#> class: SummarizedExperiment 
#> dim: 10 10 
#> metadata(5): identifiers experiment_metadata Keys fit_parameters
#>   .internal
#> assays(5): RawTreated Controls Normalized Averaged Metrics
#> rownames(10): G00002_drug_002_moa_A_72 G00003_drug_003_moa_A_72 ...
#>   G00010_drug_010_moa_A_72 G00011_drug_011_moa_B_72
#> rowData names(4): Gnumber DrugName drug_moa Duration
#> colnames(10): CL00011_cellline_BA_tissue_x_26
#>   CL00012_cellline_CA_tissue_x_30 ... CL00019_cellline_JB_tissue_z_58
#>   CL00020_cellline_KB_tissue_z_62
#> colData names(4): clid CellLineName Tissue ReferenceDivisionTime
```

The `Averaged` assay is the input to custom fit functions. Each cell of
the BumpyMatrix contains one data.table per (drug × cell line) pair:

``` r
avg_cell <- BumpyMatrix::unsplitAsDataFrame(
  assay(sa_se, "Averaged"),
  row.field = "row", column.field = "column"
)
head(avg_cell[avg_cell$row == avg_cell$row[1] &
              avg_cell$column == avg_cell$column[1], ])
#> DataFrame with 6 rows and 6 columns
#>                         row                 column normalization_type
#>                 <character>            <character>           <factor>
#> 1_RV G00002_drug_002_moa_.. CL00011_cellline_BA_..                 RV
#> 1_GR G00002_drug_002_moa_.. CL00011_cellline_BA_..                 GR
#> 2_RV G00002_drug_002_moa_.. CL00011_cellline_BA_..                 RV
#> 2_GR G00002_drug_002_moa_.. CL00011_cellline_BA_..                 GR
#> 3_RV G00002_drug_002_moa_.. CL00011_cellline_BA_..                 RV
#> 3_GR G00002_drug_002_moa_.. CL00011_cellline_BA_..                 GR
#>      Concentration         x      x_std
#>          <numeric> <numeric>  <numeric>
#> 1_RV    0.00100000  0.924967 0.01066130
#> 1_GR    0.00100000  0.944433 0.00805874
#> 2_RV    0.00316228  0.739100 0.01484621
#> 2_GR    0.00316228  0.793100 0.01302152
#> 3_RV    0.01000000  0.435933 0.03071503
#> 3_GR    0.01000000  0.481333 0.03811841
```

Each cell has four columns:

| Column               | Description                                                   |
|----------------------|---------------------------------------------------------------|
| `normalization_type` | `"GR"` or `"RV"`                                              |
| `Concentration`      | drug concentration (µM)                                       |
| `x`                  | averaged normalized response (GR value or relative viability) |
| `x_std`              | standard deviation across replicates                          |

------------------------------------------------------------------------

## The fit_fn contract

A fit function must satisfy this interface:

    fit_fn(avg_dt) -> named list (or single-row data.frame / data.table)

`avg_dt` is a **data.table** containing the rows for **one** combination
of slicing column values (by default: one `normalization_type`) for one
(drug × cell line) cell.

The returned named list becomes one row in the output assay. Column
names are the names of the list elements. `fit_source` is stamped
automatically by the generic layer — do not include it in the return
value.

### Minimal example

``` r
# The simplest possible fit_fn: compute mean and SD of the response
summary_fn <- function(avg_dt) {
  list(
    x_mean = mean(avg_dt$x, na.rm = TRUE),
    x_sd = sd(avg_dt$x, na.rm = TRUE),
    n = NROW(avg_dt)
  )
}
```

------------------------------------------------------------------------

## Single-agent: apply_fit()

Apply the summary function to every (drug × cell line ×
normalization_type) triplet and write results to a custom assay named
`"custom_summary"`.

``` r
sa_out <- apply_fit(
  sa_se,
  fit_fn = summary_fn,
  data_type = "single-agent",
  output_assay = "custom_summary",
  fit_source = "demo"
)
assayNames(sa_out)
#> [1] "RawTreated"     "Controls"       "Normalized"     "Averaged"      
#> [5] "Metrics"        "custom_summary"
```

``` r
summary_df <- BumpyMatrix::unsplitAsDataFrame(
  assay(sa_out, "custom_summary"),
  row.field = "row", column.field = "column"
)
head(summary_df)
#> DataFrame with 6 rows and 7 columns
#>                      row                 column    x_mean      x_sd         n
#>              <character>            <character> <numeric> <numeric> <integer>
#> 1 G00002_drug_002_moa_.. CL00011_cellline_BA_..  0.491452  0.221336         9
#> 2 G00002_drug_002_moa_.. CL00011_cellline_BA_..  0.466411  0.214770         9
#> 3 G00003_drug_003_moa_.. CL00011_cellline_BA_..  0.417704  0.273562         9
#> 4 G00003_drug_003_moa_.. CL00011_cellline_BA_..  0.419170  0.256418         9
#> 5 G00004_drug_004_moa_.. CL00011_cellline_BA_..  0.843222  0.107986         9
#> 6 G00004_drug_004_moa_.. CL00011_cellline_BA_..  0.803893  0.133827         9
#>    fit_source normalization_type
#>   <character>        <character>
#> 1        demo                 GR
#> 2        demo                 RV
#> 3        demo                 GR
#> 4        demo                 RV
#> 5        demo                 GR
#> 6        demo                 RV
```

One row per (drug × cell line × normalization_type) triplet. The native
`"Metrics"` assay is untouched — because we chose a custom assay name.

### Writing to the Metrics assay

You can also write directly to `"Metrics"` — for example when replacing
the standard gDR Hill fit with your own model:

``` r
# SE straight from the pipeline — already has a native "Metrics" assay
# (fit_source = "gDR")
"Metrics" %in% assayNames(sa_se)
#> [1] TRUE

# Apply a custom fit to the same assay, coexisting alongside gDR rows
custom_hill <- apply_fit(
  sa_se,
  fit_fn = fit_drug_response_metrics,
  data_type = "single-agent",
  output_assay = "Metrics",
  fit_source = "custom_hill"   # distinct key keeps native "gDR" rows intact
)

metrics_df <- BumpyMatrix::unsplitAsDataFrame(
  assay(custom_hill, "Metrics"),
  row.field = "row", column.field = "column"
)
unique(metrics_df$fit_source)   # both "gDR" and "custom_hill"
#> [1] "gDR"         "custom_hill"
```

With `merge = "merge"` (default) the upsert key is
`fit_source + normalization_type`, so native rows (`fit_source = "gDR"`)
are preserved. Use `merge = "replace"` only if you intend to overwrite
the whole assay.

### Idempotent merge

Calling again with the same `fit_source` **replaces** those rows rather
than appending:

``` r
n_before <- NROW(BumpyMatrix::unsplitAsDataFrame(
  assay(sa_out, "custom_summary"),
  row.field = "row", column.field = "column"
))

# Call again — same fit_source, same data
sa_out2 <- apply_fit(
  sa_out,
  fit_fn = summary_fn,
  data_type = "single-agent",
  output_assay = "custom_summary",
  fit_source = "demo"
)
n_after <- NROW(BumpyMatrix::unsplitAsDataFrame(
  assay(sa_out2, "custom_summary"),
  row.field = "row", column.field = "column"
))

stopifnot(n_before == n_after)   # no duplicate rows
message("Row count before: ", n_before, " — after: ", n_after, " (no change)")
```

### Coexisting fit sources

Different `fit_source` values live **side by side** in the same assay:

``` r
extra_fn <- function(avg_dt) list(x_max = max(avg_dt$x, na.rm = TRUE))

sa_two <- sa_out |>
  apply_fit(extra_fn, "single-agent",
                   output_assay = "custom_summary",
                   fit_source = "extremes")

sources <- unique(BumpyMatrix::unsplitAsDataFrame(
  assay(sa_two, "custom_summary"),
  row.field = "row", column.field = "column"
)$fit_source)
message("fit_source values in assay: ", paste(sources, collapse = ", "))
```

------------------------------------------------------------------------

## Reference Hill fit: fit_drug_response_metrics()

[`fit_drug_response_metrics()`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics.md)
is a reference single-agent fit function that replicates the standard
[`fit_SE()`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)
/
[`logisticFit()`](https://gdrplatform.github.io/gDRstyle/reference/logisticFit.html)
output exactly.

| Property         | Value                                                                                                                                                                                                   |
|------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Model            | 3-parameter log-logistic ([`drc::LL.3u`](https://rdrr.io/pkg/drc/man/LL.3.html), `x_0` fixed at 1)                                                                                                      |
| Equivalent to    | [`fit_SE()`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md) / [`gDRutils::logisticFit()`](https://gdrplatform.github.io/gDRstyle/reference/logisticFit.html) |
| `x_mean`         | Predicted from fitted curve (matches `logisticFit` behaviour)                                                                                                                                           |
| `fit_type` value | `"DRC3pHillFitModelFixS0"` or `"DRCConstantFitResult"`                                                                                                                                                  |
| Output columns   | Full `Metrics` assay schema including `p_value`, `rss`, `x_AOC_range`, `x_max`, `x_sd_avg`                                                                                                              |

### Numerical equivalence with fit_SE()

[`fit_drug_response_metrics()`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics.md)
is **numerically identical** to
[`fit_SE()`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md).
The fit is deterministic —
[`drc::drm`](https://rdrr.io/pkg/drc/man/drm.html) does not use random
number generation, so [`set.seed()`](https://rdrr.io/r/base/Random.html)
has no effect. Differences between the two are purely algorithmic:

- Same [`drc::LL.3u`](https://rdrr.io/pkg/drc/man/LL.3.html) model,
  `x_0 = 1`, identical priors and concentration bounds
- Same `x_mean`: predicted from the fitted curve over the observed
  concentration range
- Same `pcutoff = 0.05` fallback: when the F-test gives
  `p_value ≥ pcutoff`, the result is replaced by a flat
  `DRCConstantFitResult` (same as
  [`logisticFit()`](https://gdrplatform.github.io/gDRstyle/reference/logisticFit.html))
- Same `n_point_cutoff = 4`: fewer unique concentrations → constant fit
  without attempting the sigmoidal model
- All output columns present in the native `"Metrics"` assay

These parameters can be overridden if your data requires it:

``` r
# Loosen the significance threshold or force the sigmoidal fit
fit_drug_response_metrics(avg_dt, pcutoff = 0.1)
fit_drug_response_metrics(avg_dt, force_fit = TRUE)
# Different range for x_AOC_range computation
fit_drug_response_metrics(avg_dt, range_conc = c(1e-3, 10))
```

``` r
hill_out <- apply_fit(
  sa_se,
  fit_fn = fit_drug_response_metrics,
  data_type = "single-agent",
  output_assay = "custom_hill",
  fit_source = "hill_ref"
)

hill_df <- BumpyMatrix::unsplitAsDataFrame(
  assay(hill_out, "custom_hill"),
  row.field = "row", column.field = "column"
)
head(hill_df[, c("row", "column", "normalization_type",
                  "ec50", "xc50", "h", "r2", "fit_type")])
#> DataFrame with 6 rows and 8 columns
#>                      row                 column normalization_type       ec50
#>              <character>            <character>        <character>  <numeric>
#> 1 G00002_drug_002_moa_.. CL00011_cellline_BA_..                 GR 0.00451256
#> 2 G00002_drug_002_moa_.. CL00011_cellline_BA_..                 RV 0.00379589
#> 3 G00003_drug_003_moa_.. CL00011_cellline_BA_..                 GR 0.00471408
#> 4 G00003_drug_003_moa_.. CL00011_cellline_BA_..                 RV 0.00396299
#> 5 G00004_drug_004_moa_.. CL00011_cellline_BA_..                 GR 0.01715280
#> 6 G00004_drug_004_moa_.. CL00011_cellline_BA_..                 RV 0.01654638
#>         xc50         h        r2               fit_type
#>    <numeric> <numeric> <numeric>            <character>
#> 1 0.00895622   1.91190  0.988441 DRC3pHillFitModelFixS0
#> 2 0.00725317   1.82703  0.992470 DRC3pHillFitModelFixS0
#> 3 0.00659147   2.27009  0.990529 DRC3pHillFitModelFixS0
#> 4 0.00565558   2.34960  0.995851 DRC3pHillFitModelFixS0
#> 5        Inf   3.14371  0.991428 DRC3pHillFitModelFixS0
#> 6        Inf   3.15994  0.991887 DRC3pHillFitModelFixS0
```

`ec50` is the raw model parameter (concentration at half-maximal
effect); `xc50` is the capped version (capped at
`capping_fold × max(Concentration)`) as reported in the native
`"Metrics"` assay.

------------------------------------------------------------------------

## summary_fn: cell-level aggregation

An optional `summary_fn` is called **once per (drug × cell line) cell**
on all rows produced by `fit_fn` for that cell (one row per
normalization type). It is the right place for metrics that aggregate
across normalization types — for example, whether the result is
synergistic across both GR and RV.

``` r
# Aggregate: mean xc50 and flag whether both norm types fitted successfully
hill_summary_fn <- function(fit_dt) {
  list(
    mean_xc50 = mean(fit_dt$xc50, na.rm = TRUE),
    mean_r2 = mean(fit_dt$r2, na.rm = TRUE),
    all_converged = all(fit_dt$fit_type == "DRC3pHillFitModelFixS0", na.rm = TRUE)
  )
}

hill_with_summary <- apply_fit(
  sa_se,
  fit_fn = fit_drug_response_metrics,
  data_type = "single-agent",
  output_assay = "custom_hill",
  summary_fn = hill_summary_fn,
  summary_assay = "custom_hill_summary",
  fit_source = "hill_ref"
)
assayNames(hill_with_summary)
#> [1] "RawTreated"          "Controls"            "Normalized"         
#> [4] "Averaged"            "Metrics"             "custom_hill"        
#> [7] "custom_hill_summary"
```

``` r
sumdf <- BumpyMatrix::unsplitAsDataFrame(
  assay(hill_with_summary, "custom_hill_summary"),
  row.field = "row", column.field = "column"
)
head(sumdf[, c("row", "column", "mean_xc50", "mean_r2", "all_converged")])
#> DataFrame with 6 rows and 5 columns
#>                      row                 column  mean_xc50   mean_r2
#>              <character>            <character>  <numeric> <numeric>
#> 1 G00002_drug_002_moa_.. CL00011_cellline_BA_.. 0.00810470  0.990455
#> 2 G00003_drug_003_moa_.. CL00011_cellline_BA_.. 0.00612352  0.993190
#> 3 G00004_drug_004_moa_.. CL00011_cellline_BA_..        Inf  0.991658
#> 4 G00005_drug_005_moa_.. CL00011_cellline_BA_..        Inf  0.977308
#> 5 G00006_drug_006_moa_.. CL00011_cellline_BA_..        Inf  0.997368
#> 6 G00007_drug_007_moa_.. CL00011_cellline_BA_..        Inf  0.850211
#>   all_converged
#>       <logical>
#> 1          TRUE
#> 2          TRUE
#> 3          TRUE
#> 4          TRUE
#> 5          TRUE
#> 6          TRUE
```

One row per (drug × cell line) — regardless of how many normalization
types were fitted.

------------------------------------------------------------------------

## Combination data: synergy scores

### What fit_SE.combinations() does internally

Before describing the new extension API, it helps to understand what the
standard
[`fit_SE.combinations()`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
function does. It is a single loop over each (drug-combo × cell-line)
pair that executes **five sequential steps**, each writing to a separate
assay:

| Step | Key functions                                                                                                                                                                                                                                                                                   | Output assay                     | Description                                                                                                                                  |
|------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| 1    | `fit_combo_cotreatments()`, `fit_combo_codilutions()`                                                                                                                                                                                                                                           | `Metrics`                        | Fit SA dose-response curves at each co-treatment concentration; produces `ec50`, `h`, `x_inf`, `x_0` per SA series                           |
| 2    | [`map_ids_to_fits()`](https://gdrplatform.github.io/gDRcore/reference/map_ids_to_fits.md)                                                                                                                                                                                                       | *(internal)*                     | Predict **smooth** single-agent responses at every combo concentration using the SA fits from step 1; average col/row/codilution predictions |
| 3    | [`calculate_HSA()`](https://gdrplatform.github.io/gDRcore/reference/calculate_matrix_metric.md), [`calculate_Bliss()`](https://gdrplatform.github.io/gDRcore/reference/calculate_matrix_metric.md), [`calculate_excess()`](https://gdrplatform.github.io/gDRcore/reference/calculate_excess.md) | `excess`                         | Compute expected response (HSA = min of SAs; Bliss = product/GR formula); compute per-point excess = expected − observed                     |
| 4    | `calculate_Loewe()`                                                                                                                                                                                                                                                                             | `isobolograms`, `all_iso_points` | Compute combination index (CI) via isobologram analysis; CI \< 1 = synergy                                                                   |
| 5    | [`calculate_score()`](https://gdrplatform.github.io/gDRcore/reference/calculate_score.md)                                                                                                                                                                                                       | `scores`                         | Reduce per-point excess to a scalar: mean of top-10-percentile values → `bliss_score`, `hsa_score`, `CIScore_50`, `CIScore_80`               |

**The key insight:** step 2 (smooth) depends on step 1 (SA fits), and
steps 3–5 all depend on step 2. The boundaries between steps are clean,
which means they can be extracted as independent public functions.

**Roadmap (GDR-3486):** a future refactoring will expose each step as a
standalone `apply_combo_*()` function, making
[`fit_SE.combinations()`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
a thin wrapper over them.
[`apply_combo_scores()`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_scores.md)
(step 5) is already available as part of GDR-3352.

------------------------------------------------------------------------

### Two approaches for combination scoring in the new API

Two approaches are available, depending on whether you have SA fits or
only raw averaged data.

### Approach 1: apply_combo_scores() — replicate fit_SE.combinations exactly

[`apply_combo_scores()`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_scores.md)
is a high-level function that reproduces the Bliss and HSA scoring logic
of
[`fit_SE.combinations()`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
exactly, using fitted SA curves from the `Metrics` assay to generate
smooth single-agent predictions.

| Property      | Value                                                                                                                  |
|---------------|------------------------------------------------------------------------------------------------------------------------|
| Requires      | `Averaged` + `Metrics` assay (with `dilution_drug`, `ec50`, `h`, …)                                                    |
| Equivalent to | [`fit_SE.combinations()`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md) Bliss and HSA scores |
| Accuracy      | cor \> 0.998 with `fit_SE.combinations` on real data                                                                   |
| Use when      | Replacing or extending `fit_SE.combinations` for standard data                                                         |

``` r
# Use the small synthetic combo dataset which has both Averaged and Metrics
combo_mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
combo_name <- gDRutils::get_supported_experiments("combo")
combo_se_full <- combo_mae[[combo_name]]

# combo_se_full already has Metrics from fit_SE.combinations
combo_scored <- apply_combo_scores(combo_se_full)
assayNames(combo_scored)
#> [1] "RawTreated"     "Controls"       "Normalized"     "Averaged"      
#> [5] "excess"         "all_iso_points" "isobolograms"   "scores"        
#> [9] "Metrics"
```

``` r
scores_df <- BumpyMatrix::unsplitAsDataFrame(
  assay(combo_scored, "scores"),
  row.field = "row", column.field = "column"
)
scores_df[, c("row", "column", "normalization_type", "bliss_score", "hsa_score")]
#> DataFrame with 24 rows and 5 columns
#>                        row                 column normalization_type
#>                <character>            <character>        <character>
#> 1   G00004_drug_004_moa_.. CL00016_cellline_GB_..                 GR
#> 2   G00004_drug_004_moa_.. CL00016_cellline_GB_..                 RV
#> 3   G00004_drug_004_moa_.. CL00016_cellline_GB_..                 GR
#> 4   G00004_drug_004_moa_.. CL00016_cellline_GB_..                 RV
#> 5   G00005_drug_005_moa_.. CL00016_cellline_GB_..                 GR
#> ...                    ...                    ...                ...
#> 20  G00005_drug_005_moa_.. CL00017_cellline_HB_..                 RV
#> 21  G00006_drug_006_moa_.. CL00017_cellline_HB_..                 GR
#> 22  G00006_drug_006_moa_.. CL00017_cellline_HB_..                 RV
#> 23  G00006_drug_006_moa_.. CL00017_cellline_HB_..                 GR
#> 24  G00006_drug_006_moa_.. CL00017_cellline_HB_..                 RV
#>     bliss_score   hsa_score
#>       <numeric>   <numeric>
#> 1   0.008339057 0.008339057
#> 2   0.002358668 0.002358668
#> 3   0.000674575 0.000674575
#> 4   0.000298478 0.000298478
#> 5   0.009101157 0.009101157
#> ...         ...         ...
#> 20  0.000151856 0.000151856
#> 21  0.031880296 0.031880296
#> 22  0.012769276 0.012769276
#> 23  0.024087057 0.024087057
#> 24  0.005341883 0.005341883
```

`bliss_score > 0` and `hsa_score > 0` indicate synergy.

### Approach 2: bliss_fit_fn / hss_fit_fn — simplified, no SA fits needed

[`bliss_fit_fn()`](https://gdrplatform.github.io/gDRcore/reference/bliss_fit_fn.md)
and
[`hss_fit_fn()`](https://gdrplatform.github.io/gDRcore/reference/hss_fit_fn.md)
are lower-level fit functions for
[`apply_fit()`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
that compute synergy scores directly from the raw Averaged data, without
requiring prior SA curve fits.

| Property    | Value                                                       |
|-------------|-------------------------------------------------------------|
| Requires    | `Averaged` assay only                                       |
| SA response | Raw single-agent edge points (no curve smoothing)           |
| Use when    | Prototyping, custom models, or when SA fits are unavailable |

``` r
# Build a minimal synthetic combination SE (Averaged only, no Metrics needed)
combo_dt <- data.table::CJ(
  row = c("DrugA", "DrugB"),
  column = "CellLine1",
  normalization_type = c("GR", "RV"),
  Concentration = c(0, 0.1, 1.0),
  Concentration_2 = c(0, 0.1, 1.0)
)
set.seed(42L)
combo_dt[, x := pmax(0.05,
  1 - 0.3 * Concentration / (Concentration + 0.5) -
      0.2 * Concentration_2 / (Concentration_2 + 0.5) +
      rnorm(.N, 0, 0.03))]

data_cols <- setdiff(names(combo_dt), c("row", "column"))
combo_bumpy <- BumpyMatrix::splitAsBumpyMatrix(
  combo_dt[, data_cols, with = FALSE],
  row = combo_dt$row, col = combo_dt$column
)
combo_se <- SummarizedExperiment(assays = list(Averaged = combo_bumpy))
```

``` r
bliss_out <- apply_fit(
  combo_se,
  fit_fn = bliss_fit_fn,
  data_type = "combination",
  output_assay = "custom_bliss",
  fit_source = "bliss"
)

bliss_df <- BumpyMatrix::unsplitAsDataFrame(
  assay(bliss_out, "custom_bliss"),
  row.field = "row", column.field = "column"
)
bliss_df[, c("row", "column", "normalization_type",
             "bliss_score", "bliss_excess_mean", "n_combo_points")]
#> DataFrame with 4 rows and 6 columns
#>           row      column normalization_type bliss_score bliss_excess_mean
#>   <character> <character>        <character>   <numeric>         <numeric>
#> 1       DrugA   CellLine1                 GR   0.0404576        0.02194040
#> 2       DrugA   CellLine1                 RV   0.1790833        0.07174050
#> 3       DrugB   CellLine1                 GR   0.1083000        0.02410753
#> 4       DrugB   CellLine1                 RV   0.0473357        0.00199227
#>   n_combo_points
#>        <integer>
#> 1              4
#> 2              4
#> 3              4
#> 4              4
```

``` r
hss_out <- apply_fit(
  combo_se,
  fit_fn = hss_fit_fn,
  data_type = "combination",
  output_assay = "custom_hss",
  fit_source = "hss"
)

hss_df <- BumpyMatrix::unsplitAsDataFrame(
  assay(hss_out, "custom_hss"),
  row.field = "row", column.field = "column"
)
hss_df[, c("row", "column", "normalization_type",
           "hss_score", "hss_excess_mean")]
#> DataFrame with 4 rows and 5 columns
#>           row      column normalization_type hss_score hss_excess_mean
#>   <character> <character>        <character> <numeric>       <numeric>
#> 1       DrugA   CellLine1                 GR  0.118126       0.0711634
#> 2       DrugA   CellLine1                 RV  0.232106       0.0971858
#> 3       DrugB   CellLine1                 GR  0.197907       0.0726057
#> 4       DrugB   CellLine1                 RV  0.166576       0.0480229
```

------------------------------------------------------------------------

## Efficient multi-fit: apply_fits()

When several fit functions operate on the **same input assay**, use
[`apply_fits()`](https://gdrplatform.github.io/gDRcore/reference/apply_fits.md)
to traverse each BumpyMatrix cell **once** and apply all functions in
that single pass.

``` r
combo_multi <- apply_fits(
  combo_se,
  fit_fns = list(
    custom_bliss = bliss_fit_fn,
    custom_hss = hss_fit_fn
  ),
  data_type = "combination",
  fit_source = "synergy_panel"
)
assayNames(combo_multi)
#> [1] "Averaged"     "custom_bliss" "custom_hss"
```

Both assays are written in one traversal — equivalent to chaining two
[`apply_fit()`](https://gdrplatform.github.io/gDRcore/reference/apply_fit.md)
calls but without the overhead of a second unsplit + iteration.

### Shared pre-computation pattern

When two metrics share an expensive intermediate (e.g. fitted
single-agent curves), a single fit function can return a **named list of
named lists** to populate multiple assays from one computation:

``` r
# Each top-level name maps to an output assay; inner lists are the rows
bliss_and_hss_combined <- function(dt) {
  # Expensive step done ONCE per cell
  sa1 <- dt[dt$Concentration_2 == 0 & dt$Concentration > 0, ]
  sa2 <- dt[dt$Concentration == 0   & dt$Concentration_2 > 0, ]

  list(
    custom_bliss = list(bliss_score = mean(sa1$x) - mean(sa2$x)), # simplified
    custom_hss = list(hss_score = min(c(sa1$x, sa2$x)))
  )
}

apply_fits(
  combo_se,
  fit_fns = list(custom_bliss = bliss_and_hss_combined,
                    custom_hss = bliss_and_hss_combined),
  data_type = "combination",
  fit_source = "shared"
)
```

------------------------------------------------------------------------

## Chaining with pipe

The functions are pipe-friendly — each call returns the updated SE:

``` r
result_se <- combo_se |>
  apply_fit(bliss_fit_fn, "combination",
                   output_assay = "custom_bliss",
                   fit_source = "bliss") |>
  apply_fit(hss_fit_fn, "combination",
                   output_assay = "custom_hss",
                   fit_source = "hss") |>
  apply_fit(
    function(dt) list(n_obs = NROW(dt)),
    "combination",
    output_assay = "combo_diagnostics",
    fit_source = "qc"
  )

assayNames(result_se)
#> [1] "Averaged"          "custom_bliss"      "custom_hss"       
#> [4] "combo_diagnostics"
```

------------------------------------------------------------------------

## Error handling

By default, a failed cell emits a warning and is skipped
(`on_error = "warn"`). Use `on_error = "stop"` to halt immediately and
propagate the error — useful when debugging a new fit function.

``` r
buggy_fn <- function(dt) {
  if (dt$normalization_type[1] == "GR") stop("GR not supported")
  list(x_rv = mean(dt$x, na.rm = TRUE))
}

se_partial <- withCallingHandlers(
  apply_fit(
    sa_se, buggy_fn, "single-agent",
    output_assay = "rv_only",
    fit_source = "rv_fn",
    on_error = "warn"
  ),
  warning = function(w) {
    message("[caught] ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
# Only RV rows are written; GR cells were skipped with a warning
rv_df <- BumpyMatrix::unsplitAsDataFrame(
  assay(se_partial, "rv_only"),
  row.field = "row", column.field = "column"
)
unique(rv_df$normalization_type)
#> [1] "RV"
```

------------------------------------------------------------------------

## Quick reference

### Function signatures

``` r
# Single fit → one output assay
apply_fit(
  se,
  fit_fn,
  data_type = "single-agent", # or "combination", "time-course"
  slicing_cols = NULL, # NULL → data_type default
  slicing_values = NULL, # NULL → all unique values found
  input_assay = NULL, # NULL → data_type default ("Averaged")
  output_assay, # REQUIRED — your assay name
  summary_fn = NULL, # optional cell-level aggregator
  summary_assay = NULL,
  merge = "merge", # or "replace"
  on_error = "warn", # or "stop"
  fit_source                         # REQUIRED — upsert key tag
)

# Multiple fits → one BumpyMatrix pass
apply_fits(
  se,
  fit_fns, # named list: name = output assay, value = fit function
  data_type = "single-agent",
  fit_source,
  ...
)
```

### fit_fn contract

|                     |                                                                                                         |
|---------------------|---------------------------------------------------------------------------------------------------------|
| **Input**           | `data.table` — one BumpyMatrix cell, filtered to one `slicing_cols` value                               |
| **Output**          | Named list → one row in `output_assay`; or named list of named lists for multi-assay pattern            |
| **fit_source**      | Stamped by the generic layer — do **not** include in the return value                                   |
| **slicing columns** | The caller’s value (e.g. `normalization_type`) IS available in the data.table — no need to filter again |

### summary_fn contract

|                 |                                                                                           |
|-----------------|-------------------------------------------------------------------------------------------|
| **Input**       | `data.table` — all rows written by `fit_fn` for one (row × column) cell                   |
| **Output**      | Named list → one row in `summary_assay`                                                   |
| **When to use** | Aggregated metrics that span normalization types (e.g. “converged in at least one type?”) |

### Reference implementations

#### Single-agent fit functions (for `apply_fit()`)

| Function                                                                                                            | Model               | Equivalent to                                                                                                                                                                                 | Key output columns                                                                  |
|---------------------------------------------------------------------------------------------------------------------|---------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|
| [`fit_drug_response_metrics()`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics.md)       | 3p LL.3u, `x_0 = 1` | [`fit_SE()`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md) / [`logisticFit()`](https://gdrplatform.github.io/gDRstyle/reference/logisticFit.html) | `ec50`, `xc50`, `h`, `r2`, `x_mean`, `x_AOC`, `fit_type = "DRC3pHillFitModelFixS0"` |
| [`fit_drug_response_metrics_4p()`](https://gdrplatform.github.io/gDRcore/reference/fit_drug_response_metrics_4p.md) | 4p LL.4, `x_0` free | — (extended variant)                                                                                                                                                                          | `ec50`, `xc50`, `h`, `r2`, `x_0`, `x_mean`, `fit_type = "DRC4pHillFitModel"`        |

#### Combination scoring

| Function                                                                            | Level                  | Requires               | Equivalent to                                                                                                 | Key output columns                                   |
|-------------------------------------------------------------------------------------|------------------------|------------------------|---------------------------------------------------------------------------------------------------------------|------------------------------------------------------|
| `apply_combo_scores(se)`                                                            | SE-level (recommended) | `Averaged` + `Metrics` | [`fit_SE.combinations()`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md) Bliss & HSA | `bliss_score`, `hsa_score`                           |
| [`bliss_fit_fn()`](https://gdrplatform.github.io/gDRcore/reference/bliss_fit_fn.md) | triplet `fit_fn`       | `Averaged` only        | — (simplified, no SA fits)                                                                                    | `bliss_score`, `bliss_excess_mean`, `n_combo_points` |
| [`hss_fit_fn()`](https://gdrplatform.github.io/gDRcore/reference/hss_fit_fn.md)     | triplet `fit_fn`       | `Averaged` only        | — (simplified, no SA fits)                                                                                    | `hss_score`, `hss_excess_mean`, `n_combo_points`     |

**When to use which:** - Use
[`apply_combo_scores()`](https://gdrplatform.github.io/gDRcore/reference/apply_combo_scores.md)
when you have a fully fitted SE (after
[`fit_SE.combinations()`](https://gdrplatform.github.io/gDRcore/reference/fit_SE.combinations.md)
or
[`apply_fit_to_se()`](https://gdrplatform.github.io/gDRcore/reference/apply_fit_to_se.md))
and want scores numerically consistent with the standard gDR pipeline. -
Use
[`bliss_fit_fn()`](https://gdrplatform.github.io/gDRcore/reference/bliss_fit_fn.md)
/
[`hss_fit_fn()`](https://gdrplatform.github.io/gDRcore/reference/hss_fit_fn.md)
when prototyping a new scoring approach, when SA fits are unavailable,
or when embedding score computation inside a larger custom `fit_fn`.

------------------------------------------------------------------------

## SessionInfo

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] data.table_1.18.4           BumpyMatrix_1.20.0         
#>  [3] SummarizedExperiment_1.42.0 Biobase_2.72.0             
#>  [5] GenomicRanges_1.64.0        Seqinfo_1.2.0              
#>  [7] IRanges_2.46.0              S4Vectors_0.50.1           
#>  [9] BiocGenerics_0.58.1         generics_0.1.4             
#> [11] MatrixGenerics_1.24.0       matrixStats_1.5.0          
#> [13] gDRutils_1.10.0             gDRtestData_1.10.0         
#> [15] gDRcore_1.11.8              BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] xfun_0.60                   bslib_0.12.0               
#>  [3] htmlwidgets_1.6.4           lattice_0.22-9             
#>  [5] tools_4.6.1                 sandwich_3.1-3             
#>  [7] parallel_4.6.1              drc_3.0-1                  
#>  [9] Matrix_1.7-5                checkmate_2.3.4            
#> [11] RColorBrewer_1.1-3          desc_1.4.3                 
#> [13] RcppParallel_6.2.0          lifecycle_1.0.5            
#> [15] farver_2.1.2                compiler_4.6.1             
#> [17] textshaping_1.0.5           codetools_0.2-20           
#> [19] carData_3.0-6               htmltools_0.5.9            
#> [21] sass_0.4.10                 yaml_2.3.12                
#> [23] Formula_1.2-6               pkgdown_2.2.1              
#> [25] car_3.1-5                   jquerylib_0.1.4            
#> [27] MASS_7.3-65                 BiocParallel_1.46.0        
#> [29] DelayedArray_0.38.2         cachem_1.1.0               
#> [31] abind_1.4-8                 multcomp_1.4-31            
#> [33] qs2_0.2.2                   gtools_3.9.5               
#> [35] digest_0.6.39               mvtnorm_1.4-2              
#> [37] bookdown_0.47               splines_4.6.1              
#> [39] fastmap_1.2.0               grid_4.6.1                 
#> [41] cli_3.6.6                   SparseArray_1.12.2         
#> [43] S4Arrays_1.12.0             survival_3.8-6             
#> [45] TH.data_1.1-5               scales_1.4.0               
#> [47] backports_1.5.1             plotrix_3.8-14             
#> [49] rmarkdown_2.31              XVector_0.52.0             
#> [51] otel_0.2.0                  zoo_1.9-0                  
#> [53] ragg_1.5.2                  stringfish_0.19.2          
#> [55] evaluate_1.0.5              knitr_1.51                 
#> [57] MultiAssayExperiment_1.38.0 rlang_1.3.0                
#> [59] Rcpp_1.1.2                  glue_1.8.1                 
#> [61] BiocManager_1.30.27         jsonlite_2.0.0             
#> [63] R6_2.6.1                    systemfonts_1.3.2          
#> [65] fs_2.1.0
```
