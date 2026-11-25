# fit_SE for combination screens

Perform fittings for combination screens.

## Usage

``` r
fit_SE.combinations(
  se,
  data_type = gDRutils::get_supported_experiments("combo"),
  series_identifiers = NULL,
  normalization_types = c("GR", "RV"),
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  score_FUN = calculate_score,
  fit_args = list()
)
```

## Arguments

- se:

  `SummarizedExperiment` object with a BumpyMatrix assay containing
  averaged data.

- data_type:

  single-agent vs combination

- series_identifiers:

  character vector of the column names in the nested `DFrame`
  corresponding to nested identifiers.

- normalization_types:

  character vector of normalization types used for calculating combo
  matrix.

- averaged_assay:

  string of the name of the averaged assay to use as input. in the `se`.

- metrics_assay:

  string of the name of the metrics assay to output in the returned
  SummarizedExperiment. whose combination represents a unique series for
  which to fit curves.

- score_FUN:

  function used to calculate score for HSA and Bliss

- fit_args:

  named list of arguments for the single-agent curve fits of the
  combination, passed to
  [`gDRutils::fit_drug_response_metrics()`](https://gdrplatform.github.io/gDRstyle/reference/fit_drug_response_metrics.html),
  e.g. `list(solver = my_solver)` to fit with another optimiser. Each
  series keeps its own `x_0`, `cap` and `force_fit`.

## Value

A `SummarizedExperiment` object with an additional assay containing the
combination metrics.

## Details

This function assumes that the combination is set up with both
concentrations nested in the assay.

## Examples

``` r
fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")

se1 <- fmae_cms[[gDRutils::get_supported_experiments("combo")]]
SummarizedExperiment::assays(se1) <-
  SummarizedExperiment::assays(se1)["Averaged"]
fit_SE.combinations(se1[1, 1])
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> class: SummarizedExperiment 
#> dim: 1 1 
#> metadata(4): identifiers experiment_metadata Keys provenance
#> assays(6): Averaged Metrics ... all_iso_points scores
#> rownames(1): G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#> rowData names(7): Gnumber DrugName ... drug_moa_2 Duration
#> colnames(1): CL00016_cellline_GB_tissue_y_46
#> colData names(4): clid CellLineName Tissue ReferenceDivisionTime
```
