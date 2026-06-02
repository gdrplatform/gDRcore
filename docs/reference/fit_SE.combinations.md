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
  score_FUN = calculate_score
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
#> Warning: overriding original x_0 argument '1' with '1' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '1' with '1' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.937720959525358' with '0.9563' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.411661143403392' with '0.4075' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '-0.466087101282658' with '-0.4678' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '-0.638711813664539' with '-0.5972' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '-0.65250302581797' with '-0.6296' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '-0.653485583857411' with '-0.692' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '-0.653555270837867' with '-0.7039' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '-0.653560191563928' with '-0.7046' (only 1 normalized value detected, setting constant fit)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: overriding original x_0 argument '1' with '1' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '1' with '1' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.960100165903769' with '0.966' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.578859775899146' with '0.577' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.123503648078634' with '0.1259' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.0689925257550951' with '0.0814' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.0658336504167383' with '0.0714' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.0656619818287103' with '0.0535' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.0656526460997796' with '0.0503' (only 1 normalized value detected, setting constant fit)
#> Warning: overriding original x_0 argument '0.0656521405542176' with '0.0501' (only 1 normalized value detected, setting constant fit)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> class: SummarizedExperiment 
#> dim: 1 1 
#> metadata(3): identifiers experiment_metadata Keys
#> assays(6): Averaged excess ... scores Metrics
#> rownames(1): G00004_drug_004_moa_A_G00021_drug_021_moa_D_72
#> rowData names(7): Gnumber DrugName ... drug_moa_2 Duration
#> colnames(1): CL00016_cellline_GB_tissue_y_46
#> colData names(4): clid CellLineName Tissue ReferenceDivisionTime
```
