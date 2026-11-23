# Add the curve fits delivered by PRISM to the Metrics assay

Adds the curves PRISM fitted (read with
[`gDRimport::read_PRISM_drc()`](https://gdrplatform.github.io/gDRimport/reference/read_PRISM_drc.html))
to the `Metrics` assay of processed PRISM data, with
`fit_source = "PRISM"`, next to the gDR fits. The curves are not
refitted, so the values match the PRISM delivery exactly.

## Usage

``` r
add_PRISM_metrics(
  se,
  prism_curves,
  drug_map = NULL,
  normalization_type = "RV",
  metrics_assay = "Metrics",
  averaged_assay = "Averaged"
)
```

## Arguments

- se:

  `SummarizedExperiment` of single-agent PRISM data, processed by gDR.

- prism_curves:

  `data.table` of PRISM curves, as returned by
  [`gDRimport::read_PRISM_drc()`](https://gdrplatform.github.io/gDRimport/reference/read_PRISM_drc.html):
  the columns of the DRC table, the cell line identifier and Duration.

- drug_map:

  named character vector mapping PRISM's `pert_name` to the drug
  identifier in `se`, for data imported with other drug identifiers. By
  default the `pert_name` is the drug identifier, as in the PRISM LEVEL5
  import.

- normalization_type:

  string, the normalization type of the PRISM curves.

- metrics_assay:

  string, name of the assay with the fits.

- averaged_assay:

  string, name of the assay used to count the concentrations (`N_conc`).

## Value

`se` with the PRISM curves in `metrics_assay`. The fit is recorded in
the provenance as fit source `"PRISM"`, and the step
`"prism_drc_import"` holds the input files of `prism_curves` (name and
md5) and the number of curves imported and left out. Re-running the
recipe needs the same `prism_curves`.

## Details

PRISM fits
`upper_limit + (lower_limit - upper_limit) / (1 + (dose / inflection)^slope)`,
which is the gDR curve with `x_0 = upper_limit`, `x_inf = lower_limit`,
`ec50 = inflection` and `h = -slope`. `r2` is `frac_var_explained` and
`fit_type` is PRISM's `fit_name`. PRISM defines AUC and IC50 differently
from gDR (the curve is clamped between 0 and 1, an IC50 outside the
tested range is `NA`), so those values are kept under their own names
(`prism_auc`, `prism_ic50` and so on) and the gDR columns `x_AOC`,
`xc50`, `x_mean` and `x_max` stay `NA`.

A PRISM screen often has several curves for one cell line and compound
(different pools, barcodes or compound plates), which gDR averages into
one cell. Each of them becomes its own row, with `depmap_id`, `pert_id`,
`pool_id`, `lua`, `cell_set` and `pert_plate` to tell them apart. Curves
with no matching cell in `se` are left out with a warning. Running the
function again replaces the PRISM rows of the earlier run.

Run it after
[`fit_SE`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md),
which replaces the whole `Metrics` assay.

## Examples

``` r
if (requireNamespace("gDRimport", quietly = TRUE) &&
    exists("read_PRISM_drc", asNamespace("gDRimport"))) {
  mae <- gDRutils::get_synthetic_data("finalMAE_prism")
  se <- MultiAssayExperiment::experiments(mae)[["single-agent"]]
  curves <- gDRimport::read_PRISM_drc(
    drc_path = system.file("testdata/prism_drc.csv", package = "gDRimport"),
    meta_data_path = system.file("testdata/prism_drc_model.csv", package = "gDRimport")
  )
  se <- add_PRISM_metrics(se, curves)
  metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  metrics[metrics$fit_source == "PRISM", ]
}
#>                                            rId
#>                                         <char>
#> 1: G00106_GDC-8025_TEAD1|TEAD2|TEAD3|TEAD4_120
#> 2: G00106_GDC-8025_TEAD1|TEAD2|TEAD3|TEAD4_120
#> 3: G00106_GDC-8025_TEAD1|TEAD2|TEAD3|TEAD4_120
#> 4: G00106_GDC-8025_TEAD1|TEAD2|TEAD3|TEAD4_120
#>                                                                                                                   cId
#>                                                                                                                <char>
#> 1:           HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_Myeloid_unknown_unknown_NA
#> 2:   HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_Myeloid_unknown_unknown_NA
#> 3:   HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_Myeloid_unknown_unknown_NA
#> 4: MONOMAC6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_MONOMAC6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_Myeloid_unknown_unknown_NA
#>    x_mean x_AOC x_AOC_range  xc50 x_max  ec50 x_inf   x_0     h    r2 p_value
#>     <num> <num>       <num> <num> <num> <num> <num> <num> <num> <num>   <num>
#> 1:     NA    NA          NA    NA    NA  0.85  0.21  0.98   1.4  0.91      NA
#> 2:     NA    NA          NA    NA    NA  2.40  0.35  1.01   0.9  0.77      NA
#> 3:     NA    NA          NA    NA    NA  1.70  0.29  0.96   1.1  0.83      NA
#> 4:     NA    NA          NA    NA    NA    NA    NA    NA    NA    NA      NA
#>      rss maxlog10Concentration N_conc x_sd_avg                fit_type
#>    <num>                 <num>  <int>    <num>                  <char>
#> 1:    NA                     1      8       NA     drc_drm_constrained
#> 2:    NA                     1      8       NA dr4pl_initL_constrained
#> 3:    NA                     1      8       NA   drc_drm_unconstrained
#> 4:    NA                     1      8       NA                    <NA>
#>    normalization_type fit_source prism_auc prism_log2_auc prism_auc_riemann
#>                <char>     <char>     <num>          <num>             <num>
#> 1:                 RV      PRISM    0.6012        -0.7341              0.58
#> 2:                 RV      PRISM    0.7208        -0.4723              0.71
#> 3:                 RV      PRISM    0.6705        -0.5767              0.66
#> 4:                 RV      PRISM        NA             NA              0.97
#>    prism_ic50 prism_log2_ic50 prism_mse prism_mad prism_successful_fit
#>         <num>           <num>     <num>     <num>               <lgcl>
#> 1:     1.0750          0.1043     0.011      0.08                 TRUE
#> 2:         NA              NA     0.014      0.09                 TRUE
#> 3:     1.9378          0.9543     0.012      0.07                 TRUE
#> 4:         NA              NA        NA        NA                FALSE
#>    prism_minimum_dose prism_maximum_dose  depmap_id              pert_id
#>                 <num>              <int>     <char>               <char>
#> 1:            0.00046                 10 ACH-000004 PRC-000000001-001-01
#> 2:            0.00046                 10 ACH-000005 PRC-000000001-001-01
#> 3:            0.00046                 10 ACH-000005 PRC-000000001-001-01
#> 4:            0.00046                 10 ACH-000006 PRC-000000001-001-01
#>    pool_id      lua cell_set pert_plate x_project_id      screen Gnumber
#>     <char>   <char>   <char>     <char>       <char>      <char>  <char>
#> 1:      P1 LUA-0001      CS1     PLATE1 TEST_PROJECT TEST_SCREEN  G00106
#> 2:      P1 LUA-0002      CS1     PLATE1 TEST_PROJECT TEST_SCREEN  G00106
#> 3:      P2 LUA-0003      CS1     PLATE2 TEST_PROJECT TEST_SCREEN  G00106
#> 4:      P1 LUA-0004      CS1     PLATE1 TEST_PROJECT TEST_SCREEN  G00106
#>    DrugName                drug_moa Duration
#>      <char>                  <char>    <num>
#> 1: GDC-8025 TEAD1|TEAD2|TEAD3|TEAD4      120
#> 2: GDC-8025 TEAD1|TEAD2|TEAD3|TEAD4      120
#> 3: GDC-8025 TEAD1|TEAD2|TEAD3|TEAD4      120
#> 4: GDC-8025 TEAD1|TEAD2|TEAD3|TEAD4      120
#>                                           clid
#>                                         <char>
#> 1:      HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
#> 2:  HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
#> 3:  HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
#> 4: MONOMAC6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
#>                                   CellLineName  Tissue parental_identifier
#>                                         <char>  <char>              <char>
#> 1:      HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE Myeloid             unknown
#> 2:  HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE Myeloid             unknown
#> 3:  HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE Myeloid             unknown
#> 4: MONOMAC6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE Myeloid             unknown
#>    subtype ReferenceDivisionTime
#>     <char>                <lgcl>
#> 1: unknown                    NA
#> 2: unknown                    NA
#> 3: unknown                    NA
#> 4: unknown                    NA
```
