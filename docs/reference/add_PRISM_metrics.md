# Add the curve fits delivered by PRISM to the Metrics assay

Adds the curves PRISM fitted (read with `gDRimport::read_PRISM_drc()`)
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
  `gDRimport::read_PRISM_drc()`: the columns of the DRC table, the cell
  line identifier and Duration.

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
```
