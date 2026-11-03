# Read the SummarizedExperiment saved after a pipeline step

Run with `data_dir`,
[`runDrugResponseProcessingPipeline`](https://gdrplatform.github.io/gDRcore/reference/runDrugResponseProcessingPipelineFxns.md)
saves the SummarizedExperiment of every experiment after every step as
`<data_dir>/<experiment>__<step>.qs2`. This reads one of them back.

## Usage

``` r
read_intermediate_data(path, step, experiment)
```

## Arguments

- path:

  string with the directory passed as `data_dir`

- step:

  string with the step name: `"create_SE"`, `"normalize_SE"`,
  `"average_SE"` or `"fit_SE"`

- experiment:

  string with the experiment name, e.g. `"single-agent"`

## Value

SummarizedExperiment

## Examples

``` r
p_dir <- file.path(tempdir(), "read_step")
dir.create(p_dir)
td <- gDRimport::get_test_data()
l_tbl <- gDRimport::load_data(
  manifest_file = gDRimport::manifest_path(td),
  df_template_files = gDRimport::template_path(td),
  results_file = gDRimport::result_path(td)
)
#> INFO [2026-11-03 14:10:32] Manifest loaded successfully
#> INFO [2026-11-03 14:10:32] Reading Template_7daytreated.xlsx with load_templates_xlsx
#> INFO [2026-11-03 14:10:32] Reading Template_Untreated.xlsx with load_templates_xlsx
#> INFO [2026-11-03 14:10:32] Loading Template_7daytreated.xlsx
#> INFO [2026-11-03 14:10:32] Loading Template_Untreated.xlsx
#> INFO [2026-11-03 14:10:32] Templates loaded successfully!
#> INFO [2026-11-03 14:10:32] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day0.xlsx, sheet Readout_0077vs0068_day7
#> New names:
#> • `` -> `...1`
#> • `` -> `...2`
#> • `` -> `...3`
#> • `` -> `...4`
#> • `` -> `...5`
#> • `` -> `...6`
#> • `` -> `...7`
#> • `` -> `...8`
#> • `` -> `...9`
#> • `` -> `...10`
#> • `` -> `...11`
#> • `` -> `...12`
#> • `` -> `...13`
#> • `` -> `...14`
#> • `` -> `...15`
#> • `` -> `...16`
#> • `` -> `...17`
#> • `` -> `...18`
#> • `` -> `...19`
#> • `` -> `...20`
#> • `` -> `...21`
#> • `` -> `...22`
#> • `` -> `...23`
#> • `` -> `...24`
#> • `` -> `...25`
#> INFO [2026-11-03 14:10:32] Plate 201904190a read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904190b read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904190c read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904190d read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904190e read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904190f read; 384 wells
#> INFO [2026-11-03 14:10:32] File done
#> INFO [2026-11-03 14:10:32] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day7.xlsx, sheet Readout_0077vs0068_day7
#> New names:
#> • `` -> `...1`
#> • `` -> `...2`
#> • `` -> `...3`
#> • `` -> `...4`
#> • `` -> `...5`
#> • `` -> `...6`
#> • `` -> `...7`
#> • `` -> `...8`
#> • `` -> `...9`
#> • `` -> `...10`
#> • `` -> `...11`
#> • `` -> `...12`
#> • `` -> `...13`
#> • `` -> `...14`
#> • `` -> `...15`
#> • `` -> `...16`
#> • `` -> `...17`
#> • `` -> `...18`
#> • `` -> `...19`
#> • `` -> `...20`
#> • `` -> `...21`
#> • `` -> `...22`
#> • `` -> `...23`
#> • `` -> `...24`
#> • `` -> `...25`
#> INFO [2026-11-03 14:10:32] Plate 201904197a read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904197b read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904197c read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904197d read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904197e read; 384 wells
#> INFO [2026-11-03 14:10:32] Plate 201904197f read; 384 wells
#> INFO [2026-11-03 14:10:32] File done
imported_data <- merge_data(l_tbl$manifest, l_tbl$treatments, l_tbl$data)
#> INFO [2026-11-03 14:10:32] Merging data
#> INFO [2026-11-03 14:10:32] Merging the metadata (manifest and treatment files)
#> WARN [2026-11-03 14:10:32] 4608 well loaded, 768 wells discarded for lack of annotation,
#>     3840 data point selected
mae <- runDrugResponseProcessingPipeline(imported_data, data_dir = p_dir)
#> Warning: 'Plate' nested confounder(s) is/are not present in the data.
#>     Switching into 'Barcode' nested confounder(s).
#> Processing combination
#> Warning: mapping original concentration '0.00457247142398638' to '0.00437'
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> mapping original concentration '0.00457247142398638' to '0.00437'
#> not enough data points (1 < 4) to perform fitting
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> mapping original concentration '0.00457247142398638' to '0.00437'
#> not enough data points (1 < 4) to perform fitting
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> mapping original concentration '0.00457247142398638' to '0.00437'
#> not enough data points (1 < 4) to perform fitting
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> mapping original concentration '0.00457247142398638' to '0.00437'
#> not enough data points (1 < 4) to perform fitting
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitt
#> Processing single-agent
#> Warning: method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
read_intermediate_data(p_dir, "normalize_SE", "single-agent")
#> class: SummarizedExperiment 
#> dim: 3 6 
#> metadata(4): identifiers experiment_metadata Keys provenance
#> assays(3): RawTreated Controls Normalized
#> rownames(3): G00002_drug_002_moa_A_168 G00004_drug_004_moa_A_168
#>   G00011_drug_011_moa_B_168
#> rowData names(4): Gnumber DrugName drug_moa Duration
#> colnames(6): CL00011_cellline_BA_breast_cellline_BA_unknown_26
#>   CL00012_cellline_CA_breast_cellline_CA_unknown_30 ...
#>   CL00015_cellline_FA_breast_cellline_FA_unknown_42
#>   CL00018_cellline_IB_breast_cellline_IB_unknown_54
#> colData names(6): clid CellLineName ... subtype ReferenceDivisionTime
```
