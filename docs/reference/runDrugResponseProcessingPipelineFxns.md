# Run drug response processing pipeline

Run different components of the gDR drug response processing pipeline.
Either: create a SummarizedExperiment and normalize raw treated and
control data (create_and_normalize_SE), average data (average_SE), or
fit the processed data (fit_SE). See details for more in-depth
explanations.

## Usage

``` r
average_SE(
  se,
  data_type,
  series_identifiers = NULL,
  normalized_assay = "Normalized",
  averaged_assay = "Averaged"
)

create_SE(
  df_,
  data_type,
  readout = "ReadoutValue",
  nested_identifiers = NULL,
  nested_confounders = intersect(names(df_), gDRutils::get_env_identifiers("barcode")),
  override_untrt_controls = NULL
)

fit_SE(
  se,
  data_type = "single-agent",
  nested_identifiers = NULL,
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  n_point_cutoff = 4,
  range_conc = c(0.005, 5),
  force_fit = FALSE,
  pcutoff = 0.05,
  cap = 0.1,
  curve_type = c("GR", "RV")
)

normalize_SE(
  se,
  data_type,
  nested_identifiers = NULL,
  nested_confounders = gDRutils::get_SE_identifiers(se, "barcode", simplify = TRUE),
  control_mean_fxn = function(x) {
     mean(x, trim = 0.25)
 },
  control_assay = "Controls",
  raw_treated_assay = "RawTreated",
  normalized_assay = "Normalized",
  ndigit_rounding = 4
)

create_and_normalize_SE(
  df_,
  data_type,
  readout = "ReadoutValue",
  control_mean_fxn = function(x) {
     mean(x, trim = 0.25)
 },
  nested_identifiers = NULL,
  nested_confounders = intersect(names(df_), gDRutils::get_env_identifiers("barcode")),
  override_untrt_controls = NULL,
  ndigit_rounding = 4,
  control_assay = "Controls",
  raw_treated_assay = "RawTreated",
  normalized_assay = "Normalized"
)

runDrugResponseProcessingPipeline(
  x,
  readout = "ReadoutValue",
  control_mean_fxn = function(x) {
     mean(x, trim = 0.25)
 },
  nested_identifiers_l = NULL,
  nested_confounders = gDRutils::get_env_identifiers("barcode"),
  override_untrt_controls = NULL,
  ndigit_rounding = 4,
  n_point_cutoff = 4,
  control_assay = "Controls",
  raw_treated_assay = "RawTreated",
  normalized_assay = "Normalized",
  averaged_assay = "Averaged",
  metrics_assay = "Metrics",
  split_data = TRUE,
  data_dir = NULL,
  partial_run = FALSE,
  start_from = get_pipeline_steps()[1],
  selected_experiments = NULL
)
```

## Arguments

- se:

  `SummarizedExperiment` object.

- data_type:

  single-agent vs combination

- series_identifiers:

  character vector of identifiers in `measured` or `metric` which define
  a unique data point.

- normalized_assay:

  string of the assay name containing the normalized data. Defaults to
  `"Normalized"`.

- averaged_assay:

  string of the name of the averaged assay in the SummarizedExperiment.
  Defaults to `"Averaged"`.

- df\_:

  data.table of raw drug response data containing both treated and
  untreated values. If a column called `"BackgroundValue"` exists in
  `df_`, it will be removed from the `readout` column.

- readout:

  string of the name containing the cell viability readout values.

- nested_identifiers:

  character vector with the nested_identifiers for the given SE with a
  given data_type

- nested_confounders:

  Character vector of the nested_confounders for a given assay.
  nested_keys is character vector of column names to include in the
  data.tables in the assays of the resulting `SummarizedExperiment`
  object. Defaults to the `nested_identifiers` and `nested_confounders`
  if passed through `create_and_normalize_SE` or
  `runDrugResponseProcessingPipeline`.

- override_untrt_controls:

  named list containing defining factors in the treatments. Defaults to
  `NULL`.

- metrics_assay:

  string of the name of the metrics assay to output in the returned
  SummarizedExperiment Defaults to `"Metrics"`.

- n_point_cutoff:

  integer of how many points should be considered the minimum required
  to try to fit a curve. Defaults to `4`.

- range_conc:

  vector of concetrations range values.

- force_fit:

  boolean indicating whether or not to force the fit.

- pcutoff:

  numeric cutoff value.

- cap:

  numeric value representing the value to cap the highest allowed
  relative viability at.

- curve_type:

  vector of curve type values.

- control_mean_fxn:

  function indicating how to average controls. Defaults to
  `mean(x, trim = 0.25)`.

- control_assay:

  string containing the name of the assay representing the controls in
  the `se`. Defaults to `"Controls"`.

- raw_treated_assay:

  string containing the name of the assay representing the raw treated
  data in the `se`. Defaults to `"RawTreated"`.

- ndigit_rounding:

  integer indicating number of digits to round to in calculations.
  Defaults to `4`.

- x:

  data.table of MAE with drug response data

- nested_identifiers_l:

  list with the nested_identifiers(character v ectors) for
  `single-agent` and (optionally) for `combination` data

- split_data:

  boolean indicating whether data provided as the MultiAssayExperiment
  should be split again into appropriate data types

- data_dir:

  string with the path to the directory with intermediate data of
  experiments (qs2 files). If set to NULL (default) intermediate data is
  not saved/read in.

- partial_run:

  logical flag indicating if the pipeline should be run partially (from
  the step defined with `start_from`)

- start_from:

  string indicating the pipeline step from which partial run should be
  launched

- selected_experiments:

  character vector with experiments for which pipeline should be run.
  This option works only for the pipeline being run partially (i.e. with
  `partial_run` flag set to `TRUE`)

## Value

MAE object

## Details

`runDrugResponseProcessingPipeline` is made up of 3 separate steps:

- "create_and_normalize_SE"

- "average_SE"

- "fit_SE"

For create_and_normalize_SE, this creates a SummarizedExperiment object
from a data.table, where the data.table contains treatments on rows, and
conditions on columns. A SummarizedExperiment object containing two
asssays is created: treated readouts will live in an assay called
`"RawTreated"`, and reference readouts live in an assay called
`"Controls"`. Subsequently, the treated and control elements will be
normalized to output two metrics:

For average_SE, take the normalized assay and average the nested
`DataFrame`s across unique`nested_identifiers`.

For fit_SE, take the averaged assay and fit curves to obtain metrics,
one set of metrics for each normalization type set.

Pipeline can be run partially with `partial_run` flag set to TRUE. The
`start_from` string defines the step from which the pipeline will be
launched. However, partial run of the pipeline is possible only if the
whole pipeline was launched at least once with defined `data_dir` and
intermediate data was saved as qs2 files into `data_dir`.

Pipeline can be run for the selected experiments by changing the default
value of `selected_experiments` param. This scenario only works when
`partial_run` is enabled.

## Examples

``` r
d <- rep(seq(0.1, 0.9, 0.1), each = 4)
v <- rep(seq(0.1, 0.4, 0.1), 9)
df <- S4Vectors::DataFrame(
  Concentration = d,
  normalization_type = rep(c("GR", "RV"), length(v) * 2),
  x = rep(v, 2)
)
normalized <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = df)

keys <- list(Trt = "Concentration")
assays <- list("Normalized" = normalized)
se <- SummarizedExperiment::SummarizedExperiment(assays = assays)
se <- gDRutils::set_SE_keys(se, keys)
se <- gDRutils::set_SE_identifiers(se, gDRutils::get_env_identifiers())
se1 <- average_SE(
  se,
  data_type = "single-agent",
  normalized_assay = "Normalized",
  averaged_assay = "Averaged"
)
#> Loading required namespace: testthat


td <- gDRimport::get_test_data()
l_tbl <- gDRimport::load_data(
  manifest_file = gDRimport::manifest_path(td),
  df_template_files = gDRimport::template_path(td),
  results_file = gDRimport::result_path(td)
)
#> INFO [2026-06-02 11:39:16] Manifest loaded successfully
#> INFO [2026-06-02 11:39:16] Reading Template_7daytreated.xlsx with load_templates_xlsx
#> INFO [2026-06-02 11:39:16] Reading Template_Untreated.xlsx with load_templates_xlsx
#> INFO [2026-06-02 11:39:16] Loading Template_7daytreated.xlsx
#> INFO [2026-06-02 11:39:16] Loading Template_Untreated.xlsx
#> INFO [2026-06-02 11:39:16] Templates loaded successfully!
#> INFO [2026-06-02 11:39:16] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day0.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-02 11:39:16] Plate 201904190a read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904190b read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904190c read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904190d read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904190e read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904190f read; 384 wells
#> INFO [2026-06-02 11:39:16] File done
#> INFO [2026-06-02 11:39:16] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day7.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-02 11:39:16] Plate 201904197a read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904197b read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904197c read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904197d read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904197e read; 384 wells
#> INFO [2026-06-02 11:39:16] Plate 201904197f read; 384 wells
#> INFO [2026-06-02 11:39:16] File done
imported_data <- merge_data(
  l_tbl$manifest,
  l_tbl$treatments,
  l_tbl$data
)
#> INFO [2026-06-02 11:39:16] Merging data
#> INFO [2026-06-02 11:39:16] Merging the metadata (manifest and treatment files)
#> WARN [2026-06-02 11:39:16] 4608 well loaded, 768 wells discarded for lack of annotation,
#>     3840 data point selected

se <- purrr::quietly(create_SE)(imported_data, data_type = "single-agent")


td <- gDRimport::get_test_data()
l_tbl <- gDRimport::load_data(
  manifest_file = gDRimport::manifest_path(td),
  df_template_files = gDRimport::template_path(td),
  results_file = gDRimport::result_path(td)
)
#> INFO [2026-06-02 11:39:17] Manifest loaded successfully
#> INFO [2026-06-02 11:39:17] Reading Template_7daytreated.xlsx with load_templates_xlsx
#> INFO [2026-06-02 11:39:17] Reading Template_Untreated.xlsx with load_templates_xlsx
#> INFO [2026-06-02 11:39:17] Loading Template_7daytreated.xlsx
#> INFO [2026-06-02 11:39:17] Loading Template_Untreated.xlsx
#> INFO [2026-06-02 11:39:17] Templates loaded successfully!
#> INFO [2026-06-02 11:39:17] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day0.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-02 11:39:17] Plate 201904190a read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904190b read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904190c read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904190d read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904190e read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904190f read; 384 wells
#> INFO [2026-06-02 11:39:17] File done
#> INFO [2026-06-02 11:39:17] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day7.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-02 11:39:17] Plate 201904197a read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904197b read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904197c read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904197d read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904197e read; 384 wells
#> INFO [2026-06-02 11:39:17] Plate 201904197f read; 384 wells
#> INFO [2026-06-02 11:39:17] File done
imported_data <- merge_data(
  l_tbl$manifest,
  l_tbl$treatments,
  l_tbl$data
)
#> INFO [2026-06-02 11:39:17] Merging data
#> INFO [2026-06-02 11:39:17] Merging the metadata (manifest and treatment files)
#> WARN [2026-06-02 11:39:17] 4608 well loaded, 768 wells discarded for lack of annotation,
#>     3840 data point selected

inl <- prepare_input(imported_data)
#> Warning: 'Plate' nested confounder(s) is/are not present in the data.
#>     Switching into 'Barcode' nested confounder(s).
se <- create_SE(
 inl$df_list[["single-agent"]],
 data_type = "single-agent",
 nested_confounders = inl$nested_confounders)

normalize_SE(se, data_type = "single-agent")
#> class: SummarizedExperiment 
#> dim: 3 6 
#> metadata(3): identifiers experiment_metadata Keys
#> assays(3): RawTreated Controls Normalized
#> rownames(3): G00002_drug_002_moa_A_168 G00004_drug_004_moa_A_168
#>   G00011_drug_011_moa_B_168
#> rowData names(4): Gnumber DrugName drug_moa Duration
#> colnames(6): CL00011_cellline_BA_breast_cellline_BA_unknown_26
#>   CL00012_cellline_CA_breast_cellline_CA_unknown_30 ...
#>   CL00015_cellline_FA_breast_cellline_FA_unknown_42
#>   CL00018_cellline_IB_breast_cellline_IB_unknown_54
#> colData names(6): clid CellLineName ... subtype ReferenceDivisionTime
p_dir <- file.path(tempdir(), "pcheck")
dir.create(p_dir)
td <- gDRimport::get_test_data()
l_tbl <- gDRimport::load_data(
  manifest_file = gDRimport::manifest_path(td),
  df_template_files = gDRimport::template_path(td),
  results_file = gDRimport::result_path(td)
)
#> INFO [2026-06-02 11:39:18] Manifest loaded successfully
#> INFO [2026-06-02 11:39:18] Reading Template_7daytreated.xlsx with load_templates_xlsx
#> INFO [2026-06-02 11:39:18] Reading Template_Untreated.xlsx with load_templates_xlsx
#> INFO [2026-06-02 11:39:18] Loading Template_7daytreated.xlsx
#> INFO [2026-06-02 11:39:18] Loading Template_Untreated.xlsx
#> INFO [2026-06-02 11:39:18] Templates loaded successfully!
#> INFO [2026-06-02 11:39:18] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day0.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-02 11:39:19] Plate 201904190a read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904190b read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904190c read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904190d read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904190e read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904190f read; 384 wells
#> INFO [2026-06-02 11:39:19] File done
#> INFO [2026-06-02 11:39:19] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day7.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-02 11:39:19] Plate 201904197a read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904197b read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904197c read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904197d read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904197e read; 384 wells
#> INFO [2026-06-02 11:39:19] Plate 201904197f read; 384 wells
#> INFO [2026-06-02 11:39:19] File done
imported_data <- merge_data(
  l_tbl$manifest,
  l_tbl$treatments,
  l_tbl$data
)
#> INFO [2026-06-02 11:39:19] Merging data
#> INFO [2026-06-02 11:39:19] Merging the metadata (manifest and treatment files)
#> WARN [2026-06-02 11:39:19] 4608 well loaded, 768 wells discarded for lack of annotation,
#>     3840 data point selected
runDrugResponseProcessingPipeline(
  imported_data,
  data_dir = p_dir
)
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
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> overriding original x_0 argument '1' with '1.08355555555556' (fit is not statistically significant (p=1.00), setting constant fit)
#> overriding original x_0 argument '1' with '1.1' (only 1 normalized value detected, setting constant fit)
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> overriding original x_0 argument '1' with '1.09306666666667' (fit is not statistically significant (p=1.00), setting constant fit)
#> overriding original x_0 argument '1' with '1.1' (only 1 normalized value detected, setting constant fit)
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> NaNs produced
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
#> not enough data points (1 < 4) to perform fitting
#> not enough data points (1 < 4) to perform fitting
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] combination: SummarizedExperiment with 2 rows and 6 columns
#>  [2] single-agent: SummarizedExperiment with 3 rows and 6 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
