# Prepare input data common for all experiments

Current steps

- refining nested confounders

- refining nested identifiers

- splitting df\_ into (per experiment) df_list

## Usage

``` r
# S3 method for class 'MultiAssayExperiment'
prepare_input(
  x,
  nested_confounders = gDRutils::get_SE_identifiers(x[[1]], "barcode"),
  nested_identifiers_l = .get_default_nested_identifiers(x[[1]]),
  raw_data_field = "experiment_raw_data",
  split_data = TRUE,
  ...
)
```

## Arguments

- x:

  MAE object with dose-response data

- nested_confounders:

  Character vector of the nested_confounders for a given assay.
  nested_keys is character vector of column names to include in the
  data.tables in the assays of the resulting `SummarizedExperiment`
  object. Defaults to the `nested_identifiers` and `nested_confounders`
  if passed through

- nested_identifiers_l:

  list with the nested_identifiers(character vectors) for `single-agent`
  and (optionally) for `combination` data

- raw_data_field:

  metadata field with raw data

- split_data:

  Boolean indicating need of splitting the data into experiment types

- ...:

  additional parameters

## Value

list of input data
