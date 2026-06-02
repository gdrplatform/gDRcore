# identify_keys

Group columns from a data.table that correspond to different

## Usage

``` r
identify_keys(
  df_,
  nested_keys = NULL,
  override_untrt_controls = NULL,
  identifiers = gDRutils::get_env_identifiers()
)
```

## Arguments

- df\_:

  a data.table to identify keys for.

- nested_keys:

  character vector of keys to exclude from the returned list. The keys
  discarded should be identical to the keys in the third dimension of
  the SummarizedExperiment. Defaults to the `"Barcode"` and the `masked`
  identifier.

- override_untrt_controls:

  named list containing defining factors in the treatments. Defaults to
  `NULL`.

- identifiers:

  named list containing all identifiers to use during processing. By
  default, this value will be obtained by the environment.

## Value

named list of key types and their corresponding key values.

## Details

This is most likely to be used for provenance tracking and will be
placed on the SummarizedExperiment metadata for downstream analyses to
reference.

## See also

map_df, create_SE

## Examples

``` r
n <- 64
md_df <- data.table::data.table(
  Gnumber = rep(c("vehicle", "untreated", paste0("G", seq(2))), each = 16),
  DrugName = rep(c("vehicle", "untreated", paste0("GN", seq(2))), each = 16),
  clid = paste0("C", rep_len(seq(4), n)),
  CellLineName = paste0("N", rep_len(seq(4), n)),
  replicates = rep_len(paste0("R", rep(seq(4), each = 4)), 64),
  drug_moa = "inhibitor",
  ReferenceDivisionTime = rep_len(c(120, 60), n),
  Tissue = "Lung",
  parental_identifier = "CL12345",
  Duration = 160
)
md_df <- unique(md_df)
ref <- md_df$Gnumber %in% c("vehicle", "untreated")
trt_df <- md_df[!ref, ]
identify_keys(trt_df)
#> $Trt
#> [1] "Gnumber"             "DrugName"            "clid"               
#> [4] "CellLineName"        "replicates"          "drug_moa"           
#> [7] "Tissue"              "parental_identifier" "Duration"           
#> 
#> $ref_Endpoint
#> [1] "clid"                "CellLineName"        "replicates"         
#> [4] "Tissue"              "parental_identifier" "Duration"           
#> 
#> $untrt_Endpoint
#> [1] "clid"                "CellLineName"        "replicates"         
#> [4] "Tissue"              "parental_identifier" "Duration"           
#> 
#> $Day0
#> [1] "clid"                "CellLineName"        "replicates"         
#> [4] "Tissue"              "parental_identifier"
#> 
#> $nested_keys
#> NULL
#> 
#> $masked_tag
#> [1] "masked"
#> 
#> $cellline_name
#> [1] "CellLineName"
#> 
#> $cellline_ref_div_time
#> [1] "ReferenceDivisionTime"
#> 
#> $duration
#> [1] "Duration"
#> 
#> $untreated_tag
#> [1] "vehicle"   "untreated"
#> 
```
