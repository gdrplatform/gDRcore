# Get the count of untreated tags per row

Get the count of untreated tags per row

## Usage

``` r
.get_untreated_tag_count(
  mat_elem,
  drug_identifier_keys = c("drug_name", "drug_name2", "drug_name3")
)
```

## Arguments

- mat_elem:

  data.table input data frame to evaluate.

- drug_identifier_keys:

  character vector of keys used to look up drug column names in the
  `gDRutils` environment.

## Value

list containing ntag, num_cols, and valid_cols
