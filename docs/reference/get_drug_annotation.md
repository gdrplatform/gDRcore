# get_drug_annotation

Get drug annotation data table

## Usage

``` r
get_drug_annotation(
  data,
  fname = "drugs.csv",
  fill = "unknown",
  annotation_package = if ("gDRinternal" %in% .packages(all.available = TRUE)) {
    
    "gDRinternal"
 } else {
     "gDRtestData"
 }
)
```

## Arguments

- data:

  data.table with drug identifiers to be matched

- fname:

  string with file name containing the annotation

- fill:

  string indicating how unknown drugs should be filled in the DB

- annotation_package:

  string indicating name of the package containing drug annotation

## Value

data.table with drug annotations

## Examples

``` r
data <- data.table::data.table(Gnumber = c("drug1", "drug2", "drug3"))
drug_annotation <- get_drug_annotation(data)
```
