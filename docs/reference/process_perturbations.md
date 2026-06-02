# Cleanup additional perturbations in the data.table

This function processes drug and concentration columns in a data.table.
It checks if there is only one unique drug (excluding a specified
untreated tag) and if there are exactly two doses (one of which is 0).
If these conditions are met, it creates a new column named after the
drug and fills it with the doses, then removes the original drug and
concentration columns.

## Usage

``` r
process_perturbations(
  dt,
  drugs_cotrt_ids,
  conc_cotrt_ids,
  untreated_tag = "vehicle"
)
```

## Arguments

- dt:

  A data.table containing the data.

- drugs_cotrt_ids:

  A vector of column names related to drugs.

- conc_cotrt_ids:

  A vector of column names related to concentrations.

- untreated_tag:

  A string representing the untreated tag (default is "vehicle").

## Value

A modified data.table with new columns for the drugs and removed
original drug and concentration columns.

## Examples

``` r
dt <- data.table::data.table(
  drug1 = c("vehicle", "drugA", "drugA"),
  conc1 = c(0, 10, 0),
  drug2 = c("vehicle", "drugB", "drugB"),
  conc2 = c(0, 20, 0)
)
drugs_cotrt_ids <- c("drug1", "drug2")
conc_cotrt_ids <- c("conc1", "conc2")
dt <- process_perturbations(dt, drugs_cotrt_ids, conc_cotrt_ids)
print(dt)
#>    drugA drugB
#>    <num> <num>
#> 1:     0     0
#> 2:    10    20
#> 3:     0     0
```
