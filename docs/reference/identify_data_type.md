# Identify type of data

Identify type of data

## Usage

``` r
identify_data_type(dt, codilution_conc = 2, matrix_conc = 1)
```

## Arguments

- dt:

  data.table of raw drug response data containing both treated and
  untreated values

- codilution_conc:

  integer of maximum number of concentration ratio of co-treatment to
  classify as codilution data type; defaults to `2`

- matrix_conc:

  integer of minimum number of concentration pairs of co-treatment to
  classify as co-treatment or matrix data type; defaults to `1`

## Value

data.table of raw drug response data with additional column `type` with
the info of data type for a given row of data.table

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
conc <- rep(seq(0, 0.3, 0.1), 2)
ctrl_dt <- S4Vectors::DataFrame(
  ReadoutValue = c(2, 2, 1, 1, 2, 1),
  Concentration = rep(0, 6),
  masked = FALSE,
  DrugName = rep(c("DRUG_10", "vehicle", "DRUG_8"), 2),
  CellLineName = "CELL1"
)

trt_dt <- S4Vectors::DataFrame(
  ReadoutValue = rep(seq(1, 4, 1), 2),
  Concentration = conc,
  masked = rep(FALSE, 8),
  DrugName = c("DRUG_10", "DRUG_8"),
  CellLineName = "CELL1"
)
input_dt <- data.table::as.data.table(rbind(ctrl_dt, trt_dt))
input_dt$Duration <- 72
input_dt$CorrectedReadout2 <- input_dt$ReadoutValue
identify_data_type(input_dt)
#>     ReadoutValue Concentration masked DrugName CellLineName Duration
#>            <num>         <num> <lgcl>   <char>       <char>    <num>
#>  1:            2           0.0  FALSE  DRUG_10        CELL1       72
#>  2:            2           0.0  FALSE  vehicle        CELL1       72
#>  3:            1           0.0  FALSE   DRUG_8        CELL1       72
#>  4:            1           0.0  FALSE  DRUG_10        CELL1       72
#>  5:            2           0.0  FALSE  vehicle        CELL1       72
#>  6:            1           0.0  FALSE   DRUG_8        CELL1       72
#>  7:            1           0.0  FALSE  DRUG_10        CELL1       72
#>  8:            2           0.1  FALSE   DRUG_8        CELL1       72
#>  9:            3           0.2  FALSE  DRUG_10        CELL1       72
#> 10:            4           0.3  FALSE   DRUG_8        CELL1       72
#> 11:            1           0.0  FALSE  DRUG_10        CELL1       72
#> 12:            2           0.1  FALSE   DRUG_8        CELL1       72
#> 13:            3           0.2  FALSE  DRUG_10        CELL1       72
#> 14:            4           0.3  FALSE   DRUG_8        CELL1       72
#>     CorrectedReadout2 record_id         type
#>                 <num>     <int>       <char>
#>  1:                 2         1      control
#>  2:                 2         2      control
#>  3:                 1         3      control
#>  4:                 1         4      control
#>  5:                 2         5      control
#>  6:                 1         6      control
#>  7:                 1         7      control
#>  8:                 2         8 single-agent
#>  9:                 3         9 single-agent
#> 10:                 4        10 single-agent
#> 11:                 1        11      control
#> 12:                 2        12 single-agent
#> 13:                 3        13 single-agent
#> 14:                 4        14 single-agent
```
