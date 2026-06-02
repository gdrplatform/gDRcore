# Split raw data into list based on the data types

Split raw data into list based on the data types

## Usage

``` r
split_raw_data(dt, type_col = "type")
```

## Arguments

- dt:

  data.table of raw drug response data containing both treated and
  untreated values with column specified in `type_col` argument.

- type_col:

  string with column names in `dt` with info about data type. Defaults
  to `"type"`.

## Value

list with split data based on its data type

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
cell_lines <- gDRtestData::create_synthetic_cell_lines()
drugs <- gDRtestData::create_synthetic_drugs()
dt_layout <- drugs[4:6, as.list(cell_lines[7:8, ]), names(drugs)]
dt_layout <- gDRtestData::add_data_replicates(dt_layout)
dt_layout <- gDRtestData::add_concentration(
  dt_layout,
  concentrations = 10 ^ (seq(-3, .5, .5))
)

dt_2 <-
  drugs[c(21, 26), as.list(cell_lines[which(cell_lines$clid %in% dt_layout$clid)]), names(drugs)]
dt_2 <- gDRtestData::add_data_replicates(dt_2)
dt_2 <- gDRtestData::add_concentration(
  dt_2,
  concentrations = 10 ^ (seq(-3, .5, .5))
)
colnames(dt_2)[colnames(dt_2) %in% c(colnames(drugs), "Concentration")] <-
  paste0(
    colnames(dt_2)[colnames(dt_2) %in% c(colnames(drugs), "Concentration")],
    "_2"
  )
dt_layout_2 <- dt_layout[dt_2, on = intersect(names(dt_layout), names(dt_2)),
                        allow.cartesian = TRUE]
dt_merged_data <- gDRtestData::generate_response_data(dt_layout_2, 0)
dt <- identify_data_type(dt_merged_data)
split_raw_data(dt)
#> $combination
#>       Barcode Gnumber DrugName drug_moa    clid CellLineName   Tissue
#>        <char>  <char>   <char>   <char>  <char>       <char>   <char>
#>    1: plate_1  G00004 drug_004    moa_A CL00016  cellline_GB tissue_y
#>    2: plate_1  G00005 drug_005    moa_A CL00016  cellline_GB tissue_y
#>    3: plate_1  G00006 drug_006    moa_A CL00016  cellline_GB tissue_y
#>    4: plate_1  G00004 drug_004    moa_A CL00016  cellline_GB tissue_y
#>    5: plate_1  G00005 drug_005    moa_A CL00016  cellline_GB tissue_y
#>   ---                                                                
#> 3596: plate_3  G00026 drug_026    moa_E CL00017  cellline_HB tissue_y
#> 3597: plate_3  G00026 drug_026    moa_E CL00017  cellline_HB tissue_y
#> 3598: plate_3  G00026 drug_026    moa_E CL00017  cellline_HB tissue_y
#> 3599: plate_3  G00026 drug_026    moa_E CL00017  cellline_HB tissue_y
#> 3600: plate_3  G00026 drug_026    moa_E CL00017  cellline_HB tissue_y
#>       ReferenceDivisionTime Concentration Gnumber_2 DrugName_2 drug_moa_2
#>                       <num>         <num>    <char>     <char>     <char>
#>    1:                    46   0.001000000    G00021   drug_021      moa_D
#>    2:                    46   0.001000000    G00021   drug_021      moa_D
#>    3:                    46   0.001000000    G00021   drug_021      moa_D
#>    4:                    46   0.003162278    G00021   drug_021      moa_D
#>    5:                    46   0.003162278    G00021   drug_021      moa_D
#>   ---                                                                    
#> 3596:                    50   3.162277660   vehicle    vehicle    vehicle
#> 3597:                    50   3.162277660   vehicle    vehicle    vehicle
#> 3598:                    50   3.162277660   vehicle    vehicle    vehicle
#> 3599:                    50   3.162277660   vehicle    vehicle    vehicle
#> 3600:                    50   3.162277660   vehicle    vehicle    vehicle
#>       Concentration_2 ReadoutValue BackgroundValue Duration record_id
#>                 <num>        <num>           <num>    <num>     <int>
#>    1:           0.001     99.89992               0       72       727
#>    2:           0.001     96.69992               0       72       728
#>    3:           0.001     99.79992               0       72       729
#>    4:           0.001     98.69900               0       72       730
#>    5:           0.001     63.29936               0       72       731
#>   ---                                                                
#> 3596:           0.000    100.00000               0       72      3572
#> 3597:           0.000    100.00000               0       72      3573
#> 3598:           0.000    100.00000               0       72      3574
#> 3599:           0.000    100.00000               0       72      3575
#> 3600:           0.000    100.00000               0       72      3576
#> 
#> $`single-agent`
#>       Barcode Gnumber DrugName drug_moa    clid CellLineName   Tissue
#>        <char>  <char>   <char>   <char>  <char>       <char>   <char>
#>    1: plate_1  G00004 drug_004    moa_A CL00016  cellline_GB tissue_y
#>    2: plate_1  G00005 drug_005    moa_A CL00016  cellline_GB tissue_y
#>    3: plate_1  G00006 drug_006    moa_A CL00016  cellline_GB tissue_y
#>    4: plate_1  G00004 drug_004    moa_A CL00016  cellline_GB tissue_y
#>    5: plate_1  G00005 drug_005    moa_A CL00016  cellline_GB tissue_y
#>   ---                                                                
#> 1292: plate_3 vehicle  vehicle  vehicle CL00017  cellline_HB tissue_y
#> 1293: plate_3 vehicle  vehicle  vehicle CL00017  cellline_HB tissue_y
#> 1294: plate_3 vehicle  vehicle  vehicle CL00017  cellline_HB tissue_y
#> 1295: plate_3 vehicle  vehicle  vehicle CL00017  cellline_HB tissue_y
#> 1296: plate_3 vehicle  vehicle  vehicle CL00017  cellline_HB tissue_y
#>       ReferenceDivisionTime Concentration ReadoutValue BackgroundValue Duration
#>                       <num>         <num>        <num>           <num>    <num>
#>    1:                    46   0.001000000     99.89992               0       72
#>    2:                    46   0.001000000     96.69992               0       72
#>    3:                    46   0.001000000     99.79992               0       72
#>    4:                    46   0.003162278     98.69900               0       72
#>    5:                    46   0.003162278     63.29936               0       72
#>   ---                                                                          
#> 1292:                    50   0.000000000    100.00000               0       72
#> 1293:                    50   0.000000000    100.00000               0       72
#> 1294:                    50   0.000000000    100.00000               0       72
#> 1295:                    50   0.000000000    100.00000               0       72
#> 1296:                    50   0.000000000    100.00000               0       72
#>       record_id
#>           <int>
#>    1:         7
#>    2:         8
#>    3:         9
#>    4:        10
#>    5:        11
#>   ---          
#> 1292:       692
#> 1293:       693
#> 1294:       694
#> 1295:       695
#> 1296:       696
#> 

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
split_dt <- identify_data_type(input_dt)
split_raw_data(split_dt)
#> $`single-agent`
#>     ReadoutValue Concentration masked DrugName CellLineName Duration
#>            <num>         <num> <lgcl>   <char>       <char>    <num>
#>  1:            2           0.1  FALSE   DRUG_8        CELL1       72
#>  2:            3           0.2  FALSE  DRUG_10        CELL1       72
#>  3:            4           0.3  FALSE   DRUG_8        CELL1       72
#>  4:            2           0.1  FALSE   DRUG_8        CELL1       72
#>  5:            3           0.2  FALSE  DRUG_10        CELL1       72
#>  6:            4           0.3  FALSE   DRUG_8        CELL1       72
#>  7:            2           0.0  FALSE  DRUG_10        CELL1       72
#>  8:            2           0.0  FALSE  vehicle        CELL1       72
#>  9:            1           0.0  FALSE   DRUG_8        CELL1       72
#> 10:            1           0.0  FALSE  DRUG_10        CELL1       72
#> 11:            2           0.0  FALSE  vehicle        CELL1       72
#> 12:            1           0.0  FALSE   DRUG_8        CELL1       72
#> 13:            1           0.0  FALSE  DRUG_10        CELL1       72
#> 14:            1           0.0  FALSE  DRUG_10        CELL1       72
#>     CorrectedReadout2 record_id
#>                 <num>     <int>
#>  1:                 2         8
#>  2:                 3         9
#>  3:                 4        10
#>  4:                 2        12
#>  5:                 3        13
#>  6:                 4        14
#>  7:                 2         1
#>  8:                 2         2
#>  9:                 1         3
#> 10:                 1         4
#> 11:                 2         5
#> 12:                 1         6
#> 13:                 1         7
#> 14:                 1        11
#> 
```
