# Transform mae into raw data

Transform mae into raw data

## Usage

``` r
convert_mae_to_raw_data(mae)
```

## Arguments

- mae:

  MultiAssayExperiment object with SummarizedExperiments containing
  "RawTreated" and "Controls" assays

## Value

data.table with raw data

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_small")
convert_mae_to_raw_data(mae)
#> Loading required namespace: BumpyMatrix
#>       Barcode Concentration ReadoutValue Gnumber DrugName drug_moa Duration
#>        <char>         <num>        <num>  <char>   <char>   <char>    <num>
#>    1: plate_1             0         95.7 vehicle  vehicle  vehicle       72
#>    2: plate_1             0        100.2 vehicle  vehicle  vehicle       72
#>    3: plate_1             0        102.6 vehicle  vehicle  vehicle       72
#>    4: plate_1             0        101.6 vehicle  vehicle  vehicle       72
#>    5: plate_1             0         99.9 vehicle  vehicle  vehicle       72
#>   ---                                                                      
#> 3296: plate_3            10         57.7  G00011 drug_011    moa_B       72
#> 3297: plate_3            10         37.7  G00011 drug_011    moa_B       72
#> 3298: plate_3            10         28.6  G00011 drug_011    moa_B       72
#> 3299: plate_3            10         29.6  G00011 drug_011    moa_B       72
#> 3300: plate_3            10         11.0  G00011 drug_011    moa_B       72
#>          clid CellLineName   Tissue ReferenceDivisionTime
#>        <char>       <char>   <char>                 <num>
#>    1: CL00011  cellline_BA tissue_x                    26
#>    2: CL00012  cellline_CA tissue_x                    30
#>    3: CL00013  cellline_DA tissue_x                    34
#>    4: CL00014  cellline_EA tissue_x                    38
#>    5: CL00015  cellline_FA tissue_x                    42
#>   ---                                                    
#> 3296: CL00016  cellline_GB tissue_y                    46
#> 3297: CL00017  cellline_HB tissue_y                    50
#> 3298: CL00018  cellline_IB tissue_y                    54
#> 3299: CL00019  cellline_JB tissue_z                    58
#> 3300: CL00020  cellline_KB tissue_z                    62
```
