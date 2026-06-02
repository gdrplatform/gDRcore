# Transform se into raw_data

Transform se into raw_data

## Usage

``` r
convert_se_to_raw_data(se)
```

## Arguments

- se:

  SummarizedExperiment object with "RawTreated" and "Controls" assays

## Value

data.table with raw data

## Examples

``` r
mae <- gDRutils::get_synthetic_data("finalMAE_small")
se <- mae[[1]]
convert_se_to_raw_data(se)
#>       Barcode Concentration BackgroundValue record_id ReadoutValue Gnumber
#>        <char>         <num>           <num>     <int>        <num>  <char>
#>    1: plate_1   0.001000000               0       601         93.5  G00002
#>    2: plate_1   0.003162278               0       901         74.8  G00002
#>    3: plate_1   0.010000000               0      1201         40.1  G00002
#>    4: plate_1   0.031622777               0      1501         33.2  G00002
#>    5: plate_1   0.100000000               0      1801         31.5  G00002
#>   ---                                                                     
#> 8696: plate_2   0.000000000               0       110        104.0 vehicle
#> 8697: plate_2   0.000000000               0       190        104.1 vehicle
#> 8698: plate_1   0.000000000               0        80        104.4 vehicle
#> 8699: plate_3   0.000000000               0       560        104.6 vehicle
#> 8700: plate_3   0.000000000               0       570        104.7 vehicle
#>       DrugName drug_moa Duration    clid CellLineName   Tissue
#>         <char>   <char>    <num>  <char>       <char>   <char>
#>    1: drug_002    moa_A       72 CL00011  cellline_BA tissue_x
#>    2: drug_002    moa_A       72 CL00011  cellline_BA tissue_x
#>    3: drug_002    moa_A       72 CL00011  cellline_BA tissue_x
#>    4: drug_002    moa_A       72 CL00011  cellline_BA tissue_x
#>    5: drug_002    moa_A       72 CL00011  cellline_BA tissue_x
#>   ---                                                         
#> 8696:  vehicle  vehicle       72 CL00020  cellline_KB tissue_z
#> 8697:  vehicle  vehicle       72 CL00020  cellline_KB tissue_z
#> 8698:  vehicle  vehicle       72 CL00020  cellline_KB tissue_z
#> 8699:  vehicle  vehicle       72 CL00020  cellline_KB tissue_z
#> 8700:  vehicle  vehicle       72 CL00020  cellline_KB tissue_z
#>       ReferenceDivisionTime
#>                       <num>
#>    1:                    26
#>    2:                    26
#>    3:                    26
#>    4:                    26
#>    5:                    26
#>   ---                      
#> 8696:                    62
#> 8697:                    62
#> 8698:                    62
#> 8699:                    62
#> 8700:                    62
```
