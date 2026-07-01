# merge_data

Merge all the input data into a single data.table

## Usage

``` r
merge_data(
  manifest,
  treatments,
  data,
  cell_line_annotation = NULL,
  drug_annotation = NULL
)
```

## Arguments

- manifest:

  a data.table with a manifest info

- treatments:

  a data.table with a treaatments info

- data:

  a data.table with a raw data info

- cell_line_annotation:

  optional data.table with cell line annotations; if NULL (default),
  annotations are looked up from gDRinternal or gDRtestData

- drug_annotation:

  optional data.table with drug annotations; if NULL (default),
  annotations are looked up from gDRinternal or gDRtestData

## Value

a data.table with merged data and metadata.

## Examples

``` r
td <- gDRimport::get_test_data()
l_tbl <- gDRimport::load_data(
  manifest_file = gDRimport::manifest_path(td),
  df_template_files = gDRimport::template_path(td),
  results_file = gDRimport::result_path(td)
)
#> INFO [2026-07-01 22:57:28] Manifest loaded successfully
#> INFO [2026-07-01 22:57:28] Reading Template_7daytreated.xlsx with load_templates_xlsx
#> INFO [2026-07-01 22:57:28] Reading Template_Untreated.xlsx with load_templates_xlsx
#> INFO [2026-07-01 22:57:28] Loading Template_7daytreated.xlsx
#> INFO [2026-07-01 22:57:28] Loading Template_Untreated.xlsx
#> INFO [2026-07-01 22:57:28] Templates loaded successfully!
#> INFO [2026-07-01 22:57:28] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day0.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-07-01 22:57:28] Plate 201904190a read; 384 wells
#> INFO [2026-07-01 22:57:28] Plate 201904190b read; 384 wells
#> INFO [2026-07-01 22:57:28] Plate 201904190c read; 384 wells
#> INFO [2026-07-01 22:57:28] Plate 201904190d read; 384 wells
#> INFO [2026-07-01 22:57:28] Plate 201904190e read; 384 wells
#> INFO [2026-07-01 22:57:28] Plate 201904190f read; 384 wells
#> INFO [2026-07-01 22:57:28] File done
#> INFO [2026-07-01 22:57:28] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day7.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-07-01 22:57:29] Plate 201904197a read; 384 wells
#> INFO [2026-07-01 22:57:29] Plate 201904197b read; 384 wells
#> INFO [2026-07-01 22:57:29] Plate 201904197c read; 384 wells
#> INFO [2026-07-01 22:57:29] Plate 201904197d read; 384 wells
#> INFO [2026-07-01 22:57:29] Plate 201904197e read; 384 wells
#> INFO [2026-07-01 22:57:29] Plate 201904197f read; 384 wells
#> INFO [2026-07-01 22:57:29] File done
merge_data(
  l_tbl$manifest,
  l_tbl$treatments,
  l_tbl$data
)
#> INFO [2026-07-01 22:57:29] Merging data
#> INFO [2026-07-01 22:57:29] Merging the metadata (manifest and treatment files)
#> WARN [2026-07-01 22:57:29] 4608 well loaded, 768 wells discarded for lack of annotation,
#>     3840 data point selected
#>       CellLineName Tissue Duration DrugName Concentration DrugName_2
#>             <char> <char>    <num>   <char>         <num>     <char>
#>    1:  cellline_BA breast        0  vehicle             0    vehicle
#>    2:  cellline_BA breast        0  vehicle             0    vehicle
#>    3:  cellline_BA breast        0  vehicle             0    vehicle
#>    4:  cellline_BA breast        0  vehicle             0    vehicle
#>    5:  cellline_BA breast        0  vehicle             0    vehicle
#>   ---                                                               
#> 3836:  cellline_IB breast      168  vehicle             0    vehicle
#> 3837:  cellline_IB breast      168  vehicle             0    vehicle
#> 3838:  cellline_IB breast      168  vehicle             0    vehicle
#> 3839:  cellline_IB breast      168  vehicle             0    vehicle
#> 3840:  cellline_IB breast      168  vehicle             0    vehicle
#>       Concentration_2 drug_moa drug_moa_2 parental_identifier subtype
#>                 <num>   <char>     <char>              <char>  <char>
#>    1:               0  vehicle    vehicle         cellline_BA unknown
#>    2:               0  vehicle    vehicle         cellline_BA unknown
#>    3:               0  vehicle    vehicle         cellline_BA unknown
#>    4:               0  vehicle    vehicle         cellline_BA unknown
#>    5:               0  vehicle    vehicle         cellline_BA unknown
#>   ---                                                                
#> 3836:               0  vehicle    vehicle         cellline_IB unknown
#> 3837:               0  vehicle    vehicle         cellline_IB unknown
#> 3838:               0  vehicle    vehicle         cellline_IB unknown
#> 3839:               0  vehicle    vehicle         cellline_IB unknown
#> 3840:               0  vehicle    vehicle         cellline_IB unknown
#>          Barcode                  Template ReadoutValue BackgroundValue
#>           <char>                    <char>        <num>           <num>
#>    1: 201904190a   Template_Untreated.xlsx        91452             570
#>    2: 201904190a   Template_Untreated.xlsx       126448             570
#>    3: 201904190a   Template_Untreated.xlsx        91461             570
#>    4: 201904190a   Template_Untreated.xlsx       126449             570
#>    5: 201904190a   Template_Untreated.xlsx        91459             570
#>   ---                                                                  
#> 3836: 201904197f Template_7daytreated.xlsx       788743             395
#> 3837: 201904197f Template_7daytreated.xlsx       359748             395
#> 3838: 201904197f Template_7daytreated.xlsx       405491             395
#> 3839: 201904197f Template_7daytreated.xlsx       575063             395
#> 3840: 201904197f Template_7daytreated.xlsx       854686             395
#>       ReferenceDivisionTime    clid Gnumber Gnumber_2 WellRow WellColumn
#>                       <num>  <char>  <char>    <char>  <char>     <char>
#>    1:                    26 CL00011 vehicle   vehicle       A          3
#>    2:                    26 CL00011 vehicle   vehicle       B          3
#>    3:                    26 CL00011 vehicle   vehicle       C          3
#>    4:                    26 CL00011 vehicle   vehicle       D          3
#>    5:                    26 CL00011 vehicle   vehicle       E          3
#>   ---                                                                   
#> 3836:                    54 CL00018 vehicle   vehicle       D         22
#> 3837:                    54 CL00018 vehicle   vehicle       I         22
#> 3838:                    54 CL00018 vehicle   vehicle       J         22
#> 3839:                    54 CL00018 vehicle   vehicle       K         22
#> 3840:                    54 CL00018 vehicle   vehicle       L         22
```
