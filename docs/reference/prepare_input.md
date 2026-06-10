# Prepare input data common for all experiments

Current steps

- refining nested confounders

- refining nested identifiers

- splitting df\_ into (per experiment) df_list

## Usage

``` r
prepare_input(x, ...)
```

## Arguments

- x:

  data.table with raw data or MAE object with dose-response data

- ...:

  additional parameters

## Value

list of input data

## Examples

``` r
td <- gDRimport::get_test_data()
l_tbl <- gDRimport::load_data(
  manifest_file = gDRimport::manifest_path(td),
  df_template_files = gDRimport::template_path(td),
  results_file = gDRimport::result_path(td)
)
#> INFO [2026-06-10 10:01:28] Manifest loaded successfully
#> INFO [2026-06-10 10:01:28] Reading Template_7daytreated.xlsx with load_templates_xlsx
#> INFO [2026-06-10 10:01:28] Reading Template_Untreated.xlsx with load_templates_xlsx
#> INFO [2026-06-10 10:01:28] Loading Template_7daytreated.xlsx
#> INFO [2026-06-10 10:01:28] Loading Template_Untreated.xlsx
#> INFO [2026-06-10 10:01:28] Templates loaded successfully!
#> INFO [2026-06-10 10:01:28] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day0.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-10 10:01:28] Plate 201904190a read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904190b read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904190c read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904190d read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904190e read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904190f read; 384 wells
#> INFO [2026-06-10 10:01:28] File done
#> INFO [2026-06-10 10:01:28] Reading file /home/runner/work/_temp/Library/gDRimport/extdata/data1/RawData_day7.xlsx, sheet Readout_0077vs0068_day7
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
#> INFO [2026-06-10 10:01:28] Plate 201904197a read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904197b read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904197c read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904197d read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904197e read; 384 wells
#> INFO [2026-06-10 10:01:28] Plate 201904197f read; 384 wells
#> INFO [2026-06-10 10:01:28] File done
df_ <- merge_data(
  l_tbl$manifest,
  l_tbl$treatments,
  l_tbl$data
)
#> INFO [2026-06-10 10:01:28] Merging data
#> INFO [2026-06-10 10:01:28] Merging the metadata (manifest and treatment files)
#> WARN [2026-06-10 10:01:28] 4608 well loaded, 768 wells discarded for lack of annotation,
#>     3840 data point selected
nested_confounders = intersect(
  names(df_),
  gDRutils::get_env_identifiers("barcode")
)
prepare_input(df_, nested_confounders, NULL)
#> $df_
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
#>       record_id    type
#>           <int>  <char>
#>    1:         1 control
#>    2:         2 control
#>    3:         3 control
#>    4:         4 control
#>    5:         5 control
#>   ---                  
#> 3836:      3836 control
#> 3837:      3837 control
#> 3838:      3838 control
#> 3839:      3839 control
#> 3840:      3840 control
#> 
#> $df_list
#> $df_list$combination
#>       CellLineName Tissue Duration DrugName Concentration DrugName_2
#>             <char> <char>    <num>   <char>         <num>     <char>
#>    1:  cellline_BA breast      168 drug_002   0.001524158   drug_011
#>    2:  cellline_BA breast      168 drug_002   0.001524158   drug_011
#>    3:  cellline_BA breast      168 drug_002   0.001524158   drug_011
#>    4:  cellline_BA breast      168 drug_002   0.001524158   drug_011
#>    5:  cellline_BA breast      168 drug_002   0.001524158   drug_011
#>   ---                                                               
#> 3836:  cellline_IB breast      168 drug_011   0.149999911    vehicle
#> 3837:  cellline_IB breast      168 drug_011   0.149999911    vehicle
#> 3838:  cellline_IB breast      168 drug_011   0.149999911    vehicle
#> 3839:  cellline_IB breast      168 drug_011   0.149999911    vehicle
#> 3840:  cellline_IB breast      168 drug_011   0.149999911    vehicle
#>       Concentration_2 drug_moa drug_moa_2 parental_identifier subtype
#>                 <num>   <char>     <char>              <char>  <char>
#>    1:       0.1499999    moa_A      moa_B         cellline_BA unknown
#>    2:       0.1499999    moa_A      moa_B         cellline_BA unknown
#>    3:       0.1499999    moa_A      moa_B         cellline_BA unknown
#>    4:       0.1499999    moa_A      moa_B         cellline_BA unknown
#>    5:       0.1499999    moa_A      moa_B         cellline_BA unknown
#>   ---                                                                
#> 3836:       0.0000000    moa_B    vehicle         cellline_IB unknown
#> 3837:       0.0000000    moa_B    vehicle         cellline_IB unknown
#> 3838:       0.0000000    moa_B    vehicle         cellline_IB unknown
#> 3839:       0.0000000    moa_B    vehicle         cellline_IB unknown
#> 3840:       0.0000000    moa_B    vehicle         cellline_IB unknown
#>          Barcode                  Template ReadoutValue BackgroundValue
#>           <char>                    <char>        <num>           <num>
#>    1: 201904197a Template_7daytreated.xlsx       102301             570
#>    2: 201904197a Template_7daytreated.xlsx        76966             570
#>    3: 201904197a Template_7daytreated.xlsx       461220             570
#>    4: 201904197a Template_7daytreated.xlsx       497047             570
#>    5: 201904197a Template_7daytreated.xlsx        64611             570
#>   ---                                                                  
#> 3836: 201904197f Template_7daytreated.xlsx       383366             395
#> 3837: 201904197f Template_7daytreated.xlsx       133207             395
#> 3838: 201904197f Template_7daytreated.xlsx       204959             395
#> 3839: 201904197f Template_7daytreated.xlsx       323669             395
#> 3840: 201904197f Template_7daytreated.xlsx       387380             395
#>       ReferenceDivisionTime    clid Gnumber Gnumber_2 WellRow WellColumn
#>                       <num>  <char>  <char>    <char>  <char>     <char>
#>    1:                    26 CL00011  G00002    G00011       E         19
#>    2:                    26 CL00011  G00002    G00011       F         19
#>    3:                    26 CL00011  G00002    G00011       G         19
#>    4:                    26 CL00011  G00002    G00011       H         19
#>    5:                    26 CL00011  G00002    G00011       E         20
#>   ---                                                                   
#> 3836:                    54 CL00018  G00011   vehicle       H         22
#> 3837:                    54 CL00018  G00011   vehicle       M         22
#> 3838:                    54 CL00018  G00011   vehicle       N         22
#> 3839:                    54 CL00018  G00011   vehicle       O         22
#> 3840:                    54 CL00018  G00011   vehicle       P         22
#>       record_id
#>           <int>
#>    1:       321
#>    2:       322
#>    3:       323
#>    4:       324
#>    5:       325
#>   ---          
#> 3836:      3820
#> 3837:      3821
#> 3838:      3822
#> 3839:      3823
#> 3840:      3824
#> 
#> $df_list$`single-agent`
#>       CellLineName Tissue Duration DrugName Concentration drug_moa
#>             <char> <char>    <num>   <char>         <num>   <char>
#>    1:  cellline_BA breast      168 drug_002   0.001524158    moa_A
#>    2:  cellline_BA breast      168 drug_002   0.001524158    moa_A
#>    3:  cellline_BA breast      168 drug_002   0.001524158    moa_A
#>    4:  cellline_BA breast      168 drug_002   0.001524158    moa_A
#>    5:  cellline_BA breast      168 drug_002   0.001524158    moa_A
#>   ---                                                             
#> 2972:  cellline_IB breast      168  vehicle   0.000000000  vehicle
#> 2973:  cellline_IB breast      168  vehicle   0.000000000  vehicle
#> 2974:  cellline_IB breast      168  vehicle   0.000000000  vehicle
#> 2975:  cellline_IB breast      168  vehicle   0.000000000  vehicle
#> 2976:  cellline_IB breast      168  vehicle   0.000000000  vehicle
#>       parental_identifier subtype    Barcode                  Template
#>                    <char>  <char>     <char>                    <char>
#>    1:         cellline_BA unknown 201904197a Template_7daytreated.xlsx
#>    2:         cellline_BA unknown 201904197a Template_7daytreated.xlsx
#>    3:         cellline_BA unknown 201904197a Template_7daytreated.xlsx
#>    4:         cellline_BA unknown 201904197a Template_7daytreated.xlsx
#>    5:         cellline_BA unknown 201904197a Template_7daytreated.xlsx
#>   ---                                                                 
#> 2972:         cellline_IB unknown 201904197f Template_7daytreated.xlsx
#> 2973:         cellline_IB unknown 201904197f Template_7daytreated.xlsx
#> 2974:         cellline_IB unknown 201904197f Template_7daytreated.xlsx
#> 2975:         cellline_IB unknown 201904197f Template_7daytreated.xlsx
#> 2976:         cellline_IB unknown 201904197f Template_7daytreated.xlsx
#>       ReadoutValue BackgroundValue ReferenceDivisionTime    clid Gnumber
#>              <num>           <num>                 <num>  <char>  <char>
#>    1:       159679             570                    26 CL00011  G00002
#>    2:       165488             570                    26 CL00011  G00002
#>    3:      1169641             570                    26 CL00011  G00002
#>    4:      1346753             570                    26 CL00011  G00002
#>    5:       112168             570                    26 CL00011  G00002
#>   ---                                                                   
#> 2972:       788743             395                    54 CL00018 vehicle
#> 2973:       359748             395                    54 CL00018 vehicle
#> 2974:       405491             395                    54 CL00018 vehicle
#> 2975:       575063             395                    54 CL00018 vehicle
#> 2976:       854686             395                    54 CL00018 vehicle
#>       WellRow WellColumn record_id
#>        <char>     <char>     <int>
#>    1:       A         19       329
#>    2:       B         19       330
#>    3:       C         19       331
#>    4:       D         19       332
#>    5:       A         20       333
#>   ---                             
#> 2972:       D         22      3836
#> 2973:       I         22      3837
#> 2974:       J         22      3838
#> 2975:       K         22      3839
#> 2976:       L         22      3840
#> 
#> 
#> $nested_confounders
#> [1] "Barcode"
#> 
#> $nested_identifiers_l
#> $nested_identifiers_l$`single-agent`
#> [1] "Concentration"
#> 
#> $nested_identifiers_l$combination
#> [1] "Concentration"   "Concentration_2"
#> 
#> $nested_identifiers_l$`time-course`
#> [1] "Concentration"   "Concentration_2"
#> 
#> 
#> $exps
#> $exps$combination
#> NULL
#> 
#> $exps$`single-agent`
#> NULL
#> 
#> 
```
