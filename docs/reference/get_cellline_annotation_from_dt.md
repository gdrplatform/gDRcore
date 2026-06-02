# Retrieve the cell line annotation from the annotated dt input

Retrieve the cell line annotation from the annotated dt input

## Usage

``` r
get_cellline_annotation_from_dt(dt)
```

## Arguments

- dt:

  annotated data.table

## Value

data.table with cell line annotation

## Examples

``` r
dt <- data.table::data.table(Gnumber = "A",
clid = "CL123",
CellLineName = "cl name",
Tissue = "Bone",
parental_identifier = "some cl",
subtype = "cortical",
ReferenceDivisionTime = 5)
get_cellline_annotation_from_dt(dt)
#>      clid CellLineName Tissue parental_identifier  subtype
#>    <char>       <char> <char>              <char>   <char>
#> 1:  CL123      cl name   Bone             some cl cortical
#>    ReferenceDivisionTime
#>                    <num>
#> 1:                     5
```
