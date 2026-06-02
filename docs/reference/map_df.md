# Map treated conditions to their respective references.

Map treated conditions to their respective Day0, untreated, or
single-agent references using condition metadata.

## Usage

``` r
map_df(
  trt_md,
  ref_md,
  override_untrt_controls = NULL,
  ref_cols,
  ref_type = c("Day0", "untrt_Endpoint")
)
```

## Arguments

- trt_md:

  data.table of treated metadata.

- ref_md:

  data.table of untreated metadata.

- override_untrt_controls:

  named list indicating what treatment metadata fields should be used as
  a control. Defaults to `NULL`.

- ref_cols:

  character vector of the names of reference columns to include. Likely
  obtained from
  [`identify_keys()`](https://gdrplatform.github.io/gDRcore/reference/identify_keys.md).

- ref_type:

  string of the reference type to map to. Should be one of
  `c("Day0", "untrt_Endpoint")`.

## Value

named list mapping treated metadata to untreated metadata.

## Details

If `override_untrt_controls` is specified, the values in the named list
will supersede the values in `trt_md` during the matching process. This
is useful for mapping treatments to specific "standard" untreated
controls.

## See also

identify_keys

## Examples

``` r
# Standard Endpoint Mapping
trt_dt <- data.table::data.table(
  clid = c("C1", "C2"),
  Duration = 72,
  rn = c("T1", "T2")
)
ref_dt <- data.table::data.table(
  clid = c("C1", "C2"),
  Duration = 72,
  rn = c("R1", "R2")
)
map_df(trt_dt, ref_dt, ref_cols = "clid", ref_type = "untrt_Endpoint")
#> $T1
#> [1] "R1"
#> 
#> $T2
#> [1] "R2"
#> 
```
