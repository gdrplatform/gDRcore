# Value Matching

Returns a lookup table or list of the positions of ALL matches of its
first argument in its second and vice versa. Similar to
[`match`](https://rdrr.io/r/base/match.html), though that function only
returns the first match.

## Usage

``` r
grr_matches(
  x,
  y,
  all.x = TRUE,
  all.y = TRUE,
  list = FALSE,
  indexes = TRUE,
  nomatch = NA
)
```

## Arguments

- x:

  vector. The values to be matched. Long vectors are not currently
  supported.

- y:

  vector. The values to be matched. Long vectors are not currently
  supported.

- all.x:

  logical; if `TRUE`, then each value in `x` will be included even if it
  has no matching values in `y`

- all.y:

  logical; if `TRUE`, then each value in `y` will be included even if it
  has no matching values in `x`

- list:

  logical. If `TRUE`, the result will be returned as a list of vectors,
  each vector being the matching values in y. If `FALSE`, result is
  returned as a data.table with repeated values for each match.

- indexes:

  logical. Whether to return the indices of the matches or the actual
  values.

- nomatch:

  the value to be returned in the case when no match is found. If not
  provided and `indexes=TRUE`, items with no match will be represented
  as `NA`. If set to `NULL`, items with no match will be set to an index
  value of `length+1`. If `indexes=FALSE`, they will default to `NA`.

## Value

data.table

## Details

This behavior can be imitated by using joins to create lookup tables,
but `matches` is simpler and faster: usually faster than the best joins
in other packages and thousands of times faster than the built in
[`merge`](https://rdrr.io/r/base/merge.html).

`all.x/all.y` correspond to the four types of database joins in the
following way:

- left:

  `all.x=TRUE`, `all.y=FALSE`

- right:

  `all.x=FALSE`, `all.y=TRUE`

- inner:

  `all.x=FALSE`, `all.y=FALSE`

- full:

  `all.x=TRUE`, `all.y=TRUE`

Note that `NA` values will match other `NA` values.

Source of the function: https://github.com/cran/grr/blob/master/R/grr.R

## Examples

``` r
mat_elem <- data.table::data.table(
  DrugName = rep(c("untreated", "drugA", "drugB", "untreated"), 2),
  DrugName_2 = rep(c("untreated", "vehicle", "drugA", "drugB"), 2),
  clid = rep(c("C1", "C2"), each = 4)
)
untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
ref_idx <- which(
  mat_elem$DrugName %in% untreated_tag |
   mat_elem$DrugName_2 %in% untreated_tag
)
ref <- mat_elem[ref_idx, ]
treated <- mat_elem[-ref_idx, ]
valid <- c("DrugName", "DrugName_2")
trt <- lapply(valid, function(x) {
  colnames <- c("clid", x)
  treated[, colnames, with = FALSE]
})
trt <- do.call(paste,
  do.call(rbind, lapply(trt, function(x) setNames(x, names(trt[[1]]))))
)
ref <- lapply(valid, function(x) {
  colnames <- c("clid", x)
  ref[, colnames, with = FALSE]
})
ref <- do.call(paste,
  do.call(rbind, lapply(ref, function(x) setNames(x, names(ref[[1]]))))
)
grr_matches(trt, ref, list = FALSE, all.y = FALSE)
#>        x     y
#>    <int> <int>
#> 1:     3     2
#> 2:     1     9
#> 3:     4     5
#> 4:     2    12
```
