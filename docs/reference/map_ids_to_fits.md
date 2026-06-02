# Get predicted values for a given fit and input.

Map fittings to identifiers and compute the predicted values for
corresponding fits.

## Usage

``` r
map_ids_to_fits(pred, match_col, fittings, fitting_id_col)
```

## Arguments

- pred:

  numeric vector for which you want predictions.

- match_col:

  vector to match on `fittings` to get the correct fit.

- fittings:

  data.table of fit metrics.

- fitting_id_col:

  string of the column name in `fittings` that should be used to match
  with `match_col` .

## Value

Numeric vector of predicted values given `pred` inputs and `fittings`
values.

## Examples

``` r
pred <- c(1, 5, 5)
match_col <- c(1, 1, 2)
fitting_id_col <- "match_on_me"

fit1 <- data.table::data.table(h = 2.09, x_inf = 0.68, x_0 = 1, ec50 = 0.003)
fit2 <- data.table::data.table(h = 0.906, x_inf = 0.46, x_0 = 1, ec50 = 0.001)
fittings <- do.call(rbind, list(fit1, fit2))
fittings[[fitting_id_col]] <- c(1, 2)

map_ids_to_fits(pred, match_col, fittings, fitting_id_col)
#> [1] 0.6800017 0.6800001 0.4602404
```
