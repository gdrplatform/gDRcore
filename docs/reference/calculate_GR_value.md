# Calculate a GR value.

Calculate a GR value for a given set of dose response values.

## Usage

``` r
calculate_GR_value(
  rel_viability,
  corrected_readout,
  day0_readout,
  untrt_readout,
  ndigit_rounding,
  duration,
  ref_div_time,
  cap = 1.25
)

calculate_time_dep_GR_value(
  corrected_readout,
  day0_readout,
  untrt_readout,
  ndigit_rounding
)

calculate_endpt_GR_value(
  rel_viability,
  duration,
  ref_div_time,
  cap = 1.25,
  ndigit_rounding
)
```

## Arguments

- rel_viability:

  numeric vector representing the Relative Viability.

- corrected_readout:

  numeric vector containing the corrected readout.

- day0_readout:

  numeric vector containing the day 0 readout.

- untrt_readout:

  numeric vector containing the untreated readout.

- ndigit_rounding:

  integer specifying the number of digits to use for calculation
  rounding.

- duration:

  numeric value specifying the length of time the cells were treated (in
  hours).

- ref_div_time:

  numeric value specifying the reference division time for the cell line
  in the experiment.

- cap:

  numeric value representing the value to cap the highest allowed
  relative viability at.

## Value

numeric vector containing GR values, one value for each element of the
input vectors.

## Details

Note that this function expects that all numeric vectors are of the same
length. `calculate_GR_value` will try to greedily calculate a GR value.
If no day 0 readouts are available, the `duration` and `ref_div_time`
will be used to try to back-calculate a day 0 value in order to produce
a GR value.

In the case of calculating the reference GR value from multiple
reference readout values, the vectorized calculation is performed and
then the resulting vector should be averaged outside of this function.

Note that it is expected that the `ref_div_time` and `duration` are
reported in the same units.

## See also

normalize_SE2

## Examples

``` r
duration <- 144
rv <- seq(0.1, 1, 0.1)
corrected <- seq(41000, 50000, 1000)
day0 <- seq(91000, 95500, 500)
untrt <- rep(c(115000, 118000), 5)

calculate_GR_value(
  rel_viability = rv,
  corrected_readout = corrected,
  day0_readout = day0,
  untrt_readout = untrt,
  ndigit_rounding = 4,
  duration = duration,
  ref_div_time = duration / 2
)
#>  [1] -0.9057 -0.8802 -0.9058 -0.8794 -0.9065 -0.8791 -0.9077 -0.8793 -0.9095
#> [10] -0.8800

readouts <- rep(10000, 5)
calculate_time_dep_GR_value(readouts, readouts * 1.32, readouts * 2, 2)
#> [1] -0.37 -0.37 -0.37 -0.37 -0.37

readouts <- rep(10000, 5)
calculate_endpt_GR_value(readouts, 72, 1, ndigit_rounding = 2)
#> [1] 1.01 1.01 1.01 1.01 1.01
```
