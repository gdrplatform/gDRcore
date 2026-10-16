# Barcodes of plates imaged more than once

Barcodes of plates imaged more than once

## Usage

``` r
.time_course_plates(dt)
```

## Arguments

- dt:

  data.table of raw drug response data

## Value

list with `time_course`, the plates imaged more than once that have a
reading at time zero and at least three time points, and `no_time_zero`,
the plates imaged more than once without a reading at time zero
