# Flag rows of plates imaged more than once

A well read at several time points has several records, which is how
[`data_model()`](https://gdrplatform.github.io/gDRcore/reference/data_model.md)
recognizes time-course data. A plate also needs a reading at time zero,
the reference of the time-course normalization, and at least three time
points, so that an end-point plate read twice, or read at time zero and
at the end point, keeps its type. The check is made per plate.

## Usage

``` r
.is_time_course_plate(dt)
```

## Arguments

- dt:

  data.table of raw drug response data

## Value

logical vector, one element per row of `dt`
