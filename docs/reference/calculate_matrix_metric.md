# Calculate a metric for combination data.

Calculate a metric based off of single-agent values in combination
screens.

## Usage

``` r
calculate_HSA(sa1, series_id1, sa2, series_id2, metric)

calculate_Bliss(
  sa1,
  series_id1,
  sa2,
  series_id2,
  metric,
  measured_col = "smooth"
)

.calculate_matrix_metric(
  sa1,
  series_id1,
  sa2,
  series_id2,
  metric,
  FXN,
  measured_col = "x"
)
```

## Arguments

- sa1:

  data.table containing single agent data where entries in `series_id2`
  are all `0`. Columns of the data.table include identifiers and the
  `metric` of interest. Metric is stored in the 'x' column.

- series_id1:

  String representing the column within `sa1` that represents id1.

- sa2:

  data.table containing single agent data where entries in `series_id1`
  are all `0`. Columns of the data.table include identifiers and the
  `metric` of interest.n Metric is stored in the 'x' column.

- series_id2:

  String representing the column within `sa2` that represents id2.

- metric:

  String specifying the metric of interest. Usually either 'GRvalue' or
  'RelativeViability'.

- measured_col:

  String specyfying the measured colname.

- FXN:

  Function to apply to the single-agent fits to calculate a metric.

## Value

data.table containing a single row for every unique combination of the
two series identifiers and the corresponding calculated metric for each
row.

## Details

`calculate_HSA` takes the minimum of the two single agents readouts.
`calculate_Bliss` performs Bliss additivity calculation based on the
single agent effects, defined as `1-x` for the corresponding
normalization. See
https://www.sciencedirect.com/science/article/pii/S1359644619303460?via%3Dihub#tb0005
for more details.

## Examples

``` r
n <- 10
sa1 <- data.table::data.table(conc = seq(n), conc2 = rep(0, n), smooth = seq(n))
sa2 <- data.table::data.table(conc = rep(0, n), conc2 = seq(n), smooth = seq(n))
calculate_HSA(sa1, "conc", sa2, "conc2", "smooth")
#>       conc conc2 metric1 metric2 metric
#>      <int> <int>   <int>   <int>  <int>
#>   1:     1     1       1       1      1
#>   2:     2     1       2       1      1
#>   3:     3     1       3       1      1
#>   4:     4     1       4       1      1
#>   5:     5     1       5       1      1
#>   6:     6     1       6       1      1
#>   7:     7     1       7       1      1
#>   8:     8     1       8       1      1
#>   9:     9     1       9       1      1
#>  10:    10     1      10       1      1
#>  11:     1     2       1       2      1
#>  12:     2     2       2       2      2
#>  13:     3     2       3       2      2
#>  14:     4     2       4       2      2
#>  15:     5     2       5       2      2
#>  16:     6     2       6       2      2
#>  17:     7     2       7       2      2
#>  18:     8     2       8       2      2
#>  19:     9     2       9       2      2
#>  20:    10     2      10       2      2
#>  21:     1     3       1       3      1
#>  22:     2     3       2       3      2
#>  23:     3     3       3       3      3
#>  24:     4     3       4       3      3
#>  25:     5     3       5       3      3
#>  26:     6     3       6       3      3
#>  27:     7     3       7       3      3
#>  28:     8     3       8       3      3
#>  29:     9     3       9       3      3
#>  30:    10     3      10       3      3
#>  31:     1     4       1       4      1
#>  32:     2     4       2       4      2
#>  33:     3     4       3       4      3
#>  34:     4     4       4       4      4
#>  35:     5     4       5       4      4
#>  36:     6     4       6       4      4
#>  37:     7     4       7       4      4
#>  38:     8     4       8       4      4
#>  39:     9     4       9       4      4
#>  40:    10     4      10       4      4
#>  41:     1     5       1       5      1
#>  42:     2     5       2       5      2
#>  43:     3     5       3       5      3
#>  44:     4     5       4       5      4
#>  45:     5     5       5       5      5
#>  46:     6     5       6       5      5
#>  47:     7     5       7       5      5
#>  48:     8     5       8       5      5
#>  49:     9     5       9       5      5
#>  50:    10     5      10       5      5
#>  51:     1     6       1       6      1
#>  52:     2     6       2       6      2
#>  53:     3     6       3       6      3
#>  54:     4     6       4       6      4
#>  55:     5     6       5       6      5
#>  56:     6     6       6       6      6
#>  57:     7     6       7       6      6
#>  58:     8     6       8       6      6
#>  59:     9     6       9       6      6
#>  60:    10     6      10       6      6
#>  61:     1     7       1       7      1
#>  62:     2     7       2       7      2
#>  63:     3     7       3       7      3
#>  64:     4     7       4       7      4
#>  65:     5     7       5       7      5
#>  66:     6     7       6       7      6
#>  67:     7     7       7       7      7
#>  68:     8     7       8       7      7
#>  69:     9     7       9       7      7
#>  70:    10     7      10       7      7
#>  71:     1     8       1       8      1
#>  72:     2     8       2       8      2
#>  73:     3     8       3       8      3
#>  74:     4     8       4       8      4
#>  75:     5     8       5       8      5
#>  76:     6     8       6       8      6
#>  77:     7     8       7       8      7
#>  78:     8     8       8       8      8
#>  79:     9     8       9       8      8
#>  80:    10     8      10       8      8
#>  81:     1     9       1       9      1
#>  82:     2     9       2       9      2
#>  83:     3     9       3       9      3
#>  84:     4     9       4       9      4
#>  85:     5     9       5       9      5
#>  86:     6     9       6       9      6
#>  87:     7     9       7       9      7
#>  88:     8     9       8       9      8
#>  89:     9     9       9       9      9
#>  90:    10     9      10       9      9
#>  91:     1    10       1      10      1
#>  92:     2    10       2      10      2
#>  93:     3    10       3      10      3
#>  94:     4    10       4      10      4
#>  95:     5    10       5      10      5
#>  96:     6    10       6      10      6
#>  97:     7    10       7      10      7
#>  98:     8    10       8      10      8
#>  99:     9    10       9      10      9
#> 100:    10    10      10      10     10
#>       conc conc2 metric1 metric2 metric
#>      <int> <int>   <int>   <int>  <int>
n <- 10
sa1 <- data.table::data.table(conc = seq(n), conc2 = rep(0, n), smooth = seq(n))
sa2 <- data.table::data.table(conc = rep(0, n), conc2 = seq(n), smooth = seq(n))
calculate_Bliss(sa1, "conc", sa2, "conc2", "smooth")
#>       conc conc2 metric1 metric2 metric
#>      <int> <int>   <int>   <int>  <int>
#>   1:     1     1       1       1      1
#>   2:     2     1       2       1      2
#>   3:     3     1       3       1      3
#>   4:     4     1       4       1      4
#>   5:     5     1       5       1      5
#>   6:     6     1       6       1      6
#>   7:     7     1       7       1      7
#>   8:     8     1       8       1      8
#>   9:     9     1       9       1      9
#>  10:    10     1      10       1     10
#>  11:     1     2       1       2      2
#>  12:     2     2       2       2      4
#>  13:     3     2       3       2      6
#>  14:     4     2       4       2      8
#>  15:     5     2       5       2     10
#>  16:     6     2       6       2     12
#>  17:     7     2       7       2     14
#>  18:     8     2       8       2     16
#>  19:     9     2       9       2     18
#>  20:    10     2      10       2     20
#>  21:     1     3       1       3      3
#>  22:     2     3       2       3      6
#>  23:     3     3       3       3      9
#>  24:     4     3       4       3     12
#>  25:     5     3       5       3     15
#>  26:     6     3       6       3     18
#>  27:     7     3       7       3     21
#>  28:     8     3       8       3     24
#>  29:     9     3       9       3     27
#>  30:    10     3      10       3     30
#>  31:     1     4       1       4      4
#>  32:     2     4       2       4      8
#>  33:     3     4       3       4     12
#>  34:     4     4       4       4     16
#>  35:     5     4       5       4     20
#>  36:     6     4       6       4     24
#>  37:     7     4       7       4     28
#>  38:     8     4       8       4     32
#>  39:     9     4       9       4     36
#>  40:    10     4      10       4     40
#>  41:     1     5       1       5      5
#>  42:     2     5       2       5     10
#>  43:     3     5       3       5     15
#>  44:     4     5       4       5     20
#>  45:     5     5       5       5     25
#>  46:     6     5       6       5     30
#>  47:     7     5       7       5     35
#>  48:     8     5       8       5     40
#>  49:     9     5       9       5     45
#>  50:    10     5      10       5     50
#>  51:     1     6       1       6      6
#>  52:     2     6       2       6     12
#>  53:     3     6       3       6     18
#>  54:     4     6       4       6     24
#>  55:     5     6       5       6     30
#>  56:     6     6       6       6     36
#>  57:     7     6       7       6     42
#>  58:     8     6       8       6     48
#>  59:     9     6       9       6     54
#>  60:    10     6      10       6     60
#>  61:     1     7       1       7      7
#>  62:     2     7       2       7     14
#>  63:     3     7       3       7     21
#>  64:     4     7       4       7     28
#>  65:     5     7       5       7     35
#>  66:     6     7       6       7     42
#>  67:     7     7       7       7     49
#>  68:     8     7       8       7     56
#>  69:     9     7       9       7     63
#>  70:    10     7      10       7     70
#>  71:     1     8       1       8      8
#>  72:     2     8       2       8     16
#>  73:     3     8       3       8     24
#>  74:     4     8       4       8     32
#>  75:     5     8       5       8     40
#>  76:     6     8       6       8     48
#>  77:     7     8       7       8     56
#>  78:     8     8       8       8     64
#>  79:     9     8       9       8     72
#>  80:    10     8      10       8     80
#>  81:     1     9       1       9      9
#>  82:     2     9       2       9     18
#>  83:     3     9       3       9     27
#>  84:     4     9       4       9     36
#>  85:     5     9       5       9     45
#>  86:     6     9       6       9     54
#>  87:     7     9       7       9     63
#>  88:     8     9       8       9     72
#>  89:     9     9       9       9     81
#>  90:    10     9      10       9     90
#>  91:     1    10       1      10     10
#>  92:     2    10       2      10     20
#>  93:     3    10       3      10     30
#>  94:     4    10       4      10     40
#>  95:     5    10       5      10     50
#>  96:     6    10       6      10     60
#>  97:     7    10       7      10     70
#>  98:     8    10       8      10     80
#>  99:     9    10       9      10     90
#> 100:    10    10      10      10    100
#>       conc conc2 metric1 metric2 metric
#>      <int> <int>   <int>   <int>  <int>
```
