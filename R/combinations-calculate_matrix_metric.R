#' Calculate a metric for combination data.
#'
#' Calculate a metric based off of single-agent values in combination screens.
#'
#' @param sa1 data.table containing single agent data where entries in 
#' \code{series_id2} are all \code{0}. Columns of the data.table include 
#' identifiers and the \code{metric} of interest. Metric is stored in the 'x' column.
#' @param series_id1 String representing the column within \code{sa1} that 
#' represents id1.
#' @param sa2 data.table containing single agent data where entries in 
#' \code{series_id1} are all \code{0}. Columns of the data.table include 
#' identifiers and the \code{metric} of interest.n Metric is stored in the 'x' column.
#' @param series_id2 String representing the column within \code{sa2} that 
#' represents id2.
#' @param metric String specifying the metric of interest. Usually either 'GRvalue'
#' or 'RelativeViability'.
#' @param measured_col String specyfying the measured colname.
#' @param FXN Function to apply to the single-agent fits to calculate a metric.
#'
#' @return data.table containing a single row for every unique combination of 
#' the two series identifiers and the corresponding calculated metric for 
#' each row.
#'
#' @name calculate_matrix_metric
#' @details
#' \code{calculate_HSA} takes the minimum of the two single agents readouts.
#' \code{calculate_Bliss} performs Bliss additivity calculation based on the 
#' single agent effects, defined as \code{1-x} for the corresponding 
#' normalization.
#' See https://www.sciencedirect.com/science/article/pii/S1359644619303460?via%3Dihub#tb0005
#' for more details.
NULL


#' @rdname calculate_matrix_metric
#' @examples
#' n <- 10
#' sa1 <- data.table::data.table(conc = seq(n), conc2 = rep(0, n), smooth = seq(n))
#' sa2 <- data.table::data.table(conc = rep(0, n), conc2 = seq(n), smooth = seq(n))
#' calculate_HSA(sa1, "conc", sa2, "conc2", "smooth")
#' @keywords combinations
#' @export
calculate_HSA <- function(sa1, series_id1, sa2, series_id2, metric) {
  .calculate_matrix_metric(sa1,
                           series_id1,
                           sa2,
                           series_id2,
                           metric,
                           FXN = pmin,
                           measured_col = "smooth")
}


#' @rdname calculate_matrix_metric
#' @examples
#' n <- 10
#' sa1 <- data.table::data.table(conc = seq(n), conc2 = rep(0, n), smooth = seq(n))
#' sa2 <- data.table::data.table(conc = rep(0, n), conc2 = seq(n), smooth = seq(n))
#' calculate_Bliss(sa1, "conc", sa2, "conc2", "smooth")
#' @keywords combinations
#' @export
calculate_Bliss <- function(sa1,
                            series_id1,
                            sa2,
                            series_id2,
                            metric,
                            measured_col = "smooth") {
  if (metric %in% c("GRvalue", "GR")) {
    lambda <- function(x, y) {
      ifelse(x < 0 | y < 0,
        pmin(x, y),
        2 ^ (log2(x + 1) * log2(y + 1)) - 1
        # formula for GR combination adapted from 
        # Holbeck et al. Cancer Res, vol.77(13), 2017 
        #   https://cancerres.aacrjournals.org/content/77/13/3564
        # growth rates are multiplicative, not GR values directly
      )
    }
  } else {
    lambda <- function(x, y) {
      ifelse(x < 0 | y < 0,
        pmin(x, y),
        x * y
        # Generalized Bliss formula for combination with potential 
        # negative values 
        #   ( Holbeck et al. Cancer Res, vol.77(13), 2017 
        #     https://cancerres.aacrjournals.org/content/77/13/3564 )
      )
    }
  }
  .calculate_matrix_metric(
    sa1, 
    series_id1, 
    sa2, 
    series_id2, 
    metric, 
    FXN = lambda,
    measured_col = measured_col
  )
}

#' @rdname calculate_matrix_metric
#' @keywords combinations
.calculate_matrix_metric <- function(sa1, 
                                     series_id1, 
                                     sa2, 
                                     series_id2, 
                                     metric,
                                     FXN,
                                     measured_col = "x") {
  checkmate::assert_data_table(sa1)
  checkmate::assert_data_table(sa2)

  checkmate::assert_true(all(sa1[[series_id2]] == 0L))
  checkmate::assert_true(all(sa2[[series_id1]] == 0L))
  
  assert_string(measured_col)

  data.table::setnames(sa1, measured_col, "metric1", skip_absent = TRUE)
  data.table::setnames(sa2, measured_col, "metric2", skip_absent = TRUE)

  u <- data.table::CJ(sa1[[series_id1]], sa2[[series_id2]])
  colnames(u) <- c(series_id1, series_id2)

  idx <- match(u[[series_id1]], sa1[[series_id1]])
  
  u <- u[sa1[, c(series_id1, "metric1"), with = FALSE], on = series_id1]
  u <- u[sa2[, c(series_id2, "metric2"), with = FALSE], on = series_id2]

  metric <- do.call(FXN, list(u$metric1, u$metric2))
  data.table::data.table(cbind(u, metric))
}
