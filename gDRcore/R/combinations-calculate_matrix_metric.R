#' Calculate a metric for combination data.
#'
#' Calculate a metric based off of single-agent values in combination screens.
#'
#' @param sa1 data.frame containing single agent data where entries in \code{series_id2} are all \code{0}.
#' Columns of the data.frame include identifiers and the \code{metric} of interest.
#' @param series_id1 String representing the column within \code{sa1} that represents id1.
#' @param sa2 data.frame containing single agent data where entries in \code{series_id1} are all \code{0}.
#' Columns of the data.frame include identifiers and the \code{metric} of interest.
#' @param series_id2 String representing the column within \code{sa2} that represents id2.
#' @param metric String of the column specifying the metric of interest. 
#' @param FXN Function to apply to the single-agent fits to calculate a metric.
#'
#' @return DataFrame containing a single row for every unique combination of the two series
#' identifiers and the corresponding calculated metric for each row.
#'
#' @name calculate_matrix_metric
#' @details
#' \code{calculate_HSA} takes the minimum of the two single agents readouts.
#' \code{calculate_Bliss} performs Bliss additivity calculation based on the single agent effects,
#' defined as \code{1-x} for the corresponding normalization.
#' See https://www.sciencedirect.com/science/article/pii/S1359644619303460?via%3Dihub#tb0005
#' for more details.
NULL


#' @rdname calculate_matrix_metric
#' @export
calculate_HSA <- function(sa1, series_id1, sa2, series_id2, metric) {
  .calculate_matrix_metric(sa1, series_id1, sa2, series_id2, metric, FXN = pmin)
}


#' @rdname calculate_matrix_metric
#' @export
calculate_Bliss <- function(sa1, series_id1, sa2, series_id2, metric) {
  if (metric %in% c("GRvalue", "GR")) {
    lambda <- function(x, y) {
      ifelse(x < 0 | y < 0,
        pmin(x, y),
        2 ^ (log2(x + 1) * log2(y + 1)) - 1
        # formula for GR combination adapted from Holbeck et al. Cancer Res, vol.77(13), 2017 
        #   https://cancerres.aacrjournals.org/content/77/13/3564
        # growth rates are multiplicative, not GR values directly
      )
    }
  } else {
    lambda <- function(x, y) {
      ifelse(x < 0 | y < 0,
        pmin(x, y),
        x * y
        # Generalized Bliss formula for combination with potential negative values 
        #   ( Holbeck et al. Cancer Res, vol.77(13), 2017 
        #     https://cancerres.aacrjournals.org/content/77/13/3564 )
      )
    }
  }
  .calculate_matrix_metric(sa1, series_id1, sa2, series_id2, metric, FXN = lambda)
}


#' @keywords internal
.calculate_matrix_metric <- function(sa1, series_id1, sa2, series_id2, metric, FXN) {
  checkmate::assert_true(all(sa1[[series_id2]] == 0L))
  checkmate::assert_true(all(sa2[[series_id1]] == 0L))

  colnames(sa1)[colnames(sa1) == metric] <- "metric1"
  colnames(sa2)[colnames(sa2) == metric] <- "metric2"

  # TODO: ensure they're unique?
  u <- expand.grid(sa1[[series_id1]], sa2[[series_id2]])
  colnames(u) <- c(series_id1, series_id2)

  idx <- match(u[[series_id1]], sa1[[series_id1]])
  u <- base::merge(u, sa1[, c(series_id1, "metric1")], by = series_id1)
  u <- base::merge(u, sa2[, c(series_id2, "metric2")], by = series_id2)

  metric <- do.call(FXN, list(u$metric1, u$metric2))
  cbind(u, metric)
}
