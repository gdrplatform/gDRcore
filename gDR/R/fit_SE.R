#' metrics_SE
#' @export
#'
metrics_SE <- function(...) {
  .Deprecated("fit_SE", package="gDR")
  fit_SE(...)
}


#' fit_SE
#'
#' Calculate metrics for DR data
#'
#' @param avgSE a SummarizedExperiment with averaged and normalized assays
#' @param studyConcThresh a numeric with study concentration threshold (4 by default)
#'
#' @return a SummarizedExperiment with additional assay with metrics
#'
#' @export
#'
fit_SE <- function(avgSE, studyConcThresh = 4) {
    .Deprecated(msg = "see fit_SE2 for similar, but not identical functionality")

    # Assertions:
    checkmate::assert_class(avgSE, "SummarizedExperiment")
    checkmate::assert_number(studyConcThresh)

    stopifnot(is.numeric(studyConcThresh))
    # this is not used as we enforce the same conditions as the input SE; not collapsing allowed
    # if (is.null(DoseRespKeys)) {
    #     if ("Keys" %in% names(metadata(avgSE))) DoseResp = metadata(avgSE)$Keys$DoseResp
    #     else DoseRespKeys = identify_keys(avgSE)$DoseResp
    # } else {
    #     metadata(avgSE)$Keys$DoseResp = DoseRespKeys
    # }

    metricsSE <- avgSE
    SummarizedExperiment::assay(metricsSE, "Metrics") <- SummarizedExperiment::assay(metricsSE, "Averaged")

    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in a foor loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    mSE_m <- SummarizedExperiment::assay(metricsSE, "Metrics")
    a_SE <- SummarizedExperiment::assay(metricsSE, "Averaged")
    aCtrl_SE <- SummarizedExperiment::assay(metricsSE, "Avg_Controls")
    for (i in rownames(metricsSE)) {
        for (j in colnames(metricsSE)) {
            df_ <- a_SE[[i, j]]
            if (!is.null(df_) && all(dim(df_) > 0)) { # studyConcThresh is embeded in RVGRfits
                mSE_m[[i, j]] <- DataFrame(gDRutils::fit_curves(df_,
                    e_0 = aCtrl_SE[[i, j]]$RefRelativeViability,
                    GR_0 = aCtrl_SE[[i, j]]$RefGRvalue,
                    n_point_cutoff = studyConcThresh))
            } else {
                out <- DataFrame(matrix(NA, 0, length(gDRutils::get_header("response_metrics"))+2))
                colnames(out) <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
                mSE_m[[i, j]] <- out
            }
        }
    }
    SummarizedExperiment::assay(metricsSE, "Metrics") <- mSE_m
    return(metricsSE)
}


#' fit_SE2
#'
#' Fit curves and obtain fit metrics from normalized, averaged drug response data.
#'
#' @param se a \linkS4class{SummarizedExperiment} with averaged and normalized assays.
#' Corresponding assay names are specified by \code{averaged_assay}, \code{ref_GR_assay}, and \code{ref_RV_assy}.
#' @param averaged_assay string of the name of the averaged assay in the \linkS4class{SummarizedExperiment}.
#' @param ref_GR_assay string of the name of the reference GR assay in the \linkS4class{SummarizedExperiment}.
#' @param ref_RV_assay string of the name of the reference Relative Viability assay in the \linkS4class{SummarizedExperiment}.
#' @param metrics_assay string of the name of the metrics assay to output.
#' @param ndigit_rounding integer indicating number of digits to round to in calculations.
#' Defaults to \code{4}. 
#'
#' @return the original \linkS4class{SummarizedExperiment} with an additional assay named \code{metrics_assay}.
#'
#' @seealso runDrugResponseProcessingPipeline2
#' @export
#'
fit_SE2 <- function(se, 
                    averaged_assay = "Averaged", 
                    ref_GR_assay = "RefGRvalue",
                    ref_RV_assay = "RefRelativeViability",
                    metrics_assay = "Metrics", 
                    ndigit_rounding = 4) {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_number(ndigit_rounding)

  metric_cols <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
  out <- vector("list", nrow(se) * ncol(se))
  avg_trt <- SummarizedExperiment::assay(se, averaged_assay)
  ref_GR <- SummarizedExperiment::assay(se, ref_GR_assay)
  ref_RV <- SummarizedExperiment::assay(se, ref_RV_assay)
  count <- 1
  for (i in seq_len(nrow(se))) {
    for (j in seq_len(ncol(se))) {
      avg_df <- avg_trt[i, j][[1]]
      fit_df <- S4Vectors::DataFrame(matrix(NA, 0, length(metric_cols)))
      colnames(fit_df) <- metric_cols

      if (!is.null(avg_df) && all(dim(avg_df) > 0)) {
	fit_df <- S4Vectors::DataFrame(gDRutils::fit_curves(avg_df,
	  e_0 = ref_RV[i, j],
	  GR_0 = ref_GR[i, j],
	  n_point_cutoff = ndigit_rounding))
      }

      if (nrow(fit_df) != 0L) {
	fit_df$row_id <- rep(rownames(se)[i], nrow(fit_df))
	fit_df$col_id <- rep(colnames(se)[j], nrow(fit_df))
      }
      out[[count]] <- fit_df
      count <- count + 1
    }
  }

  out <- S4Vectors::DataFrame(do.call("rbind", out))
  metrics <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c("row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assay(se, metrics_assay) <- metrics
  return(se)
}
