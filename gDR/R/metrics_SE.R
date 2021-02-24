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
#' Calculate metrics for DR data
#'
#' @param se a SummarizedExperiment with averaged and normalized assays
#' @param studyConcThresh a numeric with study concentration threshold (4 by default)
#'
#' @return a \code{SummarizedExperiment} with additional assay with metrics
#'
#' @export
#'
fit_SE2 <- function(se, 
                    averaged_assay = "Averaged", 
                    metrics_assay = "Metrics", 
                    studyConcThresh = 4) {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_number(studyConcThresh)

  stopifnot(is.numeric(studyConcThresh))
  # this is not used as we enforce the same conditions as the input SE; not collapsing allowed
  # if (is.null(DoseRespKeys)) {
  #     if ("Keys" %in% names(metadata(se))) DoseResp = metadata(se)$Keys$DoseResp
  #     else DoseRespKeys = identify_keys(se)$DoseResp
  # } else {
  #     metadata(se)$Keys$DoseResp = DoseRespKeys
  # }

  # TODO: create bumpy matrix? 
  # metrics <- create_empty_bm()
  avg_trt <- SummarizedExperiment::assay(se, averaged_treated_assay)
  fit_ref <- SummarizedExperiment::assay(se, fitting_ref_assay)
  for (i in rownames(se)) {
      for (j in colnames(se)) {
	  df_ <- avg_trt[[i, j]]
	  if (!is.null(df_) && all(dim(df_) > 0)) { # studyConcThresh is embeded in RVGRfits
	      metrics[[i, j]] <- DataFrame(gDRutils::fit_curves(df_,
		  e_0 = fit_ref[[i, j]]$RefRelativeViability,
		  GR_0 = fit_ref[[i, j]]$RefGRvalue,
		  n_point_cutoff = studyConcThresh))
	  } else {
	      out <- DataFrame(matrix(NA, 0, length(gDRutils::get_header("response_metrics")) + 2))
	      colnames(out) <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
	      metrics[[i, j]] <- out
	  }
      }
  }

  SummarizedExperiment::assay(se, metrics_assay) <- metrics
  return(se)
}
