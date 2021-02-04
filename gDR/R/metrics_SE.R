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
#' @param se a SummarizedExperiment with averaged and normalized assays
#' @param studyConcThresh a numeric with study concentration threshold (4 by default)
#'
#' @return a \code{SummarizedExperiment} with additional assay with metrics
#'
#' @export
#'
fit_SE <- function(se, averaged_ctrl_assay = "AvgControls", averaged_assay = "Averaged", 
  metrics_assay = "Metrics", studyConcThresh = 4) {

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
  avg_trt <- SummarizedExperiment::assay(se, "Averaged")
  avg_ctrl <- SummarizedExperiment::assay(se, "Avg_Controls")
  for (i in rownames(se)) {
      for (j in colnames(se)) {
	  df_ <- avg_trt[[i, j]]
	  if (!is.null(df_) && all(dim(df_) > 0)) { # studyConcThresh is embeded in RVGRfits
	      metrics[[i, j]] <- DataFrame(gDRutils::fit_curves(df_,
		  e_0 = avg_ctrl[[i, j]]$RefRelativeViability,
		  GR_0 = avg_ctrl[[i, j]]$RefGRvalue,
		  n_point_cutoff = studyConcThresh))
	  } else {
	      out <- DataFrame(matrix(NA, 0, length(gDRutils::get_header("response_metrics")) + 2))
	      colnames(out) <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
	      metrics[[i, j]] <- out
	  }
      }
  }
  SummarizedExperiment::assay(se, "Metrics") <- metrics
  return(se)
}
