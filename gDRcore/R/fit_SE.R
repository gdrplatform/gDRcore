#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
fit_SE <- function(se, 
                   nested_identifiers = gDRutils::get_SE_identifiers(se, "concentration", simplify = TRUE),
                   averaged_assay = "Averaged", 
                   ref_GR_assay = "RefGRvalue",
                   ref_RV_assay = "RefRelativeViability",
                   metrics_assay = "Metrics", 
                   n_point_cutoff = 4,
                   range_conc = c(5e-3, 5),
                   force_fit = FALSE,
                   pcutoff = 0.05,
                   cap = 0.1,
                   curve_type = c("GR", "RV")) {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_number(n_point_cutoff)
  checkmate::assert_numeric(range_conc)
  checkmate::assert_logical(force_fit)
  checkmate::assert_number(pcutoff)
  checkmate::assert_number(cap)
  req_assays <- c(averaged_assay, ref_GR_assay, ref_RV_assay)
  present <- req_assays %in% SummarizedExperiment::assayNames(se)
  if (!all(present)) {
    stop(sprintf("unable to find required assays: '%s'", 
      paste0(req_assays[!present], collapse = ", ")))
  }

  metric_cols <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
  out <- vector("list", nrow(se) * ncol(se))
  avg_trt <- SummarizedExperiment::assay(se, averaged_assay)
  ref_GR <- SummarizedExperiment::assay(se, ref_GR_assay)
  ref_RV <- SummarizedExperiment::assay(se, ref_RV_assay)

  count <- 0
  for (i in seq_len(nrow(se))) {
    for (j in seq_len(ncol(se))) {
      count <- count + 1
      avg_df <- avg_trt[i, j][[1]]
      if (nrow(avg_df) == 0L) {
        next
        }
      
      fit_df <- S4Vectors::DataFrame(matrix(NA, 2, length(metric_cols)))
      colnames(fit_df) <- metric_cols
      rownames(fit_df) <- c("RV", "GR")

      
      if (!is.null(avg_df) && all(dim(avg_df) > 0)) {
        fit_df <- S4Vectors::DataFrame(gDRutils::fit_curves(avg_df,
          series_identifiers = nested_identifiers,
          e_0 = ref_RV[i, j],
          GR_0 = ref_GR[i, j],
          n_point_cutoff = n_point_cutoff,
          range_conc = range_conc,
          force_fit = force_fit,
          pcutoff = pcutoff,
          cap = cap,
          normalization_type = curve_type))
      }

      if (nrow(fit_df) != 0L) {
        fit_df$row_id <- rep(rownames(se)[i], nrow(fit_df))
        fit_df$col_id <- rep(colnames(se)[j], nrow(fit_df))
      }
      out[[count]] <- fit_df
    }
  }

  out <- S4Vectors::DataFrame(do.call("rbind", out))
  metrics <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c("row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assay(se, metrics_assay) <- metrics
  se <- gDRutils::set_SE_fit_parameters(se, 
    value = list(
      n_point_cutoff = n_point_cutoff,
      range_conc = range_conc,
      force_fit = force_fit,
      pcutoff = pcutoff,
      cap = cap)
  )
  se <- gDRutils::set_SE_processing_metadata(se,
                                             value = list(
                                               date_processed = Sys.Date(),
                                               session_info = utils::sessionInfo()))
  
  se
}
