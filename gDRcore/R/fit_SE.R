#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
fit_SE <- function(se, 
                   data_type = "single-agent",
                   nested_identifiers = NULL,
                   averaged_assay = "Averaged", 
                   metrics_assay = "Metrics", 
                   n_point_cutoff = 4,
                   range_conc = c(5e-3, 5),
                   force_fit = FALSE,
                   pcutoff = 0.05,
                   cap = 0.1,
                   curve_type = c("GRvalue", "RelativeViability")) {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_number(n_point_cutoff)
  checkmate::assert_numeric(range_conc)
  checkmate::assert_logical(force_fit)
  checkmate::assert_number(pcutoff)
  checkmate::assert_number(cap)
  checkmate::assert_character(curve_type)
  
  req_assays <- c(averaged_assay)
  lapply(req_assays, function(x) gDRutils::validate_se_assay_name(se, x))

  
  if (is.null(nested_identifiers)) {
    nested_identifiers <- get_default_nested_identifiers(se, data_model(data_type))
  }
  
  
  metric_cols <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
  avg_trt <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se, averaged_assay))
  iterator <- unique(avg_trt[, c("column", "row")])
  
  out <- gDRutils::loop(seq_len(nrow(iterator)), function(row) {
    x <- iterator[row, ]
    i <- x[["row"]]
    j <- x[["column"]]
 
    avg_df <- avg_trt[avg_trt$row == i & avg_trt$column == j, ]

    fit_df <- S4Vectors::DataFrame(matrix(NA, 2, length(metric_cols)))
    colnames(fit_df) <- metric_cols
    rownames(fit_df) <- c("RV", "GR")

    
    if (!is.null(avg_df) && all(dim(avg_df) > 0)) {
      if (!all(is.na(avg_df[[gDRutils::get_env_identifiers("concentration")]]))) {
        avg_df <- avg_df[avg_df[[
          gDRutils::get_env_identifiers("concentration")]] != 0, ]
      }
      fit_df <- S4Vectors::DataFrame(gDRutils::fit_curves(avg_df,
        series_identifiers = nested_identifiers,
        e_0 = 1,
        GR_0 = 1,
        n_point_cutoff = n_point_cutoff,
        range_conc = range_conc,
        force_fit = force_fit,
        pcutoff = pcutoff,
        cap = cap,
        normalization_type = curve_type))
    }

    if (nrow(fit_df) != 0L) {
      fit_df$row_id <- i
      fit_df$col_id <- j
    }
    fit_df
  })

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
