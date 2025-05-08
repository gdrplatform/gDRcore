#' @rdname runDrugResponseProcessingPipelineFxns
#' 
#' @examples
#' 
#' d <- rep(seq(0.1, 0.9, 0.1), each = 4)
#' v <- rep(seq(0.1, 0.4, 0.1), 9)
#' df <- S4Vectors::DataFrame(
#'   Concentration = d,
#'   normalization_type = rep(c("GR", "RV"), length(v) * 2),
#'   x = rep(v, 2)
#' )
#' normalized <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = df)
#' 
#' keys <- list(Trt = "Concentration")
#' assays <- list("Normalized" = normalized)
#' se <- SummarizedExperiment::SummarizedExperiment(assays = assays)
#' se <- gDRutils::set_SE_keys(se, keys)
#' se <- gDRutils::set_SE_identifiers(se, gDRutils::get_env_identifiers())
#' se1 <- average_SE(
#'   se,
#'   data_type = "single-agent",
#'   normalized_assay = "Normalized",
#'   averaged_assay = "Averaged"
#' )
#' 
#' @export
#'
average_SE <- function(se,
                       data_type,
                       series_identifiers = NULL,
                       normalized_assay = "Normalized", 
                       averaged_assay = "Averaged") {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE)
  checkmate::assert_string(normalized_assay)
  checkmate::assert_string(averaged_assay)
  
  gDRutils::validate_se_assay_name(se, normalized_assay)

  if (is.null(series_identifiers)) {
    series_identifiers <- get_default_nested_identifiers(
      se, 
      data_model(data_type)
    )
  }
  
  trt_keys <- gDRutils::get_SE_keys(se, "Trt")

  checkmate::expect_character(trt_keys,
                              min.len = 1,
                              min.chars = 1,
                              info = "unexpected treated keys on 'se' object")
  
  gDRutils::apply_bumpy_function(se = se,
                                 FUN = average_FUN,
                                 req_assay_name = normalized_assay,
                                 out_assay_name = averaged_assay,
                                 parallelize = FALSE,
                                 series_identifiers = series_identifiers,
                                 trt_keys = trt_keys)
}


#' @keywords internal
average_FUN <- function(x,
                        series_identifiers = series_identifiers,
                        trt_keys = trt_keys) {
  
  x <- data.table::as.data.table(x)
  series_identifiers <- intersect(series_identifiers, colnames(x))
  agg_df <- x[, list(mean(x, na.rm = TRUE),
                     # set sd = 1 for single values
                     gDRutils::calc_sd(x)),
              by = c("normalization_type", series_identifiers)]
  data.table::setorderv(agg_df, c(series_identifiers, "normalization_type"))
  data.table::setnames(agg_df, c("V1", "V2"), c("x", "x_std"))
  agg_df <- S4Vectors::DataFrame(agg_df)
  
  rownames(agg_df) <- paste(
    cumsum(!duplicated(agg_df[, series_identifiers])),
    agg_df$normalization_type, 
    sep = "_"
  )
  
  if (nrow(agg_df) == 0) {
    p_trt_keys <- intersect(trt_keys, colnames(x))
    all_cols <- unique(c(series_identifiers, p_trt_keys, "x", "x_std",
                         "normalization_type"))
    norm_types <- unique(x$normalization_type)
    agg_df <- S4Vectors::DataFrame(matrix(NA, length(norm_types), length(all_cols)))
    colnames(agg_df) <- all_cols
    agg_df$normalization_type <- norm_types
    rownames(agg_df) <- paste(
      seq_len(nrow(agg_df)), 
      agg_df$normalization_type, 
      sep = "_"
    )
  }
  agg_df
}
