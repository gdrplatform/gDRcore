#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
average_SE <- function(se,
                       data_type,
                       series_identifiers = NULL,
                       override_masked = FALSE,
                       normalized_assay = "Normalized", 
                       averaged_assay = "Averaged") {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE)
  checkmate::assert_flag(override_masked)
  checkmate::assert_string(normalized_assay)
  checkmate::assert_string(averaged_assay)
  
  gDRutils::validate_se_assay_name(se, normalized_assay)

  if (is.null(series_identifiers)) {
    series_identifiers <- get_default_nested_identifiers(se, data_model(data_type))
  }

  trt_keys <- gDRutils::get_SE_keys(se, "Trt")
  masked_tag_key <- gDRutils::get_SE_keys(se, "masked_tag")

  checkmate::expect_character(trt_keys,
                              min.len = 1,
                              min.chars = 1,
                              info = "unexpected treated keys on 'se' object")
  checkmate::expect_character(
    masked_tag_key,
    min.len = 1,
    min.chars = 1,
    info = "unexpected masked_tag_key keys on 'se' object"
  )
  
  normalized <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se, normalized_assay))

  std_cols <- c("GRvalue", "RelativeViability")
  iterator <- unique(normalized[, c("column", "row")])
  
  out <- gDRutils::loop(seq_len(nrow(iterator)), function(row) {
    x <- iterator[row, ]
    i <- x[["row"]]
    j <- x[["column"]]
    norm_df <- normalized[normalized$row == i & normalized$column == j, ]
    # bypass 'masked' filter
    masked <- norm_df[[masked_tag_key]] & !override_masked
    if (sum(!masked) > 0) {
      series_identifiers <- intersect(series_identifiers, colnames(norm_df))
      if (length(missing <- setdiff(std_cols, colnames(norm_df))) > 0L) {
        stop(sprintf("missing expected columns in nested normalized dataframe: '%s'", 
          paste0(missing, collapse = ", ")))
      }

      unmasked <- norm_df[!masked, , drop = FALSE]
      avg_df <- stats::aggregate(unmasked[, std_cols],
        by = as.list(unmasked[, series_identifiers, drop = FALSE]), 
        function(x) mean(x, na.rm = TRUE))

      std_df <- stats::aggregate(unmasked[, std_cols],
        by = as.list(unmasked[, series_identifiers, drop = FALSE]), 
        function(x) stats::sd(x, na.rm = TRUE))
      colnames(std_df)[colnames(std_df) %in% std_cols] <-
        paste0("std_", colnames(std_df)[colnames(std_df) %in% std_cols])

      agg_df <- S4Vectors::DataFrame(merge(avg_df, std_df, by = series_identifiers, sort = FALSE))
    } else {
      # <= 1L unmasked values
      p_trt_keys <- intersect(trt_keys, colnames(norm_df))
      all_cols <- unique(c(series_identifiers, std_cols, p_trt_keys, paste0("std_", std_cols), "row_id", "col_id"))
      agg_df <- S4Vectors::DataFrame(matrix(NA, 1, length(all_cols)))
      colnames(agg_df) <- all_cols
    }

    if (nrow(agg_df) != 0L) {
      agg_df$row_id <- i
      agg_df$col_id <- j
    }
    agg_df
  })

  out <- S4Vectors::DataFrame(do.call("rbind", out))
  
  averaged <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c(masked_tag_key, "row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assays(se)[[averaged_assay]] <- averaged
  se
}
