#' @rdname runDrugResponseProcessingPipelineFxns
#' 
#' @examples
#' 
#' d <- rep(seq(0.1, 0.9, 0.1), each = 4)
#' v <- rep(seq(0.1, 0.4, 0.1), 9)
#' df <- S4Vectors::DataFrame(Concentration = d,
#'                            masked = rep(c(TRUE, TRUE, TRUE, FALSE), 9),
#'                            normalization_type = rep(c("GR", "RV"),
#'                                                     length(v) * 2),
#'                            x = rep(v, 2))
#' normalized <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = df)
#' 
#' keys <- list(Trt = "Concentration",
#'              "masked_tag" = "masked")
#' assays <- list("Normalized" = normalized)
#' se <- SummarizedExperiment::SummarizedExperiment(assays = assays)
#' se <- gDRutils::set_SE_keys(se, keys)
#' se <- gDRutils::set_SE_identifiers(se, gDRutils::get_env_identifiers())
#' .
#' se1 <-
#'   average_SE(
#'     se,
#'     data_type = "single-agent",
#'     override_masked = FALSE,
#'     normalized_assay = "Normalized",
#'     averaged_assay = "Averaged"
#'   )
#' 
#' 
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
      unmasked <- norm_df[!masked, , drop = FALSE]
      avg_df <- stats::aggregate(stats::as.formula(paste("x ~ ",
                                                  paste(c("normalization_type", series_identifiers),
                                                        collapse = " + "))),
                               function(x) mean(x, na.rm = TRUE),
                               data = unmasked, na.action = stats::na.pass)
      std_df <- stats::aggregate(stats::as.formula(paste("x ~ ",
                                                  paste(c("normalization_type", series_identifiers),
                                                        collapse = " + "))),
                               function(x) stats::sd(x, na.rm = TRUE),
                               data = unmasked, na.action = stats::na.pass)
    
      colnames(std_df)[colnames(std_df) == "x"] <- "x_std"
      agg_df <- S4Vectors::DataFrame(merge(avg_df, std_df,
                                         by = c(series_identifiers, "normalization_type"),
                                         sort = FALSE))
      rownames(agg_df) <- paste(cumsum(!duplicated(agg_df[, series_identifiers])),
                                agg_df$normalization_type, sep = "_")
    }

    if (all(masked) || nrow(agg_df) == 0) {
      p_trt_keys <- intersect(trt_keys, colnames(norm_df))
      all_cols <- unique(c(series_identifiers, p_trt_keys, "x", "x_std",
                           "normalization_type", "row_id", "col_id"))
      agg_df <- S4Vectors::DataFrame(matrix(NA, 1, length(all_cols)))
      colnames(agg_df) <- all_cols
      rownames(agg_df) <- paste(seq_len(nrow(agg_df)), agg_df$normalization_type, sep = "_")
    }
    agg_df$row_id <- i
    agg_df$col_id <- j
    agg_df
  })

  out <- S4Vectors::DataFrame(do.call("rbind", out))
  
  averaged <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c(masked_tag_key, "row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assays(se)[[averaged_assay]] <- averaged
  se
}
