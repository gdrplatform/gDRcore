#' Transform mae into raw_data
#'
#' @param mae MultiAssayExperiment object with SummarizedExperiments containing "RawTreated" and "Controls" assays
#'
#' @return data.frame with raw data
#' @export
convert_mae_to_raw_data <- function(mae) {
  
  checkmate::assert_class(mae, "MultiAssayExperiment")
  lapply(MAEpply(mae, SummarizedExperiment::assayNames), function(x) all(c("RawTreated", "Controls") %in% x))
  
  # Get the data from se
  trt <- gDRutils::convert_mae_assay_to_dt(mae, "RawTreated")
  ctrl <- gDRutils::convert_mae_assay_to_dt(mae, "Controls")
  
  data.table::setDT(trt)
  data.table::setDT(ctrl)
  
  # Identifiers
  conc_cols <- intersect(unlist(get_SE_identifiers(mae[[1]], c("concentration", "concentration2"),
                                                   simplify = FALSE)), names(trt))
  drug_cols <- intersect(unlist(get_SE_identifiers(mae[[1]], c("drug", "drug_name", "drug_moa",
                                                         "drug2", "drug_name2", "drug_moa2"),
                                                   simplify = FALSE)), names(trt))
  untreated_tag <- get_SE_identifiers(mae[[1]], "untreated_tag")[1]
  masked_tag <- get_SE_identifiers(mae[[1]], "masked_tag")[1]
  
  # Add required cols and correct the data
  ctrl$Duration <- ifelse(!is.na(ctrl$isDay0) & ctrl$isDay0, 0, ctrl$Duration)
  ctrl[, c("control_type", "isDay0") := NULL]
  ctrl[, c(eval(conc_cols), "BackgroundValue") := 0]
  ctrl[, eval(drug_cols) := untreated_tag]
  ctrl[, eval(masked_tag) := FALSE]
  data.table::setnames(ctrl, "CorrectedReadout", "ReadoutValue")
  
  # Merge and adjust the data
  merged_df <- data.table::rbindlist(list(trt, unique(ctrl)), fill = TRUE)
  merged_df$rId <- merged_df$cId <-  merged_df$CorrectedReadout <- NULL
  
  as.data.frame(merged_df[!is.na(merged_df$ReadoutValue), ])
}


#' Transform se into raw_data
#'
#' @param se SummarizedExperiment object with "RawTreated" and "Controls" assays
#'
#' @return data.frame with raw data
#' @export
convert_se_to_raw_data <- function(se) {
  
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_true(all(c("RawTreated", "Controls") %in% SummarizedExperiment::assayNames(se)))
  
  # Get the data from se
  trt <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
  ctrl <- gDRutils::convert_se_assay_to_dt(se, "Controls")
  
  # Identifiers
  conc_cols <- intersect(unlist(get_SE_identifiers(se, c("concentration", "concentration2"),
                                                   simplify = FALSE)), names(trt))
  drug_cols <- intersect(unlist(get_SE_identifiers(se, c("drug", "drug_name", "drug_moa",
                                                         "drug2", "drug_name2", "drug_moa2"),
                                                   simplify = FALSE)), names(trt))
  untreated_tag <- get_SE_identifiers(se, "untreated_tag")[1]
  masked_tag <- get_SE_identifiers(se, "masked_tag")[1]
  
  # Add required cols and correct the data
  ctrl$Duration <- ifelse(!is.na(ctrl$isDay0) & ctrl$isDay0, 0, ctrl$Duration)
  ctrl[, c("control_type", "isDay0") := NULL]
  ctrl[, c(eval(conc_cols), "BackgroundValue") := 0]
  ctrl[, eval(drug_cols) := untreated_tag]
  ctrl[, eval(masked_tag) := FALSE]
  data.table::setnames(ctrl, "CorrectedReadout", "ReadoutValue")
  
  # Merge and adjust the data
  merged_df <- data.table::rbindlist(list(trt, unique(ctrl)), fill = TRUE)
  merged_df$rId <- merged_df$cId <-  merged_df$CorrectedReadout <- NULL
  
  as.data.frame(merged_df[!is.na(merged_df$ReadoutValue), ])
}
