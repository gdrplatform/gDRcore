#' Transform mae into raw data
#'
#' The raw data uses the default identifiers, also for experiments built with
#' custom ones. The session identifiers are the same after the call as before.
#'
#' @param mae MultiAssayExperiment object with SummarizedExperiments containing "RawTreated" and "Controls" assays
#'
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_small")
#' convert_mae_to_raw_data(mae)
#'
#' @return data.table with raw data
#' @keywords convert_to_raw_data
#' @export
convert_mae_to_raw_data <- function(mae) {

  checkmate::assert_class(mae, "MultiAssayExperiment")
  checkmate::assert_true(all(vapply(gDRutils::MAEpply(mae, SummarizedExperiment::assayNames),
         function(x) all(c("RawTreated", "Controls") %in% x), FUN.VALUE = logical(1))))

  # the raw data comes back with the default identifiers (see convert_se_to_raw_data()),
  # so they are used for the whole call and the session is put back on exit
  session_identifiers <- gDRutils::get_env_identifiers()
  on.exit(restore_env_identifiers(session_identifiers), add = TRUE)
  gDRutils::reset_env_identifiers()

  # Get the data from se

  data <- gDRutils::MAEpply(mae, convert_se_to_raw_data)

  # Remove duplicates shared between assays to keep only original single-agent
  common_records <- Reduce(intersect, lapply(data, "[[", "record_id"))
  sa_name <- gDRutils::get_supported_experiments("sa")

  # Get combo exp name (also support extracting the obsolete name for reprocessing purposes)
  combo_name <- intersect(names(data), c(gDRutils::get_supported_experiments("combo"), "matrix"))
  # check if data contains combination data and shared single-agent records are unique (true for internal data)
  if (length(names(data)) > 1 && max(table(unique(data[[combo_name]])[record_id %in% common_records]$record_id)) == 1) {
    data[[sa_name]] <- data[[sa_name]][!record_id %in% common_records]
  } else {
    data[setdiff(names(data), sa_name)] <-
      lapply(data[setdiff(names(data), sa_name)], function(x) {
        x[!record_id %in% common_records]
      })
  }

  data_df <- data.table::rbindlist(data, fill = TRUE)

  data_df <- replace_NA_in_raw_data(data_df, mae)

  data.table::setorder(data_df)

  data_df <- data_df[!duplicated(data_df$record_id), ]
  selected_columns <- !names(data_df) %in% c("record_id", "BackgroundValue",
                                          "WellColumn", "WellRow", "Template", "swap_sa")
  data.table::setorder(data_df, record_id)
  data_df[, selected_columns, with = FALSE]
}


#' Transform se into raw_data
#'
#' The SE is first standardized, so the raw data uses the default identifiers
#' also when the SE was built with custom ones. The session identifiers are the
#' same after the call as before.
#'
#' @param se SummarizedExperiment object with "RawTreated" and "Controls" assays
#'
#' @return data.table with raw data
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' convert_se_to_raw_data(se)
#'
#' @keywords convert_to_raw_data
#' @export
convert_se_to_raw_data <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_true(all(c("RawTreated", "Controls") %in% SummarizedExperiment::assayNames(se)))

  # standardize_se() resets the session identifiers, so they are put back on exit
  session_identifiers <- gDRutils::get_env_identifiers()
  on.exit(restore_env_identifiers(session_identifiers), add = TRUE)
  # writing the default identifiers into the SE copy is the point here, so the
  # warning about overwriting its identifiers is not passed on
  se <- withCallingHandlers(gDRutils::standardize_se(se), warning = function(w) {
    if (grepl("overwriting existing metadata entry: 'identifiers'", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })

  # Get the data from se
  trt <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
  ctrl <- gDRutils::convert_se_assay_to_dt(se, "Controls")

  # Identifiers
  conc_cols1 <- intersect(unlist(gDRutils::get_env_identifiers("concentration",
                                                   simplify = FALSE)), names(trt))
  conc_cols2 <- intersect(unlist(gDRutils::get_env_identifiers("concentration2",
                                                    simplify = FALSE)), names(trt))
  conc_cols <- c(conc_cols1, conc_cols2)
  drug_cols1 <- intersect(unlist(gDRutils::get_env_identifiers(c("drug", "drug_name", "drug_moa"),
                                                   simplify = FALSE)), names(trt))
  drug_cols2 <- intersect(unlist(gDRutils::get_env_identifiers(c("drug2", "drug_name2", "drug_moa2"),
                                                   simplify = FALSE)), names(trt))
  drug_cols <- c(drug_cols1, drug_cols2)
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")[1]
  masked_tag <- gDRutils::get_env_identifiers("masked_tag")[1]
  ref_div_time <- gDRutils::get_env_identifiers("cellline_ref_div_time")
  readout_val <- "ReadoutValue"
  background_val <- "BackgroundValue"
  duration <- gDRutils::get_env_identifiers("duration")

  # Add required cols and correct the data
  ctrl[, (drug_cols) := NULL]
  ctrl[, (duration) := ifelse(!is.na(isDay0) & isDay0, 0, get(duration))]
  ctrl[, c("control_type", "isDay0") := NULL]
  ctrl[, c(eval(conc_cols), "BackgroundValue") := 0]
  ctrl[, eval(drug_cols) := untreated_tag]
  trt[get(conc_cols1) == 0, (drug_cols1) := untreated_tag]


  if (masked_tag %in% names(trt)) {
    ctrl[, eval(masked_tag) := FALSE]
  }

  if (length(conc_cols2) > 0 && length(drug_cols2) > 0) {
    trt[get(conc_cols2) == 0, (drug_cols2) := untreated_tag]
  }

  # Merge and adjust the data
  merged_df <- data.table::rbindlist(list(trt, ctrl), fill = TRUE)
  merged_df$rId <- merged_df$cId <- merged_df$ReadoutValue <- NULL
  data.table::setnames(merged_df, "CorrectedReadout", "ReadoutValue")

  exp_metadata <- gDRutils::get_SE_experiment_metadata(se)
  na.omit(merged_df, col = "ReadoutValue")
}


#' Set the session identifiers back to a saved list
#'
#' @param identifiers named list from \code{gDRutils::get_env_identifiers()}.
#'
#' @return \code{NULL}, invisibly.
#' @keywords internal
restore_env_identifiers <- function(identifiers) {
  gDRutils::reset_env_identifiers()
  for (k in names(identifiers)) {
    gDRutils::set_env_identifier(k, identifiers[[k]])
  }
  invisible(NULL)
}


#' @keywords internal
replace_NA_in_raw_data <- function(df, mae) {
  df_ <- data.table::copy(df)
  # Identifiers
  conc_cols <- intersect(unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2"),
                                                             simplify = FALSE)), names(df_))

  drug_cols <- intersect(unlist(gDRutils::get_env_identifiers(c("drug", "drug_name", "drug_moa",
                                                                         "drug2", "drug_name2", "drug_moa2"),
                                                             simplify = FALSE)), names(df_))

  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")

  df_[, (conc_cols) := lapply(.SD, function(x) replace(x, is.na(x), 0)), .SDcols = conc_cols]
  df_[, (drug_cols) := lapply(.SD, function(x) replace(x, is.na(x), untreated_tag[[1]])), .SDcols = drug_cols]
  df_
}
