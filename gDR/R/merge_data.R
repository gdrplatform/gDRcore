#' merge_data
#'
#' Merge all the input data into a single data.frame
#'
#' @param manifest a data frame with a manifest info
#' @param treatments a data frame with a treaatments info
#' @param data a data frame with a raw data info
#' @param collapse_Drugs mode drug_2 to Drug when Drug == vehicle (default is TRUE)
#'
#' @return a dataframe with merged data and metadata.
#' @export
#'
merge_data <- function(manifest, treatments, data, collapse_Drugs = TRUE) {
  # Assertions:
  stopifnot(inherits(manifest, "data.frame"))
  stopifnot(inherits(treatments, "data.frame"))
  stopifnot(inherits(data, "data.frame"))

  futile.logger::flog.info("Merging data")

  # first unify capitalization in the headers of treatments with manifest
  duplicated_col <-
    setdiff(colnames(treatments)[toupper(colnames(treatments)) %in%
                                   toupper(colnames(manifest))],
            colnames(treatments)[colnames(treatments) %in% colnames(manifest)])
  for (m_col in duplicated_col) {
    colnames(treatments)[colnames(treatments) == m_col] <-
      colnames(manifest)[toupper(m_col) == toupper(colnames(manifest))]
    futile.logger::flog.trace("Header %s in templates corrected to match case with manifest", m_col)
  }

  # merge manifest and treatment files first
  df_metadata <- merge(manifest, treatments, by = "Template")
  futile.logger::flog.info("Merging the metadata (manifest and treatment files)")

  # sort out duplicate metadata columns
  duplicated_col <-
    setdiff(intersect(colnames(manifest), colnames(treatments)), "Template")
  for (m_col in duplicated_col) {
    df_metadata[, m_col] <-
      df_metadata[, paste0(m_col, ".y")] # parse template values
    missing_idx <-
      is.na(df_metadata[, m_col]) | df_metadata[, m_col] %in% c("", "-")

    # add manifest values when missing in template
    df_metadata[missing_idx, m_col] <-
      df_metadata[missing_idx, paste0(m_col, ".x")]

    # check for conflicts
    double_idx <- !(is.na(df_metadata[, paste0(m_col, ".x")]) |
                      df_metadata[, paste0(m_col, ".x")] %in% c("", "-")) &
      !(is.na(df_metadata[, paste0(m_col, ".y")]) |
          df_metadata[, paste0(m_col, ".y")] %in% c("", "-"))
    if (any(double_idx) &&
        any(df_metadata[, paste0(m_col, ".x")] != df_metadata[, paste0(m_col, ".y")], na.rm = TRUE)) {
      futile.logger::flog.warn("Merge data: metadata field %s found in both the manifest
                               and some templates with inconsistent values;
                               values in template supersede the ones in the manifest", m_col)
    }
    df_metadata[, paste0(m_col, ".x")] <- NULL
    df_metadata[, paste0(m_col, ".y")] <- NULL
  }

  # check for the expected columns
  expected_headers <- gDRutils::get_identifier("cellline")
  headersOK <- expected_headers %in% colnames(df_metadata)
  if (any(!headersOK)) {
    stop(sprintf(
      "df_metadata does not contains all expected headers: %s required",
      paste(expected_headers[!(expected_headers %in% col_df)], collpase = " ; ")
    ))
  }

  # remove wells not labeled
  df_metadata_trimmed <-
    df_metadata[!is.na(df_metadata[, gDRutils::get_identifier("drug")]),]
  futile.logger::flog.warn("%i well loaded, %i wells discarded for lack of annotation, %i data point selected\n",
                           dim(data)[1],
                           sum(is.na(df_metadata[, gDRutils::get_identifier("drug")])),
                           dim(df_metadata_trimmed)[1]
                           )

  # clean up the metadata
  cleanedup_metadata <-
    cleanup_metadata(df_metadata_trimmed)
  stopifnot(dim(cleanedup_metadata)[1] == dim(df_metadata_trimmed)[1]) # should not happen

  df_merged <- merge(cleanedup_metadata, data, by = c("Barcode",
                                                      gDRutils::get_identifier("well_position")))
  if (dim(df_merged)[1] != dim(data)[1]) {
    # need to identify issue and output relevant warning
    futile.logger::flog.warn("merge_data: Not all results have been matched with treatments;
                             merged table is smaller than data table")
  }
  if (dim(df_merged)[1] != dim(df_metadata)[1]) {
    # need to identify issue and print relevant warning
    futile.logger::flog.warn("merge_data: Not all treatments have been matched with results;
                             merged table is smaller than metadata table")
  }

  # remove wells not labeled
  df_raw_data <-
    df_merged[!is.na(df_merged[, gDRutils::get_identifier("drug")]), ]
    
  if (collapse_Drugs && paste0(gDRutils::get_identifier("drug"), '_2') %in% colnames(df_raw_data)) {
        # Secondary drug for some combinations can primary drug when viewed as single-agent (conc1 = 0), so swap. 
        # swap the data related to Drug and Drug_2 such that it can considered as a single-agent condition
        swap_idx <- df_raw_data$Concentration_2 > 0 &
                      df_raw_data[,gDRutils::get_identifier("drug")] %in% gDRutils::get_identifier("untreated_tag")
        temp_df <- df_raw_data[swap_idx, ]
        header_names = c(gDRutils::get_identifier("drug"), gDRutils::get_identifier("drugname"), 
                          gDRutils::get_identifier("drug_moa"), 'Concentration')
        temp_df[, c(header_names, paste0(header_names, '_2'))] <-
          temp_df[, c(paste0(header_names, '_2'), header_names)]

        futile.logger::flog.warn("merge_data: swapping Drug and Drug_2 for %i rows because Concentration == 0 and Concentration_2 > 0",
          sum(swap_idx))

        df_raw_data <- rbind(df_raw_data[!swap_idx, ], temp_df)
  }

  # reorder the columns
  df_raw_data <- Order_result_df(df_raw_data)

  return(df_raw_data)
}
