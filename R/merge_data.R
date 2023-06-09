#' merge_data
#'
#' Merge all the input data into a single data.table
#'
#' @param manifest a data.table with a manifest info
#' @param treatments a data.table with a treaatments info
#' @param data a data.table with a raw data info
#' 
#' @examples 
#' td <- gDRimport::get_test_data()
#' l_tbl <- gDRimport::load_data(
#'   manifest_file = gDRimport::manifest_path(td), 
#'   df_template_files = gDRimport::template_path(td), 
#'   results_file = gDRimport::result_path(td)
#' )
#' merge_data(
#'   l_tbl$manifest, 
#'   l_tbl$treatments, 
#'   l_tbl$data
#' )
#'
#' @return a data.table with merged data and metadata.
#' @export
#'
merge_data <- function(manifest, treatments, data) {
  # Assertions:
  checkmate::assert_data_table(manifest)
  checkmate::assert_data_table(treatments)
  checkmate::assert_data_table(data)
  
  futile.logger::flog.info("Merging data")
  
  # first unify capitalization in the headers of treatments with manifest
  duplicated_col <-
    setdiff(colnames(treatments)[toupper(colnames(treatments)) %in%
                                   toupper(colnames(manifest))],
            colnames(treatments)[colnames(treatments) %in% colnames(manifest)])
  for (m_col in duplicated_col) {
    colnames(treatments)[colnames(treatments) == m_col] <-
      colnames(manifest)[toupper(m_col) == toupper(colnames(manifest))]
    futile.logger::flog.trace(
      "Header %s in templates corrected to match case with manifest", 
      m_col
    )
  }

  # merge manifest and treatment files first
  identifiers <- gDRutils::validate_identifiers(manifest, req_ids = "barcode")
  df_metadata <- manifest[treatments, on = identifiers[["template"]], allow.cartesian = TRUE]
  
  futile.logger::flog.info(
    "Merging the metadata (manifest and treatment files)"
  )

  # sort out duplicate metadata columns
  duplicated_col <- setdiff(
    intersect(colnames(manifest), colnames(treatments)), 
    identifiers[["template"]]
  )
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
        any(
          df_metadata[, paste0(m_col, ".x")] != 
          df_metadata[, paste0(m_col, ".y")], 
          na.rm = TRUE
    )) {
      futile.logger::flog.warn(
        "Merge data: metadata field %s found in both the manifest
        and some templates with inconsistent values;
        values in template supersede the ones in the manifest", m_col
      )
    }
    df_metadata[, paste0(m_col, ".x")] <- NULL
    df_metadata[, paste0(m_col, ".y")] <- NULL
  }

  # check for the expected columns
  expected_headers <- identifiers[["cellline"]]
  headersOK <- expected_headers %in% colnames(df_metadata)
  if (any(!headersOK)) {
    headers <- paste(expected_headers[!headersOK], collpase = " ; ")
    stop(
      sprintf(
        "df_metadata does not contains all expected headers: %s required",
        headers
      )
    )
  }

  # remove wells not labeled
  drug_id <- identifiers[["drug"]]
  df_metadata_trimmed <-
    df_metadata[which(!is.na(df_metadata[, drug_id, with = FALSE])), ]
  futile.logger::flog.warn(
    "%i well loaded, %i wells discarded for lack of annotation, 
    %i data point selected\n",
    nrow(data),
    sum(is.na(df_metadata[, drug_id, with = FALSE])),
    nrow(df_metadata_trimmed)
  )
  
  # clean up the metadata
  cleanedup_metadata <- cleanup_metadata(df_metadata_trimmed)
  # should not happen
  stopifnot(nrow(cleanedup_metadata) == nrow(df_metadata_trimmed))
  
  data$WellColumn <- as.character(data$WellColumn)

  df_merged <- cleanedup_metadata[data, on = c(identifiers[["barcode"]],
                                               identifiers[["well_position"]])]
  if (nrow(df_merged) != nrow(data)) {
    # need to identify issue and output relevant warning
    futile.logger::flog.warn(
      "merge_data: Not all results have been matched with treatments;
      merged table is smaller than data table"
    )
  }
  if (nrow(df_merged) != nrow(df_metadata)) {
    # need to identify issue and print relevant warning
    futile.logger::flog.warn(
      "merge_data: Not all treatments have been matched with results;
      merged table is smaller than metadata table"
    )
  }

  # remove wells not labeled
  df_raw_data <-
    df_merged[which(!is.na(df_merged[, drug_id, with = FALSE])), ]
  
  # reorder the columns
  df_raw_data <- order_result_df(df_raw_data)

  return(df_raw_data)
}
