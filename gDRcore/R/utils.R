#' @noRd
#' @keywords internal
#'
.create_mapping_factors <- function(rowdata, coldata) {
  mapping_cols <- c(colnames(rowdata), colnames(coldata))
  mapping_entries <- expand.grid(row_id = rownames(rowdata), col_id = rownames(coldata), stringsAsFactors = FALSE)

  mapping_entries <- base::merge(mapping_entries, rowdata, by.x = "row_id", by.y = 0, all.x = TRUE)
  mapping_entries <- base::merge(mapping_entries, coldata, by.x = "col_id", by.y = 0, all.x = TRUE)

  rownames(mapping_entries) <- seq_len(nrow(mapping_entries))
  mapping_entries
}


#' cleanup_metadata
#'
#' Cleanup a dataframe with metadata
#'
#' @param df_metadata a dataframe with metadata
#'
#' @return a dataframe with cleaned metadata
#'
#' @export
#'
cleanup_metadata <- function(df_metadata) {

  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))

  data.table::setDT(df_metadata)

  # clean up numberic fields
  df_metadata[[gDRutils::get_identifier("duration")]] <-
    round(as.numeric(df_metadata[[gDRutils::get_identifier("duration")]], 6))

  # identify potential numeric fields and replace NA by 0 - convert strings in factors
  for (c in setdiff(1:dim(df_metadata)[2], c(
    agrep(gDRutils::get_identifier("drug"), colnames(df_metadata)),
    agrep("Concentration", colnames(df_metadata)),
    grep(paste(
      c(
        gDRutils::get_identifier("cellline"),
        gDRutils::get_header("manifest"),
        gDRutils::get_identifier("well_position"),
        "compoundId"
      ),
      collapse = "|"
    ), colnames(df_metadata))
  ))) {
    vals <- unique(df_metadata[[c]])

    if (is.character(vals)) {
      num_vals <- as.numeric(vals)
      if (sum(is.na(num_vals)) > 2 || all(is.na(num_vals))) {
        df_metadata[[c]] <- as.character(df_metadata[[c]])
        futile.logger::flog.warn("Metadata field %s converted to strings",
                colnames(df_metadata)[c])
      } else {
        is.na(df_metadata[[c]]) <- 0
        df_metadata[[c]] <- as.numeric(df_metadata[[c]])
        futile.logger::flog.warn("Metadata field %s converted to numeric values",
                colnames(df_metadata)[c])
      }
    }
  }
    # TODO: specific to GNE database --> need to be replaced by a function
    df_metadata <- add_CellLine_annotation(df_metadata)
    # check that Gnumber_* are in the format 'G####' and add common name (or Vehicle or Untreated)

    for (i in agrep(gDRutils::get_identifier("drug"), colnames(df_metadata))) { # correct case issues
        for (w in gDRutils::get_identifier("untreated_tag")) {
            df_metadata[grep(w, df_metadata[[i]], ignore.case = TRUE), i] <- w
        }
    }
    # -----------------------

    df_metadata <- add_Drug_annotation(df_metadata)

    data.table::setDF(df_metadata)
    # clean up concentration fields
    for (i in agrep("Concentration", colnames(df_metadata))) {
        trt_n <- ifelse(regexpr("_\\d", colnames(df_metadata)[i]) > 0,
                            substr(colnames(df_metadata)[i], 15, 20), 1)
        DrugID_col <- ifelse(trt_n == 1, gDRutils::get_identifier("drug"), paste0(gDRutils::get_identifier("drug"), 
                                                                                  "_", trt_n))
        # set all untreated to 0
        df_metadata[df_metadata[, DrugID_col] %in% gDRutils::get_identifier("untreated_tag"), i] <- 0

        DrugID_0 <- setdiff(unique(df_metadata[df_metadata[, i] == 0, DrugID_col]), 
                            gDRutils::get_identifier("untreated_tag"))
        DrugID_0 <- DrugID_0[!is.na(DrugID_0)]
        if (length(DrugID_0) > 0) {
          futile.logger::flog.warn("Some concentration for %s are 0: %s",
                                   DrugID_col,
                                   paste(DrugID_0, collapse = " ; "))

        }
        df_metadata[, i] <- 10 ** round(log10(as.numeric(df_metadata[, i])), 6)
    }
  return(df_metadata)
}


#' Order_result_df
#'
#' Order a dataframe with results
#'
#' @param df_ a dataframe with results
#'
#' @return a ordered dataframe with results
#' @export
#'
Order_result_df <- function(df_) {

  # Assertions:
  stopifnot(inherits(df_, "data.frame"))

  cols <- c(gDRutils::get_header("ordered_1"),
            setdiff(colnames(df_),
                    c(
                      gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
                    )),
            gDRutils::get_header("ordered_2"))
  cols <- intersect(cols, colnames(df_))

  row_order_col <-
    intersect(
      c(
        gDRutils::get_identifier("cellline_name"),
        gDRutils::get_identifier("duration"),
        gDRutils::get_identifier("drugname"),
        "Concentration",
        paste0(c(
          paste0(gDRutils::get_identifier("drugname"), "_"), "Concentration_"
        ),
        sort(rep(2:10, 2))),
        setdiff(colnames(df_), c(
          gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
        ))
      ),
      cols
    )

  df_ <- df_[do.call(order, df_[, row_order_col]), cols]

  return(df_)
}


#' standardize_record_values
#'
#' map values to a dictionary
#'
#' @param x a named array
#' @param dictionary a named array
#'
#' @return a named array with updated names
#' @export
#'
standardize_record_values <- function(x, dictionary = DICTIONARY) {
  for (i in seq_len(length(dictionary))) {
    x[x == names(dictionary[i])] <- dictionary[[i]]
  }
  x
}
