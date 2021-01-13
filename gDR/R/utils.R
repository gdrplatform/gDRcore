#' merge_data
#'
#' Merge all the input data into a single dataframe
#'
#' @param manifest a data frame with a manifest info
#' @param treatments a data frame with a treaatments info
#' @param data a data frame with a raw data info
#'
#' @return a dataframe with merged data
#' @export

merge_data <- function(manifest, treatments, data) {
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
  futile.logger::flog.info("Merging the metadata (manifest and treatment files")

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
  futile.logger::flog.warn("%i wells discarded for lack of annotation, %i data point selected",
                           dim(df_metadata_trimmed)[1],
                           sum(is.na(df_metadata[, gDRutils::get_identifier("drug")])))

  # clean up the metadata
  cleanedup_metadata <-
    cleanup_metadata(df_metadata_trimmed)
  stopifnot(dim(cleanedup_metadata)[1] == dim(df_metadata_trimmed)[1]) # should not happen

  df_merged <- merge(cleanedup_metadata, data, by = c("Barcode",
                                                      gDRutils::get_identifier("WellPosition")))
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
  futile.logger::flog.warn("%i well loaded, %i discarded for lack of annotation, %i data point selected",
      dim(data)[1],
      sum(is.na(df_merged[, gDRutils::get_identifier("drug")])),
      dim(df_raw_data)[1]
    )

  # reorder the columns
  df_raw_data <- Order_result_df(df_raw_data)

  return(df_raw_data)
}

#' identify_keys
#'
#' Identify keys in the DR data represented by dataframe or SummarizedExperiment or MultiAssayExperiment objects
#'
#' @param df_se_mae a dataframe or SummarizedExperiment or MultiassayExperiment with keys
#'
#' @return a list of keys
#' @export
#'

identify_keys <- function(df_se_mae) {

  # Assertions:
  stopifnot(inherits(df_se_mae, c("data.frame", "MultiAssayExperiment", "SummarizedExperiment")))


    if (any(class(df_se_mae) %in% c("MultiAssayExperiment", "SummarizedExperiment"))) {
        if ("MultiAssayExperiment" %in% class(df_se_mae)) {
            # if MAE, convert to SE based on the treated SE (could be optimized)
            df_se_mae <- df_se_mae[["treated"]]
            se_untrt <-  df_se_mae[["untreated"]]
        } else se_untrt <- NULL
        all_keys <- unique(c(
            colnames(SummarizedExperiment::rowData(df_se_mae)),
            colnames(SummarizedExperiment::colData(df_se_mae)),
            unlist(lapply(SummarizedExperiment::assay(df_se_mae), colnames))))
    } else { # case of a data frame
        all_keys <- colnames(df_se_mae)
    }

    keys <- list(Trt = setdiff(all_keys, "Barcode"),
            DoseResp = setdiff(all_keys,  "Barcode"),
            ref_Endpoint = setdiff(all_keys, c("Concentration",
                                            gDRutils::get_identifier("drug"),
                                            gDRutils::get_identifier("drugname"))),
            untrt_Endpoint = all_keys[ c(-agrep("Concentration", all_keys),
                                            -agrep(gDRutils::get_identifier("drug"), all_keys),
                                            -agrep(gDRutils::get_identifier("drugname"), all_keys))])
    keys[["Day0"]] <- setdiff(keys[["untrt_Endpoint"]], gDRutils::get_identifier("duration"))
    keys <- lapply(keys, function(x) setdiff(x, c(gDRutils::get_header("raw_data"),
        gDRutils::get_header("normalized_results"), "Template", gDRutils::get_identifier("WellPosition"), gDRutils::get_header("averaged_results"),
            gDRutils::get_header("metrics_results"), "ReferenceDivisionTime"
    )))
    keys <- lapply(keys, sort)

    # check if all values of a key is NA
    for (k in keys[["untrt_Endpoint"]]) {

        if ("SummarizedExperiment" %in% class(df_se_mae)) {
            # check the metadata fields for NA
            if (k %in% colnames(SummarizedExperiment::rowData(df_se_mae))) df_ <- SummarizedExperiment::rowData(df_se_mae)
            else if (k %in% colnames(SummarizedExperiment::colData(df_se_mae))) df_ <- SummarizedExperiment::colData(df_se_mae)
            else next # not a metadata

            if (all(is.na(df_[,k]))) keys <- lapply(keys, function(x) setdiff(x, k))

            if (!is.null(se_untrt) && k %in% colnames(SummarizedExperiment::rowData(se_untrt))) {
                df_ <- SummarizedExperiment::rowData(se_untrt)
                if (all(is.na(df_[df_[,gDRutils::get_identifier("duration")] == 0, k]))) {
                    keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
                }
            }
        } else { # case of a data frame
            if (all(is.na(df_se_mae[, k]))) {
                keys <- lapply(keys, function(x) setdiff(x, k))
            }
            if (all(is.na(df_se_mae[df_se_mae[,gDRutils::get_identifier("duration")] == 0, k]))) {
                keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
            }
        }
    }
  return(keys)
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
        gDRutils::get_identifier("WellPosition"),
        "compoundId"
      ),
      collapse = "|"
    ), colnames(df_metadata))
  ))) {
    vals <- unique(df_metadata[[c]])

    if (is.character(vals)) {
      num_vals <- as.numeric(vals)
      if (sum(is.na(num_vals)) > 2 || all(is.na(num_vals))) {
        df_metadata[[c]] <- factor(df_metadata[[c]])
        futile.logger::flog.warn("Metadata field %s converted to factors",
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
            df_metadata[grep(w, df_metadata[[i]], ignore.case = TRUE),i] <- w
        }
    }
    # -----------------------

    df_metadata <- add_Drug_annotation(df_metadata)
  
    data.table::setDF(df_metadata)
    # clean up concentration fields
    for (i in agrep("Concentration", colnames(df_metadata))) {
        trt_n <- ifelse(regexpr("_\\d", colnames(df_metadata)[i]) > 0,
                            substr(colnames(df_metadata)[i], 15, 20), 1)
        DrugID_col <- ifelse(trt_n == 1, gDRutils::get_identifier("drug"), paste0(gDRutils::get_identifier("drug"), "_", trt_n))
        df_metadata[df_metadata[,DrugID_col] %in% gDRutils::get_identifier("untreated_tag"), i] <- 0 # set all untreated to 0

        DrugID_0 <- setdiff(unique(df_metadata[ df_metadata[,i] == 0, DrugID_col]), gDRutils::get_identifier("untreated_tag"))
        DrugID_0 <- DrugID_0[!is.na(DrugID_0)]
        if (length(DrugID_0) > 0) {
          futile.logger::flog.warn("Some concentration for %s are 0: %s",
                                   DrugID_col,
                                   paste(DrugID_0, collapse = " ; "))

        }
        df_metadata[,i] <- 10 ** round(log10(as.numeric(df_metadata[, i])), 6)
        # df_metadata[,i] <- round(as.numeric(df_metadata[, i]), 10) # avoid mismatch due to string truncation
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

Order_result_df <- function (df_) {

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
        gDRutils::get_header("add_clid")[1],
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

standardize_record_values <- function(x, dictionary = DICTIONARY){
  for (i in 1:length(dictionary)) {
    x[x == names(dictionary[i])] <- dictionary[[i]]
  }
  x
}
