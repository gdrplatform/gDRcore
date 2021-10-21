#' @noRd
#' @keywords internal
#'
.create_mapping_factors <- function(rowdata, coldata) {
  mapping_cols <- c(colnames(rowdata), colnames(coldata))
  mapping_entries <- expand.grid(row_id = rownames(rowdata), col_id = rownames(coldata), stringsAsFactors = FALSE)

  mapping_entries <- base::merge(mapping_entries, rowdata, by.x = "row_id", by.y = 0, all.x = TRUE)
  mapping_entries <- base::merge(mapping_entries, coldata, by.x = "col_id", by.y = 0, all.x = TRUE)

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
  df_metadata[[gDRutils::get_env_identifiers("duration")]] <-
    round(as.numeric(df_metadata[[gDRutils::get_env_identifiers("duration")]], 6))

  # identify potential numeric fields and replace NA by 0 - convert strings in factors
  for (c in setdiff(1:dim(df_metadata)[2], c(
    agrep(gDRutils::get_env_identifiers("drug"), colnames(df_metadata)),
    agrep("Concentration", colnames(df_metadata)),
    grep(paste(
      c(
        gDRutils::get_env_identifiers("cellline"),
        gDRutils::get_header("manifest"),
        gDRutils::get_env_identifiers("well_position"),
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

    for (i in agrep(gDRutils::get_env_identifiers("drug"), colnames(df_metadata))) { # correct case issues
        for (w in gDRutils::get_env_identifiers("untreated_tag")) {
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
        DrugID_col <- ifelse(trt_n == 1, gDRutils::get_env_identifiers("drug"),
                             paste0(gDRutils::get_env_identifiers("drug"), "_", trt_n))
        df_metadata[df_metadata[, DrugID_col] %in%
                      gDRutils::get_env_identifiers("untreated_tag"), i] <- 0 # set all untreated to 0

        DrugID_0 <- setdiff(unique(df_metadata[df_metadata[, i] == 0, DrugID_col]),
                            gDRutils::get_env_identifiers("untreated_tag"))
        DrugID_0 <- DrugID_0[!is.na(DrugID_0)]
        if (length(DrugID_0) > 0) {
          futile.logger::flog.warn("Some concentration for %s are 0: %s",
                                   DrugID_col,
                                   paste(DrugID_0, collapse = " ; "))

        }
        df_metadata[, i] <- 10 ^ round(log10(as.numeric(df_metadata[, i])), 6)
        # df_metadata[,i] <- round(as.numeric(df_metadata[, i]), 10) # avoid mismatch due to string truncation # nolint
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
        gDRutils::get_env_identifiers(c("cellline_name", "duration", "drugname"), simplify = FALSE),
        "Concentration",
        paste0(c(
          paste0(gDRutils::get_env_identifiers("drugname"), "_"), "Concentration_"
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


#' Detect model of data
#'
#' @param df_ data.frame of raw drug response data containing both treated and untreated values. 
#'
#' @return string with the information of the raw data follows single-agent or combo data model
#' @export
data_model <- function(df_) {
  checkmate::assert_data_frame(df_)
  if (all(.get_default_combo_identifiers() %in% colnames(df_))) {
    if (all(df_[[gDRutils::get_env_identifiers("concentration2")]]
            %in% gDRutils::get_env_identifiers("untreated_tag"))) {
      "single-agent"
    } else {
      "combo"
    }
  } else if (.get_default_single_agent_identifiers() %in% colnames(df_)) {
    "single-agent"
  } else {
    stop("Unknown data model")
  }
}

#' Get default nested identifiers
#'
#' @param x data.frame with raw data or SummarizedExperiment object with gDR assays
#' @param assayName assay name used for finding nested_identifiers in SummarizedExperiment object
#' @return vector of nested identifiers
#' @export
get_nested_default_identifiers <- function(x, ...) {
  UseMethod("get_nested_default_identifiers")
}


#' @export
#' @describeIn get_nested_default_identifiers
get_nested_default_identifiers.data.frame <- function(x) {
  checkmate::assert_data_frame(x)
  data_type <- data_model(x)
  if (data_type == "single-agent") {
    .get_default_single_agent_identifiers()
  } else {
    .get_default_combo_identifiers()
  }
}

#' @export
#' @describeIn get_nested_default_identifiers
get_nested_default_identifiers.SummarizedExperiment <- function(x,
                                                                assayName =
                                                                  tail(SummarizedExperiment::assayNames(x), 1)) {
  checkmate::assert_class(x, "SummarizedExperiment")
  intersect(.get_default_combo_identifiers(),
            names(BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(x, assayName))))
}

.get_default_single_agent_identifiers <- function() {
  gDRutils::get_env_identifiers("concentration")
}

.get_default_combo_identifiers <- function() {
  unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2"),
                                       simplify = FALSE))
}

.catch_warnings <- function(x) {
  warn  <- unlist(unique(x$warning))
  if (!is.null(warn)) {
    futile.logger::flog.warn(paste0(warn, collapse = "\n"))
  }
}

.round_concentration <- function(x) {
  10 ^ (round(log10(x), 3))
}

#' @keywords internal
#' @noRd
rbindList <- function(x) {
  S4Vectors::DataFrame(do.call(plyr::rbind.fill, x))
}
