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
        gDRutils::get_env_identifiers(c("cellline_name", "duration", "drug_name"), simplify = FALSE),
        "Concentration",
        paste0(c(
          paste0(gDRutils::get_env_identifiers("drug_name"), "_"), "Concentration_"
        ),
        sort(rep(2:10, 2))),
        setdiff(colnames(df_), c(
          gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
        ))
      ),
      cols
    )

  df_ <- df_[do.call(order, df_[, unlist(row_order_col)]), unlist(cols)]

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
  drug_ids <- unlist(gDRutils::get_env_identifiers(c("drug_name", "drug_name2"), simplify = FALSE))
  cl_id <- gDRutils::get_env_identifiers("cellline")
  conc2 <- gDRutils::get_env_identifiers("concentration2")
  if (all(.get_default_combo_identifiers() %in% colnames(df_))) {
    if (all(df_[[conc2]]
            %in% gDRutils::get_env_identifiers("untreated_tag"))) {
      "single-agent"
    } else {
      df_subset <- unique(subset(df_, select = c(drug_ids, cl_id, conc2)))
      cotrt <- all(table(subset(df_subset, select = c(drug_ids, cl_id))) < 4) # detect co-trt
      if (cotrt) {
        "single-agent"
      } else {
        "combo"
      }
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

#' Round concentration to ndigit significant digits
#'
#' @param x value to be rounded.
#' @param ndigit number of significant digits (default = 4).
#'
#' @return rounded x
#' @export
round_concentration <- function(x, ndigit = 3) {
  round(10 ^ (round(log10(x), ndigit)), ndigit - 1 - floor(log10(x)))
}


#' @keywords internal
#' @noRd
rbindList <- function(x) {
  S4Vectors::DataFrame(do.call(plyr::rbind.fill, x))
}

#' @keywords internal
#' @noRd
rbindParallelList <- function(x, name) {
  S4Vectors::DataFrame(do.call(plyr::rbind.fill, lapply(x, function(x) as.data.frame("[[" (x, name)))))
}

#' Function for detection of available cores used in parallel computing
#'
#' @return number of available cores
#' @export
detect_cores <- function() {
  x <- as.numeric(Sys.getenv("NUM_CORES"))
  if (is.na(x)) {
    x <- parallel::detectCores() - 1
  }
  x
}

#' Value Matching
#' 
#' Returns a lookup table or list of the positions of ALL matches of its first
#' argument in its second and vice versa. Similar to \code{\link{match}}, though
#' that function only returns the first match.
#' 
#' This behavior can be imitated by using joins to create lookup tables, but
#' \code{matches} is simpler and faster: usually faster than the best joins in
#' other packages and thousands of times faster than the built in
#' \code{\link{merge}}.
#' 
#' \code{all.x/all.y} correspond to the four types of database joins in the
#' following way:
#' 
#' \describe{ \item{left}{\code{all.x=TRUE}, \code{all.y=FALSE}} 
#' \item{right}{\code{all.x=FALSE}, \code{all.y=TRUE}} 
#' \item{inner}{\code{all.x=FALSE}, \code{all.y=FALSE}} 
#' \item{full}{\code{all.x=TRUE}, \code{all.y=TRUE}} }
#' 
#' Note that \code{NA} values will match other \code{NA} values.
#' 
#' @param x vector.  The values to be matched.  Long vectors are not currently
#'   supported.
#' @param y vector.  The values to be matched.  Long vectors are not currently
#'   supported.
#' @param all.x logical; if \code{TRUE}, then each value in \code{x} will be
#'   included even if it has no matching values in \code{y}
#' @param all.y logical; if \code{TRUE}, then each value in \code{y} will be
#'   included even if it has no matching values in \code{x}
#' @param list logical.  If \code{TRUE}, the result will be returned as a list
#'   of vectors, each vector being the matching values in y. If \code{FALSE},
#'   result is returned as a data frame with repeated values for each match.
#' @param indexes logical.  Whether to return the indices of the matches or the
#'   actual values.
#' @param nomatch the value to be returned in the case when no match is found.
#'   If not provided and \code{indexes=TRUE}, items with no match will be
#'   represented as \code{NA}.  If set to \code{NULL}, items with no match will
#'   be set to an index value of \code{length+1}.  If {indexes=FALSE}, they will
#'   default to \code{NA}.
#' @details Source of the function: https://github.com/cran/grr/blob/master/R/grr.R
#' @export
matches <- function(x, y, all.x = TRUE, all.y = TRUE, list = FALSE, indexes = TRUE, nomatch = NA) {
  require(gD)
  result <- .Call(getNativeSymbolInfo("matches"), x, y)
  result <- data.frame(x = result[[1]], y = result[[2]])
  if (!all.y) {
    result <- result[result$x != length(x) + 1, ]
  }
  if (!all.x) {
    result <- result[result$y != length(y) + 1, ]
  }
  if (!indexes) {
    result$x <- x[result$x]
    result$y <- y[result$y]
  } else if (!is.null(nomatch)) {
    result$x[result$x == length(x) + 1] <- nomatch
    result$y[result$y == length(y) + 1] <- nomatch
  }
  if (list) {
    result <- tapply(result$y, result$x, function(z) z[!is.na(z)])
  }
  result
}
