#' df_to_bm_assay
#'
#' Convert data.frame with dose-reponse data into a BumpyMatrix assay.
#'
#' The 'assay' is simply a BumpyMatrix object with rownames being the treatment ids,
#' colnames being the ids of the cell lines and values with
#' dose-response data for given cell lines under given conditions.
#'
#' @param data tibble or data.frame with drug-response data
#' @param data_type string type of data to be returned: all, for untreated conditions only ('untreated')
#' or for treated conditions only ('treated')
#' @param discard_keys a vector of keys that should be discarded
#'
#' @return BumpyMatrix object
#'
#' @export
df_to_bm_assay <-
  function(data,
           data_type = c("all", "treated", "untreated"),
           discard_keys = NULL) {
    # Assertions:
    stopifnot(any(inherits(data, "data.frame"), checkmate::test_character(data), inherits(data, "DataFrame")))
    checkmate::assert_character(data_type)
    data_type <- match.arg(data_type)
    checkmate::assert_character(discard_keys, null.ok = TRUE)
    ####
    data <- as(data, "DataFrame")
    allMetadata <- gDR::getMetaData(data, discard_keys = discard_keys)

    seColData <- allMetadata$colData
    cl_entries <- setdiff(colnames(seColData), c("col_id", "name_"))
    seRowData <- allMetadata$rowData
    cond_entries <-
      setdiff(colnames(seRowData), c("row_id", "name_"))
    dataCols <- allMetadata$dataCols

    complete <-
      S4Vectors::DataFrame(
        expand.grid(
          row_id = seRowData$row_id,
          col_id = seColData$col_id,
          stringsAsFactors = FALSE
        )
      )
    complete <- merge(merge(complete, seRowData, by = "row_id"),
                      seColData, by = "col_id")
    complete <- complete[ order(complete$col_id, complete$row_id), ]
    complete$factor_id <- 1:nrow(complete)
    data_assigned <-
      merge(data, complete, by = c(cond_entries, cl_entries))
    bm <- BumpyMatrix::splitToBumpyMatrix(data_assigned[, allMetadata$dataCols],
                                    row = data_assigned$name_.x,
                                    column = data_assigned$name_.y)

    if (data_type == "untreated") {
      untreatedConds <-
        .get_untreated_conditions(seRowData)
      return(bm[rownames(bm) %in% untreatedConds,])
    } else if (data_type == "treated") {
      treatedConds <- .get_treated_conditions(seRowData)
      return(bm[rownames(bm) %in% treatedConds,])
    } else if (data_type == "all") {
      return(bm)
    } else {
      stop(sprintf("bad 'data_type': ('%s')", data_type))
    }
  }
