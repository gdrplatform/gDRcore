#' df_to_assay
#'
#' Convert data.frame with dose-reponse data to the assay.
#'
#' The 'assay' object is simply a matrix with rownames being the treatment ids,
#' colnames being the ids of the cell lines and values being the DataFrames with
#' dose-response data for given cell lines under given conditions.
#'
#' @param data tibble or data.frame with drug-response data
#' @param data_type string type of data to be returned: all, for untreated conditions only ('untreated')
#' or for treated conditions only ('treated')
#' @param discard_keys a vector of keys that should be discarded
#'
#' @return matrix
#'
#' @export
df_to_assay <-
  function(data,
           data_type = c("all", "treated", "untreated"),
           discard_keys = NULL) {
    # Assertions:
    stopifnot(any(inherits(data, "data.frame"), checkmate::test_character(data), inherits(data, "DataFrame")))
    checkmate::assert_character(data_type)
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
    
    by_factor <- lapply(1:nrow(complete), function(x)
      data_assigned[data_assigned$factor_id == x, dataCols])
    names(by_factor) <- 1:nrow(complete)

    stopifnot(nrow(data) == sum(sapply(by_factor, nrow)))
    stopifnot(length(by_factor) == nrow(complete))

    # full.set <- vector("list", nrow(complete))
    # full.set[as.integer(names(by_factor))] <- by_factor
    full.set <- by_factor

    dim(full.set) <- c(nrow(seRowData), nrow(seColData))
    dimnames(full.set) <- list(seRowData$name_, seColData$name_)

    #add NAs for treatments not present in the given assay
    # ---------------
    # removed as it should be added when combining different assays
    #
    # tNotFound <- setdiff(adrug, rownames(full.set))
    # mNotFound <-
    #   matrix(nrow = length(tNotFound), ncol = ncol(full.set))
    # dimnames(mNotFound) <- list(tNotFound, colnames(full.set))
    #
    # final.set <- rbind(mNotFound, full.set)
    #
    # ---------------
    final.set <- full.set

    if (data_type == "untreated") {
      untreatedConds <-
        .get_untreated_conditions(seRowData)
      return(final.set[rownames(final.set) %in% untreatedConds, , drop = F])
    } else if (data_type == "treated") {
      treatedConds <- .get_treated_conditions(seRowData)
      return(final.set[rownames(final.set) %in% treatedConds, , drop = F])
    } else if (data_type == "all") {
      return(final.set)
    } else {
      stop(sprintf("bad 'data_type': ('%s')", data_type))
    }
  }
