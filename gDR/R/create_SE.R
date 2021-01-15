#' createSE
#' TODO: delete me for create_SE function. 
#' @export
createSE <- function(...) {
  .Deprecated("create_SE", package="gDR")
  create_SE(...)
}


#' create_SE
#'
#' Create SummarizedExperiment object from data.frame(s) with dose-reponse data
#'
#' @param df_data a dataframe with DR data
#' @param data_type string type of data to be returned: with untreated conditions only ('untreated'),
#' with treated conditions only ('treated') or all
#' @param readout a character that indiciated readout value in a df ('ReadoutValue' by default)
#' @param discard_keys a vector of keys that should be discarded
#' @param assay_type string assay type for the returned SummarizedExperiment
#'
#' @return SummarizedExperiment object with dose-reponse data
#'
#' @export
create_SE <-
  function(df_data,
           data_type = c("untreated", "treated", "all"),
           readout = 'ReadoutValue',
           discard_keys = NULL,
           assay_type = c("matrix", "BumpyMatrix")) {
    # Assertions:
    stopifnot(any(inherits(df_data, "data.frame"), inherits(df_data, "DataFrame")))
    checkmate::assert_character(data_type)
    checkmate::assert_string(readout)
    checkmate::assert_character(discard_keys, null.ok = TRUE)
    data_type <- match.arg(data_type)
    assay_type <- match.arg(assay_type)

    # stopifnot(all(names(dfList) %in% .assayNames))
    # #dfList must contain first assay (i.e. df_raw_data)
    # stopifnot(.assayNames[1] %in% names(dfList))

    
    mats <- if (assay_type == "matrix") {
      df_to_assay(df_data, data_type = data_type, discard_keys = discard_keys)
    } else if (assay_type == "BumpyMatrix") {
      df_to_bm_assay(df_data, data_type = data_type, discard_keys = discard_keys)
    } else {
      stop(sprintf("bad assay type: '%s'", assay_type))
    }
    
    allMetadata <- getMetaData(df_data, discard_keys = discard_keys)

    seColData <- allMetadata$colData
    rownames(seColData) <- seColData$name_
    seRowData <- allMetadata$rowData
    rownames(seRowData) <- seRowData$name_

    seColData <-
      seColData[colnames(mats), setdiff(colnames(seColData), c('col_id', 'name_'))]
    seRowData <- seRowData[rownames(mats),
                           setdiff(colnames(seRowData), c('row_id', 'name_'))]
    matsL <- list(mats)
    names(matsL) <- readout
    
    se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
                                                     colData = seColData,
                                                     rowData = seRowData)

  }
