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


#' addAssayToMAE
#'
#' Add assay to one of MAEs experiments (i.e. SEs) with dose-reponse data
#'
#' @param mae  MultiAssayExperiment object with dose-reponse
#' @param assay matrix with dose-response data
#' @param exp_name string name of the MAE experiment (i.e. SEs) to which add the assay
#' @param assay_name string name of the assay to be used in SE
#' @param update_assay logical allow for assay update if the assay with 'assay_name' currently exists in given 'exp_name'
#'
#' @return MultiAssayExperiment object with dose-reponse data
#'
#' @export
addAssayToMAE <-
  function(mae,
           assay,
           assay_name,
           exp_name = c("treated", "untreated"),
           update_assay = FALSE) {
    # Assertions:
    stopifnot("MultiAssayExperiment" %in% class(mae))
    stopifnot(assay_name %in% .assayNames)
    stopifnot("matrix" %in% class(assay))
    checkmate::assert_logical(update_assay)

    exp_name <- match.arg(exp_name)
    #mae must contain SE with at least first assay (i.e. df_raw_data)
    stopifnot(.assayNames[1] %in% SummarizedExperiment::assayNames(mae[[exp_name]]))

    if (assay_name %in% SummarizedExperiment::assayNames(mae[[exp_name]]) &&
        update_assay == FALSE) {
      futile.logger::flog.error(
        "The assay '%s' can't be added to experiment '%s' as it currently exists.
        Please set 'update_assay' flag to TRUE to be able to update the assay instead of adding it",
        assay_name,
        exp_name
      )
      stop()
    }

    if (!identical(dim(SummarizedExperiment::assay(mae[[exp_name]])), dim(assay))) {
      futile.logger::flog.error(
        "The assay '%s' can't be added to experiment '%s' as it has different dimensions ('%s') than the assays present in the experiment ('%s').",
        assay_name,
        exp_name,
        paste(dim(assay), collapse = "x"),
        paste(dim(SummarizedExperiment::assay(mae[[exp_name]])), collapse = "x")
      )
      stop()
    }

    SummarizedExperiment::assay(mae[[exp_name]], assay_name) <- assay
  }

