#' @import S4Vectors

#### AUXILIARY FUNCTIONS ####
.drugNameRegex <- "^DrugName$|^DrugName_[[:digit:]]+$"
.untreateDrugNameRegex <- "^untreated$|^vehicle$"
.assayNames <-
  c("df_raw_data",
    "df_normalized",
    "df_averaged",
    "df_metrics")

#' .get_untreated_conditions
#'
#' Get untreated conditions
#'
#' @param drug_data tibble or data.frame with treatment information
#'
#' @return character vector with untreated conditions
#'
.get_untreated_conditions <-
  function(drug_data) {
    as.data.frame(drug_data) %>%
      dplyr::filter(grepl(.untreateDrugNameRegex, DrugName)) %>%
      dplyr::pull("name_")
  }

#' .get_treated_conditions
#'
#' Get treated conditions
#'
#' @param drug_data tibble or data.frame with treatment information
#'
#' @return character vector with treated conditions
#'
.get_treated_conditions <-
  function(drug_data) {
    as.data.frame(drug_data) %>%
      dplyr::filter(!grepl(.untreateDrugNameRegex, DrugName)) %>%
      dplyr::pull("name_")
  }


#' aapply
#'
#' works like sapply but on each nested dataframe of the assay of an SE
#'
#' @param function: function to apply on the nested dataframes
#' @param assay_type: integer or name of the assay on which to apply the function
#'
#' @return the same SE object with updated nested dataframe
#'
aapply <-
function(SE, fx, assay_type = 1) {
    SummarizedExperiment::assay(SE, assay_type) = matrix(sapply(SummarizedExperiment::assay(SE, assay_type), fx), nrow = nrow(SE), ncol = ncol(SE))
    return(SE)
}



#' getMetaData
#'
#' Get metadata of given project
#'
#' @param data tibble or data.frame with drug-response data
#' @param cell_id string name of the column with cell line names
#'
#' @return list with two DataFrames ('colData' and 'rowData') and character vector ('dataCols')
#'
#' @export
getMetaData <- function(data,
                        cell_id = gDR::get_identifier("cellline"),
                        discard_keys = NULL) {
  data <- as(data, "DataFrame")

  # get the metadata variables
  metavars <-
    setdiff(
      colnames(data),
      c(
        gDR::get_header("raw_data"),
        gDR::get_header("normalized_results"),
        gDR::get_header("averaged_results"),
        gDR::get_header("metrics_results"),
        gDR::get_identifier("WellPosition"),
        "Barcode",
        "Template",
        # not sure how to handle these ones ....    < --------
        "Concentration"
      )
    ) # remove as it will be the third dimension
  # find all unique conditions
  conditions <- unique(data[, metavars])
  # get the metadata not directly related to cells
  nocell_metavars <- setdiff(metavars,
                             c(get_identifier("cellline"), gDR::get_header("add_clid")))
  constant_metavars <-
    setdiff(
      nocell_metavars[sapply(nocell_metavars,
                             function(x)
                               nrow(unique(conditions[, x, drop = FALSE]))) == 1],
      # protect cell line and drug name
      c(
        gDR::get_header("add_clid"),
        gDR::get_identifier("drug"),
        gDR::get_identifier("drugname")
      )
    )
  unique_metavars <- c(intersect(
    c(
      gDR::get_identifier("cellline"),
      gDR::get_header("add_clid"),
      gDR::get_identifier("drug"),
      gDR::get_identifier("drugname")
    ),
    metavars
  ),
  nocell_metavars[sapply(nocell_metavars, function(x)
    nrow(unique(conditions[, x, drop = FALSE]))) > 1])

  # find the cell lines and related data (for the columns in the SE)
  cl_entries <- cell_id
  for (j in setdiff(unique_metavars, cell_id)) {
    if (nrow(unique(conditions[, c(cell_id, j)])) ==
        nrow(unique(conditions[, cell_id, drop = FALSE]))) {
      cl_entries = c(cl_entries, j)
    }
  }
  # --> not very robust --> need testing
  cl_entries <- setdiff(cl_entries,
                        c(
                          gDR::get_identifier("drug"),
                          paste0(gDR::get_identifier("drug"), "_", 2:10),
                          gDR::get_identifier("drugname"),
                          paste0(gDR::get_identifier("drugname"), "_", 2:10)
                        ))

  # temporary removing extra column to avoid bug
  cl_entries <- setdiff(cl_entries, "ReferenceDivisionTime")
  
  #colData
  colData <- unique(conditions[, cl_entries, drop = FALSE])
  colData$col_id <- 1:nrow(colData)
  colData$name_ <-
    apply(colData, 1, function(x)
      paste(x, collapse = "_"))
  
  # get all other metadata for the rows
  cond_entries <- setdiff(unique_metavars, cl_entries)
  # temporary removing extra column to avoid bug
  cond_entries <- setdiff(cond_entries, c('ReferenceDivisionTime', discard_keys))
  rowData <- unique(conditions[, cond_entries, drop = FALSE])
  rowData$row_id <- 1:nrow(rowData)
  rowData$name_ <-
    apply(rowData, 1, function(x)
      paste(x, collapse = "_"))
  
  # get the remaining columns as data
  dataCols <- setdiff(colnames(data), setdiff(metavars, discard_keys))

  # constant metadata (useful for annotation of the experiment)
  csteData = unique(conditions[,constant_metavars,drop=F])

  return(list(
    colData = colData,
    rowData = rowData,
    dataCols = dataCols,
    csteData = csteData
  ))
}

#' df_to_assay
#'
#' Convert data.frame with dose-reponse data to the assay.
#'
#' The 'assay' object is simply a matrix with rownames being the treatment ids,
#' colnames being the ids of the cell lines and values being the DataFrames with
#' dose-response data for given cell lines under given conditions.
#'
#' @param data tibble or data.frame with drug-response data
#' @param treatment_id string id of the column with the treatment conditions data (in 'data' object)
#' @param data_type string type of data to be returned: all, for untreated conditions only ('untreated')
#' or for treated conditions only ('treated')
#'
#' @return matrix
#'
#' @export
df_to_assay <-
  function(data,
           data_type = c("all", "treated", "untreated"),
           discard_keys = NULL) {
    data <- as(data, "DataFrame")

    ####

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
    complete = complete[ order(complete$col_id, complete$row_id), ]
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
      return(final.set[rownames(final.set) %in% untreatedConds, ])
    } else if (data_type == "treated") {
      treatedConds <- .get_treated_conditions(seRowData)
      return(final.set[rownames(final.set) %in% treatedConds, ])
    } else if (data_type == "all") {
      return(final.set)
    } else {
      stop(sprintf("bad 'data_type': ('%s')", data_type))
    }
  }

#' createSE
#'
#' Create SummarizedExperiment object from data.frame(s) with dose-reponse data
#'
#' @param dfList list with data.frame(s) with dose-reponse data
#' @param data_type string type of data to be returned: with untreated conditions only ('untreated'),
#' with treated conditions only ('treated') or all
#'
#' @return SummarizedExperiment object with dose-reponse data
#'
#' @export
createSE <-
  function(df_data,
           data_type = c("untreated", "treated", "all"),
           readout = 'ReadoutValue', discard_keys = NULL) {
    data_type <- match.arg(data_type)

    # stopifnot(all(names(dfList) %in% .assayNames))
    # #dfList must contain first assay (i.e. df_raw_data)
    # stopifnot(.assayNames[1] %in% names(dfList))

    mats <- df_to_assay(df_data, data_type = data_type, discard_keys = discard_keys)

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
    stopifnot("MultiAssayExperiment" %in% class(mae))
    stopifnot(assay_name %in% .assayNames)
    stopifnot("matrix" %in% class(assay))
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

#' assay_to_df
#'
#' Convert SE assay to data.frame
#'
#' @param se  SummarizedExperiment object with dose-response data
#' @param assay_name string name of the assay
#'
#' @return data.frame with dose-reponse data
#'
#' @export
assay_to_df <- function(se, assay_name) {
  stopifnot(any("SummarizedExperiment" %in% class(se)))
  #checkmate::assertString(assay_name)
  
  # define data.frame with data from rowData/colData
  ids <- expand.grid(rownames(SummarizedExperiment::rowData(se)), rownames(SummarizedExperiment::colData(se)))
  colnames(ids) <- c("rId", "cId")
  ids[] <- lapply(ids, as.character)
  rData <- data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)
  rData$rId <- rownames(rData)
  cData <- data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  cData$cId <- rownames(cData)
  annotTbl <-
    dplyr::left_join(ids, rData, by = "rId")
  annotTbl <-
    dplyr::left_join(annotTbl, cData, by = "cId")
  
  #merge assay data with data from colData/rowData
  asL <- lapply(1:nrow(SummarizedExperiment::colData(se)), function(x) {
    myL <- SummarizedExperiment::assay(se, assay_name)[, x]
    myV <- vapply(myL, nrow, integer(1))
    rCol <- rep(names(myV), as.numeric(myV))
    
    df <- data.frame(do.call(rbind, myL))
    df$rId <- rCol
    df$cId <- rownames(SummarizedExperiment::colData(se))[x]
    full.df <- left_join(df, annotTbl, by = c("rId", "cId"))
  })
  asDf <- data.frame(do.call(rbind, asL))
  if (assay_name == "Metrics") {
    asDf$dr_metric <- c("IC", "GR")
  }
  asDf
}  
