#### AUXILIARY FUNCTIONS ####
.drugNameRegex <- "^DrugName$|^DrugName_[[:digit:]]+$"
.untreateDrugNameRegex <- "^untreated$|^vehicle$"

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
                        cell_id = "CellLineName") {
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
  cond_entries <- setdiff(cond_entries, "ReferenceDivisionTime")
  rowData <- unique(conditions[, cond_entries, drop = FALSE])
  rowData$row_id <- 1:nrow(rowData)
  rowData$name_ <-
    apply(rowData, 1, function(x)
      paste(x, collapse = "_"))
  
  # get the remaining columns as data
  dataCols <- setdiff(colnames(data), metavars)
  
  return(list(
    colData = colData,
    rowData = rowData,
    dataCols = dataCols
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
           treatment_id = "row_id",
           data_type = c("all", "treated", "untreated")) {
    data <- as(data, "DataFrame")
    
    ####
    
    allMetadata <- getMetaData(data)
    
    seColData <- allMetadata$colData
    cl_entries <- setdiff(colnames(seColData), c("col_id", "name_"))
    seRowData <- allMetadata$rowData
    cond_entries <-
      setdiff(colnames(seRowData), c("row_id", "name_"))
    dataCols <- allMetadata$dataCols
    
    complete <-
      DataFrame(
        expand.grid(
          row_id = seRowData$row_id,
          col_id = seColData$col_id,
          stringsAsFactors = FALSE
        )
      )
    complete <- merge(merge(complete, seRowData, by = "row_id"),
                      seColData, by = "col_id")
    complete$factor_id <- 1:nrow(complete)
    data_assigned <-
      merge(data, complete, by = c(cond_entries, cl_entries))
    
    by_factor <- lapply(1:nrow(complete), function(x)
      data_assigned[data_assigned$factor_id == x, dataCols])
    names(by_factor) <- 1:nrow(complete)
    
    stopifnot(nrow(data) == sum(sapply(by_factor, nrow)))
    
    full.set <- vector("list", nrow(complete))
    full.set[as.integer(names(by_factor))] <- by_factor
    
    
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