#' @import S4Vectors

#### AUXILIARY FUNCTIONS ####
.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$", 
                          gDRutils::get_identifier("drugname"), 
                          gDRutils::get_identifier("drugname"))

.untreated_tag_patterns <- vapply(gDRutils::get_identifier("untreated_tag"), function(x) {sprintf("^%s$", x)}, character(1))
.untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse="|")

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
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    drugnames <- tolower(as.data.frame(drug_data)[, gDRutils::get_identifier("drugname")])
    drug_data[grepl(.untreatedDrugNameRegex, drugnames), "name_"]
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
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    drugnames <- tolower(as.data.frame(drug_data)[, gDRutils::get_identifier("drugname")])
    drug_data[!grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }


#' aapply
#'
#' Works like sapply but on each nested dataframe of the assay of an SE
#'
#' @param SE a SummarizedExperiment object
#' @param fx any function
#' @param assay_type a name of an assay
#'
#' @return the same SE object with updated nested dataframe
#' @export
aapply <-
  function(SE, fx, assay_type = 1) {
    # Assertions:
    checkmate::assert_class(SE, "SummarizedExperiment")
    checkmate::assert_function(fx)
    checkmate::assert_scalar(assay_type)

    SummarizedExperiment::assay(SE, assay_type, withDimnames=FALSE) = matrix(sapply(SummarizedExperiment::assay(SE, assay_type), fx), nrow = nrow(SE), ncol = ncol(SE))
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
                        cell_id = gDRutils::get_identifier("cellline"),
                        discard_keys = NULL) {
  # Assertions:
  stopifnot(any(inherits(data, "data.frame"), inherits(data, "DataFrame")))
  checkmate::assert_character(cell_id)
  checkmate::assert_character(discard_keys, null.ok = TRUE)

  data <- as(data, "DataFrame")

  # get the metadata variables
  metavars <-
    setdiff(
      colnames(data),
      c(
        gDRutils::get_header("raw_data"),
        gDRutils::get_header("normalized_results"),
        gDRutils::get_header("averaged_results"),
        gDRutils::get_header("metrics_results"),
        gDRutils::get_identifier("WellPosition"),
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
                             c(gDRutils::get_identifier("cellline"), gDRutils::get_header("add_clid")))
  constant_metavars <-
    setdiff(
      nocell_metavars[sapply(nocell_metavars,
                             function(x)
                               nrow(unique(conditions[, x, drop = FALSE]))) == 1],
      # protect cell line and drug name and duration
      c(
        gDRutils::get_header("add_clid"),
        gDRutils::get_identifier("drug"),
        gDRutils::get_identifier("drugname"),
        gDRutils::get_identifier("duration")
      )
    )
  unique_metavars <- c(intersect(
    c(
      gDRutils::get_identifier("cellline"),
      gDRutils::get_header("add_clid"),
      gDRutils::get_identifier("drug"),
      gDRutils::get_identifier("drugname"),
      gDRutils::get_identifier("duration")
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
                          gDRutils::get_identifier("drug"),
                          paste0(gDRutils::get_identifier("drug"), "_", 2:10),
                          gDRutils::get_identifier("drugname"),
                          paste0(gDRutils::get_identifier("drugname"), "_", 2:10),
                          gDRutils::get_identifier("duration")
                        ))

  # # temporary removing extra column to avoid bug
  # cl_entries <- setdiff(cl_entries, "ReferenceDivisionTime")

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
