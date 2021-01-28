#' @import S4Vectors

#### AUXILIARY FUNCTIONS ####
.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$", 
                          gDRutils::get_identifier("drugname"), 
                          gDRutils::get_identifier("drugname"))

.untreated_tag_patterns <- vapply(gDRutils::get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
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
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with untreated conditions
#'
#' @export
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
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with treated conditions
#'
#' @export
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
