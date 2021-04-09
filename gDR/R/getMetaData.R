#' getMetaData
#'
#' Get metadata of given project
#'
#' @param data data.frame with drug-response data
#' @param discard_keys character vector of keys to discard from the data columns
#'
#' @return named list containing different elements of a \linkS4class{SummarizedExperiment};
#' see details.
#' 
#' @details
#' Named list containing the following elements:
#' \itemize{
#'  \item{rowData}{treatment metadata}
#'  \item{colData}{condition metadata}
#'  \item{dataCols}{all data.frame columns corresponding to non-metadata fields}
#'  \item{csteData}{metadata that is constant for all entries of the data.frame}
#' }
#'
#' @export
#'
getMetaData <- function(data, discard_keys = NULL) {
  .Deprecated(msg = "see split_SE_components for similar, but not identical functionality")

  # Assertions.
  stopifnot(any(inherits(data, "data.frame"), inherits(data, "DataFrame")))
  checkmate::assert_character(discard_keys, null.ok = TRUE)

  cell_id <- gDRutils::get_identifier("cellline")

  data <- S4Vectors::DataFrame(data)
  all_data_cols <- colnames(data)

  # Separate out metadata versus data variables.
  metavars <-
    setdiff(
      all_data_cols,
      c(
        gDRutils::get_header("raw_data"),
        gDRutils::get_header("normalized_results"),
        gDRutils::get_header("averaged_results"),
        gDRutils::get_header("metrics_results"),
        gDRutils::get_identifier("well_position"),
        "Barcode",
        "Template",
        # not sure how to handle these ones ....    < --------
        "Concentration"
      )
    ) # remove as it will be the third dimension

  conditions <- unique(data[, metavars])

  # Remove cell-related metadata.
  nocell_metavars <- setdiff(metavars,
                             c(gDRutils::get_identifier("cellline"), gDRutils::get_header("add_clid")))
  singleton_cols <- vapply(nocell_metavars,
			   function(x) {nrow(unique(conditions[, x, drop = FALSE])) == 1L},
			   logical(1))

  # Remove drug metadata and duration.
  constant_metavars <- setdiff(nocell_metavars[singleton_cols],
			       c(gDRutils::get_identifier("drug"),
				 gDRutils::get_identifier("drugname"),
				 gDRutils::get_identifier("duration")
			        ))

  unique_metavars <- c(intersect(c(gDRutils::get_identifier("cellline"),
				   gDRutils::get_header("add_clid"),
				   gDRutils::get_identifier("drug"),
				   gDRutils::get_identifier("drugname"),
				   gDRutils::get_identifier("duration")),
				   metavars),
                        nocell_metavars[!singleton_cols])

  cl_entries <- cell_id
  for (j in setdiff(unique_metavars, cell_id)) {
    if (nrow(unique(conditions[, c(cell_id, j)])) ==
        nrow(unique(conditions[, cell_id, drop = FALSE]))) {
      cl_entries <- c(cl_entries, j)
    }
  }

  pattern <- sprintf("^%s*|^%s*|^%s$", 
    gDRutils::get_identifier("drug"), 
    gDRutils::get_identifier("drugname"),
    gDRutils::get_identifier("duration"))
  cl_entries <- cl_entries[!grepl(pattern, cl_entries)]

  ## colData
  colData <- unique(conditions[, cl_entries, drop = FALSE])
  colData <- colData[sort(colnames(colData))] # TODO: can delete me later. Sorting for simple comparison. 
  colData$col_id <- seq_len(nrow(colData))

  colData$name_ <-
    apply(colData, 1, function(x)
      paste(x, collapse = "_"))

  ## rowData
  cond_entries <- setdiff(unique_metavars, c(cl_entries, discard_keys))
  rowData <- unique(conditions[, cond_entries, drop = FALSE])
  rowData <- rowData[sort(colnames(rowData))] # TODO: can delete me later. Sorting for simple comparison. 
  rowData$row_id <- seq_len(nrow(rowData))
	  
  rowData$name_ <-
    apply(rowData, 1, function(x)
      paste(x, collapse = "_"))

  ## dataCols
  dataCols <- setdiff(all_data_cols, setdiff(metavars, discard_keys))

  ## constant metadata (useful for annotation of the experiment)
  csteData <- unique(conditions[, constant_metavars, drop = FALSE])

  return(list(
    colData = colData,
    rowData = rowData,
    dataCols = dataCols,
    csteData = csteData
  ))
}
