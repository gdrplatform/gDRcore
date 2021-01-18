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
#'
getMetaData <- function(data,
                        cell_id = gDRutils::get_identifier("cellline"),
                        discard_keys = NULL) {
  # Assertions:
  stopifnot(any(inherits(data, "data.frame"), inherits(data, "DataFrame")))
  checkmate::assert_character(cell_id)
  checkmate::assert_character(discard_keys, null.ok = TRUE)

  data <- methods::as(data, "DataFrame")

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

  unique_metavars <- c(intersect(c(gDRutils::get_identifier("cellline"),
				   gDRutils::get_header("add_clid"),
				   gDRutils::get_identifier("drug"),
				   gDRutils::get_identifier("drugname"),
				   gDRutils::get_identifier("duration")
                                  ),
				  metavars),
                       nocell_metavars[vapply(nocell_metavars, function(x) {length(unique(conditions[[x]])) > 1}, logical(1))])

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
    apply(colData[, grep('parental_identifier|subtype', colnames(colData), invert = TRUE)], 1, function(x)
      paste(x, collapse = "_"))

  # get all other metadata for the rows
  cond_entries <- setdiff(unique_metavars, cl_entries)
  # temporary removing extra column to avoid bug
  cond_entries <- setdiff(cond_entries, c('ReferenceDivisionTime', discard_keys))
  rowData <- unique(conditions[, cond_entries, drop = FALSE])
  rowData$row_id <- 1:nrow(rowData)
  rowData$name_ <-
    apply(rowData[, grep('drug_moa', colnames(rowData), invert = TRUE)], 1, function(x)
      paste(x, collapse = "_"))

  # get the remaining columns as data
  dataCols <- setdiff(colnames(data), setdiff(metavars, discard_keys))

  # constant metadata (useful for annotation of the experiment)
  csteData <- unique(conditions[, constant_metavars, drop = FALSE])

  return(list(
    colData = colData,
    rowData = rowData,
    dataCols = dataCols,
    csteData = csteData
  ))
}
