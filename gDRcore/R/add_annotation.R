#' add_CellLine_annotation
#'
#' add cellline annotation to a data.frame with metadata
#'
#' @param df_metadata a data.frame with metadata
#' @param fill a string indicating how unknown cell lines should be filled in the DB
#' @details   The logic of adding celline annotation for df_metadata is based on the function
#' get_cell_lines from the gDRwrapper
#' we added additional parameter 'fill' that allows to fill the DB with clid info for these cell lines
#' that are not present in the DB.
#' Other fields are set as "UNKNOWN". If the fill is set as FALSE we add unknown cell lines
#' only to the tibble.
#' This approach will be corrected once we will implement final solution for adding cell lines.
#' @return a data.frame with metadata with annotated cell lines
#' @export
#'
add_CellLine_annotation <- function(df_metadata,
                                    fill = NULL) {

  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))
  checkmate::assert_string(fill, null.ok = TRUE)

  DB_cellid_header <- "cell_line_identifier"
  DB_cell_annotate <- c("cell_line_name", "primary_tissue", "doubling_time",
                        "parental_identifier", "subtype")
  cellline <- gDRutils::get_identifier("cellline")
  add_clid <- gDRutils::get_header("add_clid")

  # Read local cell_lines annotations
  CLs_info <- data.table::fread(system.file("data/cell_lines.csv", package = "gDRcore"))
  CLs_info <- CLs_info[, c(DB_cellid_header, DB_cell_annotate), with = FALSE]

  if (nrow(CLs_info) == 0) return(df_metadata)

  validatedCLs <- all(unique(df_metadata[[cellline]]) %in% CLs_info[[DB_cellid_header]])
  missingTblCellLines <- NULL
  if (!is.null(fill) && !validatedCLs) {
    missingTblCellLines <- data.table::data.table(parental_identifier = fill,
                                                  cell_line_name = fill,
                                                  cell_line_identifier = unlist(unique(df_metadata[,
                                                                                                   cellline, with = FALSE])),
                                                  doubling_time = NA,
                                                  primary_tissue = fill,
                                                  subtype = fill)

  }

  if (any(!df_metadata$clid %in% CLs_info[[DB_cellid_header]]) && !is.null(missingTblCellLines)) {
    CLs_info <- rbind(CLs_info, missingTblCellLines)
  }

  colnames(CLs_info)[1:4] <- c(cellline, add_clid)
  CLIDs <- unique(df_metadata[[cellline]])
  bad_CL <- CLs_info[cellline %in% CLIDs][[cellline]]
  if (any(bad_CL)) {
    futile.logger::flog.warn("Cell line ID %s not found in cell line database",
                             paste(CLIDs[bad_CL], collapse = " ; "))
    temp_CLIDs <- data.frame(CLIDs[bad_CL], CLIDs[bad_CL])
    temp_CLIDs[, 1 + (2:length(add_clid))] <- NA
    colnames(temp_CLIDs) <- c(cellline,
                              add_clid)
    CLs_info <- rbind(CLs_info, temp_CLIDs)
  }

  futile.logger::flog.info("Merge with Cell line info")
  nrows_df <- nrow(df_metadata)
  df_metadata <- base::merge(df_metadata, CLs_info, by = cellline, all.x = TRUE)
  stopifnot(nrows_df == nrow(df_metadata))
  return(df_metadata)
}


#' add_Drug_annotation
#'
#' add drug annotation to a data.frame with metadata
#'
#' @param df_metadata a data.frame with metadata
#' @param fill a string indicating how unknown cell lines should be filled in the DB
#'
#'
#' @return a data.frame with metadata with annotated drugs
#' @export
#'
add_Drug_annotation <- function(df_metadata,
                                fill = NULL) {

  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))
  checkmate::assert_string(fill, null.ok = TRUE)

  nrows_df <- nrow(df_metadata)

  DB_drug_identifier <- c("gnumber", "drug_name", "drug_moa")
  drug <- gDRutils::get_identifier("drug")
  untreated_tag <- gDRutils::get_identifier("untreated_tag")
  drugname <- gDRutils::get_identifier("drugname")
  drugsTreated <- unique(df_metadata[[drug]])

  # Read local drugs annotations

  Drug_info <- data.table::fread(system.file("data/drugs.csv", package = "gDRcore"))
  colnames(Drug_info)[1:2] <- c("drug", "drug_name")
  drugsTreated <- drugsTreated[!drugsTreated %in% untreated_tag]
  validatedDrugs <- drugsTreated %in% Drug_info[["drug"]]
  #### function should be parallelized
  missingTblDrugs <- NULL
  if (!is.null(fill) && any(!validatedDrugs)) {
    missingTblDrugs <- tibble::tibble(drug = remove_drug_batch(drugsTreated[!validatedDrugs]),
                                      drug_name = drugsTreated[!validatedDrugs],
                                      drug_moa = fill
    )
  }


  if (nrow(Drug_info) == 0) {
    df_metadata[, drugname] <- df_metadata[, drug]
    return(df_metadata)
  }

  # -----------------------
  Drug_info$drug <- remove_drug_batch(Drug_info$drug)
  Drug_info <-
    rbind(data.frame(
      drug = untreated_tag,
      drug_name = untreated_tag,
      drug_moa = untreated_tag
    ),
    Drug_info)
  Drug_info <- Drug_info[!duplicated(Drug_info[["drug"]]), ]
  DrIDs <- unique(unlist(df_metadata[, grep(drug, colnames(df_metadata)), with = FALSE]))
  if (any(!remove_drug_batch(drugsTreated) %in% Drug_info$drug) && !is.null(missingTblDrugs)) {
    Drug_info <- rbind(Drug_info, data.table::setnames(
      missingTblDrugs[!(remove_drug_batch(missingTblDrugs$drug) %in% Drug_info$drug), ],
      names(Drug_info)))
  }
  bad_DrID <- !(remove_drug_batch(DrIDs) %in% Drug_info$drug) & !is.na(DrIDs)
  if (any(bad_DrID)) {
    # G number, but not registered
    ok_DrID <- attr(regexpr("^G\\d*", DrIDs), "match.length") == 9
    if (any(ok_DrID)) {
      futile.logger::flog.warn("cleanup_metadata: Drug %s  not found in gCSI database;
                               use G# as DrugName",
                               paste(DrIDs[ok_DrID & bad_DrID], collapse = " ; "))
      default_drug_moa <- ifelse(is.null(fill), "unknown", fill)
      Drug_info <-
        rbind(Drug_info, data.frame(drug = remove_drug_batch(DrIDs[ok_DrID & bad_DrID]),
                                    drug_name = remove_drug_batch(DrIDs[ok_DrID & bad_DrID]),
                                    drug_moa = default_drug_moa))
    } else {
      futile.logger::flog.error("Drug %s not in the correct format for database",
                                paste(DrIDs[!ok_DrID], collapse = " ; "))
    }
  }
  colnames(Drug_info)[2] <- drugname
  futile.logger::flog.info("Merge with Drug_info for Drug 1")
  df_metadata[[paste0(drug, "_temp")]] <-
    remove_drug_batch(df_metadata[[drug]])
  df_metadata <- base::merge(df_metadata, Drug_info, by.x =
                               paste0(drug, "_temp"), by.y = "drug", all.x = TRUE)
  df_metadata <- df_metadata[, .SD, .SDcols = !endsWith(names(df_metadata), "temp")]
  # add info for columns Gnumber_*
  for (i in grep(paste0(drug, "_\\d"), colnames(df_metadata))) {
    df_metadata[is.na(df_metadata[[i]]), i] <- untreated_tag[1] # set missing values to Untreated
    Drug_info_ <- Drug_info
    colnames(Drug_info_)[2:length(colnames(Drug_info_))] <- paste0(colnames(Drug_info_)[c(2, 3)],
                                                                   substr(colnames(df_metadata)[i], 8, 12))
    futile.logger::flog.info("Merge with Drug_info for %s", colnames(df_metadata)[[i]])
    df_metadata[[paste0(colnames(df_metadata)[[i]], "_temp")]] <- remove_drug_batch(df_metadata[[i]])
    df_metadata <- base::merge(df_metadata, Drug_info_, by.x =
                                 remove_drug_batch(paste0(colnames(df_metadata)[[i]],
                                                          "_temp")), by.y = "drug", all.x = TRUE)
    df_metadata <- df_metadata[, .SD, .SDcols = !endsWith(names(df_metadata), "temp")]
  }
  stopifnot(nrows_df == nrow(df_metadata))

  return(df_metadata)
}