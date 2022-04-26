#' add_CellLine_annotation
#'
#' add cellline annotation to a data.frame with metadata
#'
#' @param df_metadata a data.frame with metadata
#' @param DB_cellid_header a string with colnames with cell line identifier in the annotation file
#' @param DB_cell_annotate a character vector with mandatory colnames used in the annotation file
#' @param fname a string with file name with annotation
#' @param fill a string indicating how unknown cell lines should be filled in the DB
#' @details   The logic of adding celline annotation for df_metadata based on the annotation file stored
#' in gDRtestData.
#' Other fields are set as "unknown".
#' This approach will be corrected once we will implement final solution for adding cell lines.
#' @return a data.frame with metadata with annotated cell lines
#' @export
#'
add_CellLine_annotation <- function(df_metadata,
                                    DB_cellid_header = "cell_line_identifier",
                                    DB_cell_annotate = c("cell_line_name", "primary_tissue", "doubling_time",
                                                         "parental_identifier", "subtype"),
                                    fname = "cell_lines.csv",
                                    fill = "unknown") {
  
  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))
  checkmate::assert_string(fill, null.ok = TRUE)
  data.table::setDT(df_metadata)

  cellline <- gDRutils::get_env_identifiers("cellline")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  add_clid <- gDRutils::get_header("add_clid")
  
  # Read local cell_lines annotations
  annotationPackage <- ifelse(requireNamespace("gDRinternalData", quietly = TRUE),
                              "gDRinternalData", "gDRtestData")
  CLs_info <- read.csv(system.file("data", fname, package = annotationPackage))
  CLs_info <- CLs_info[, c(DB_cellid_header, DB_cell_annotate)]
  
  if (nrow(CLs_info) == 0) return(df_metadata)
  
  validatedCLs <- unique(df_metadata[[cellline]]) %in% CLs_info[[DB_cellid_header]]
  missingTblCellLines <- NULL
  if (!is.null(fill) && !all(validatedCLs)) {
    unValidatedCellLine <- unique(df_metadata[[cellline]])[!validatedCLs]
    missingTblCellLines <- data.table::data.table(parental_identifier = unValidatedCellLine,
                                                  cell_line_name = unValidatedCellLine,
                                                  cell_line_identifier = unValidatedCellLine,
                                                  doubling_time = NA,
                                                  primary_tissue = fill,
                                                  subtype = fill)
    
  }
  
  if (any(!df_metadata$clid %in% CLs_info[[DB_cellid_header]]) && !is.null(missingTblCellLines)) {
    CLs_info <- rbind(CLs_info, missingTblCellLines)
  }
  CLs_info[is.na(CLs_info)] <- fill
  CLs_info[, "doubling_time"] <- as.numeric(CLs_info[, "doubling_time"])
  
  colnames(CLs_info) <- unlist(c(cellline, add_clid, tail(DB_cell_annotate, 2)))
  CLIDs <- unique(df_metadata[[cellline]])
  bad_CL <- CLs_info[cellline %in% CLIDs][[cellline]]
  if (any(bad_CL)) {
    futile.logger::flog.warn("Cell line ID %s not found in cell line database",
                             paste(CLIDs[bad_CL], collapse = " ; "))
    temp_CLIDs <- as.data.frame(matrix(ncol = length(c(cellline,
                                                       add_clid)),
                                       nrow = length(CLIDs[bad_CL])))
    colnames(temp_CLIDs) <- c(cellline,
                              add_clid)
    
    temp_CLIDs[, cellline] <- temp_CLIDs[, cellline_name] <- CLIDs[bad_CL]
    CLs_info <- rbind(CLs_info, temp_CLIDs, fill = TRUE)
  }
  
  futile.logger::flog.info("Merge with Cell line info")
  nrows_df <- nrow(df_metadata)
  df_metadata <- base::merge(df_metadata, CLs_info, by = cellline, all.x = TRUE)
  stopifnot(nrows_df == nrow(df_metadata))
  df_metadata
}


#' add_Drug_annotation
#'
#' add drug annotation to a data.frame with metadata
#'
#' @param df_metadata a data.frame with metadata
#' @param fname a string with file name with annotation
#' @param fill a string indicating how unknown cell lines should be filled in the DB
#' @details The logic of adding drug annotation for df_metadata based on the annotation file stored
#' in gDRtestData.
#' @return a data.frame with metadata with annotated drugs
#' @export
#'
add_Drug_annotation <- function(df_metadata,
                                fname = "drugs.csv",
                                fill = "unknown") {
  
  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))
  checkmate::assert_string(fill, null.ok = TRUE)
  data.table::setDT(df_metadata)
  nrows_df <- nrow(df_metadata)
  
  drug <- gDRutils::get_env_identifiers("drug")
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drugsTreated <- unique(df_metadata[[drug]])
  
  # Read local drugs annotations
  annotationPackage <- ifelse(requireNamespace("gDRinternalData", quietly = TRUE),
                              "gDRinternalData", "gDRtestData")
  Drug_info <- read.csv(system.file("data", fname, package = annotationPackage),
                        header = TRUE)
  Drug_info <- Drug_info[, c("gnumber", "drug_name", "drug_moa")]
  data.table::setnames(Drug_info, c("drug", "drug_name", "drug_moa"))
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
    df_metadata[, drug_name] <- df_metadata[, drug]
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
  DrIDs <- unique(unlist(df_metadata[, grep(drug, colnames(df_metadata))]))
  if (any(!remove_drug_batch(drugsTreated) %in% Drug_info$drug) && !is.null(missingTblDrugs)) {
    Drug_info <- rbind(Drug_info, data.table::setnames(
      missingTblDrugs[!(remove_drug_batch(missingTblDrugs$drug) %in% Drug_info$drug), ],
      names(Drug_info)))
  }
  DrIDs[is.na(DrIDs)] <- fill
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
  colnames(Drug_info) <- c("drug", drug_name, drug_moa)
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
    colnames(Drug_info_)[colnames(Drug_info_)
                         %in% c(drug_name, drug_moa)] <- paste0(c(drug_name, drug_moa),
                                                               gsub(".*(_\\d)", "\\1",
                                                                    colnames(df_metadata)[i]))
    futile.logger::flog.info("Merge with Drug_info for %s", colnames(df_metadata)[[i]])
    df_metadata[[paste0(colnames(df_metadata)[[i]], "_temp")]] <- remove_drug_batch(df_metadata[[i]])
    df_metadata <- base::merge(df_metadata, Drug_info_, by.x =
                                 remove_drug_batch(paste0(colnames(df_metadata)[[i]],
                                                          "_temp")), by.y = "drug", all.x = TRUE)
    df_metadata <- df_metadata[, .SD, .SDcols = !endsWith(names(df_metadata), "temp")]
  }
  stopifnot(nrows_df == nrow(df_metadata))
  df_metadata
}

#' Remove batch from Gnumber
#'
#' @return Gnumber without a batch
#' @export
remove_drug_batch <- function(drug) {
  gsub("\\.[0-9]+.*", "", drug)
}
