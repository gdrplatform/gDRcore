#' add_CellLine_annotation
#'
#' add cellline annotation to a data.table with metadata
#'
#' @param df_metadata data.table with metadata
#' @param DB_cellid_header string with colnames with cell line identifier
#'                         in the annotation file
#' @param DB_cell_annotate character vector with mandatory colnames 
#'                         used in the annotation file
#' @param fname string with file name with annotation
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @param annotationPackage string indication name of the package containing cellline annotation
#' @details
#' The logic of adding celline annotation for df_metadata based on 
#' the annotation file stored in gDRtestData. Other fields are set as "unknown".
#' This approach will be corrected once we will implement final solution 
#' for adding cell lines.
#' 
#' @examples 
#' 
#' add_CellLine_annotation(
#'   data.table::data.table(
#'     clid = "123", 
#'     CellLineName = "name of the cell line")
#' )
#'
#' @return data.table with metadata with annotated cell lines
#' @export
#'
add_CellLine_annotation <- function(
    df_metadata,
    DB_cellid_header = "cell_line_identifier",
    DB_cell_annotate = c(
      "cell_line_name", 
      "primary_tissue", 
      "doubling_time",
      "parental_identifier", 
      "subtype"
    ),
    fname = "cell_lines.csv",
    fill = "unknown",
    annotationPackage = if (requireNamespace("gDRinternalData", quietly = TRUE)) {
      "gDRinternalData"
    } else {
      "gDRtestData"
    }
) {
  
  # Assertions:
  checkmate::assert_data_table(df_metadata)
  checkmate::assert_string(fill, null.ok = TRUE)
  
  cellline <- gDRutils::get_env_identifiers("cellline")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  add_clid <- gDRutils::get_header("add_clid")
  
  CLs_info <- data.table::fread(
    system.file("annotation_data", fname, package = annotationPackage), header = TRUE
  )
  CLs_info <- CLs_info[, c(DB_cellid_header, DB_cell_annotate), with = FALSE]
  
  if (nrow(CLs_info) == 0) return(df_metadata)
  
  validatedCLs <- unique(df_metadata[[cellline]]) %in% CLs_info[[DB_cellid_header]]
  missingTblCellLines <- NULL
  if (!is.null(fill) && !all(validatedCLs)) {
    unValidatedCellLine <- unique(df_metadata[[cellline]])[!validatedCLs]
    missingTblCellLines <- data.table::data.table(
      parental_identifier = unValidatedCellLine,
      cell_line_name = unValidatedCellLine,
      cell_line_identifier = unValidatedCellLine,
      doubling_time = NA,
      primary_tissue = fill,
      subtype = fill
    )
  }
  
  if (any(!df_metadata[[cellline]] %in% CLs_info[[DB_cellid_header]]) && 
      !is.null(missingTblCellLines)) {
    CLs_info <- rbind(CLs_info, missingTblCellLines)
  }
  cols_to_fill <- names(CLs_info)[!names(CLs_info) %in% c("doubling_time")]
  for (col in cols_to_fill) {
    CLs_info[[col]][is.na(CLs_info[[col]])] <- fill
  }
  cols_num <- "doubling_time"
  CLs_info[, (cols_num) := as.numeric(get(cols_num))]
  
  
  colnames(CLs_info) <- unlist(c(cellline, add_clid, tail(DB_cell_annotate, 2)))
  
  futile.logger::flog.info("Merge with Cell line info")
  nrows_df <- nrow(df_metadata)
  df_metadata <- base::merge(df_metadata, CLs_info, by = cellline, all.x = TRUE)
  stopifnot(nrows_df == nrow(df_metadata))
  df_metadata
}


#' add_Drug_annotation
#'
#' add drug annotation to a data.table with metadata
#'
#' @param df_metadata data.table with metadata
#' @param fname string with file name with annotation
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @param annotationPackage string indication name of the package containing drug annotation
#' @details The logic of adding drug annotation for df_metadata 
#' based on the annotation file stored in gDRtestData.
#' @examples
#' add_Drug_annotation(
#'   data.table::data.table(
#'     Gnumber = "drug_id", 
#'     DrugName = "name of the drug")
#' )
#' @return data.table with metadata with annotated drugs
#' @export
#'
add_Drug_annotation <- function(
    df_metadata,
    fname = "drugs.csv",
    fill = "unknown",
    annotationPackage = if (requireNamespace("gDRinternalData", quietly = TRUE)) {
      "gDRinternalData"
    } else {
      "gDRtestData"
    }
) {
  
  # Assertions:
  checkmate::assert_data_table(df_metadata)
  checkmate::assert_string(fill, null.ok = TRUE)
  nrows_df <- nrow(df_metadata)
  
  drug <- unlist(gDRutils::get_env_identifiers(c(
    "drug", "drug2", "drug3"), simplify = FALSE))
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  drug_name <- unlist(gDRutils::get_env_identifiers(c(
    "drug_name", "drug_name2", "drug_name3"), simplify = FALSE))
  drug_moa <- unlist(gDRutils::get_env_identifiers(c(
    "drug_moa", "drug_moa2", "drug_moa3"), simplify = FALSE))
  drug_idx <- which(drug %in% names(df_metadata))
  drugsTreated <- unique(unlist(subset(df_metadata, select = drug[drug_idx])))
  
  drug_full_identifiers <- c(
    drug[drug_idx], 
    drug_name[drug_idx], 
    drug_moa[drug_idx]
  )
  drug_ids <- stringr::str_extract(names(drug_full_identifiers), "[0-9]") 
  drug_ids <- ifelse(is.na(drug_ids), 1, drug_ids)
  drug_identifiers_list <- split(drug_full_identifiers, drug_ids)
  names(drug_identifiers_list) <- drug[drug_idx]
  
  # Read local drugs annotations
  Drug_info <- data.table::fread(
    system.file("annotation_data", fname, package = annotationPackage), header = TRUE
  )
  Drug_info <- Drug_info[, c("gnumber", "drug_name", "drug_moa"), with = FALSE]
  data.table::setnames(Drug_info, c("gnumber", "drug_name", "drug_moa"),
                       c("drug", "drug_name", "drug_moa"))
  drugsTreated <- drugsTreated[!drugsTreated %in% untreated_tag]
  validatedDrugs <- 
    remove_drug_batch(drugsTreated) %in% remove_drug_batch(Drug_info[["drug"]])
  #### function should be parallelizeized
  missingTblDrugs <- NULL
  if (!is.null(fill) && any(!validatedDrugs)) {
    missingTblDrugs <- data.table::data.table(
      drug = remove_drug_batch(drugsTreated[!validatedDrugs]),
      drug_name = drugsTreated[!validatedDrugs],
      drug_moa = fill
    )
  }
  
  if (nrow(Drug_info) == 0) {
    drug_name <- intersect(drug_name, colnames(df_metadata))
    drug <- intersect(drug, colnames(df_metadata))
    df_metadata[, (drug_name) := get(drug)]
    return(df_metadata)
  }
  
  Drug_info <- rbind(Drug_info, missingTblDrugs)
  
  Drug_info$drug <- remove_drug_batch(Drug_info$drug)
  Drug_info <-
    rbind(data.table::data.table(
      drug = untreated_tag,
      drug_name = untreated_tag,
      drug_moa = untreated_tag
    ),
    Drug_info)
  Drug_info <- Drug_info[!duplicated(Drug_info[["drug"]]), ]
  if (any(!remove_drug_batch(drugsTreated) %in% Drug_info$drug) && 
      !is.null(missingTblDrugs)) {
    Drug_info <- rbind(Drug_info, stats::setNames(
      missingTblDrugs[!(remove_drug_batch(missingTblDrugs$drug) %in% 
                          Drug_info$drug), ],
      names(Drug_info)
    ))
  }
  for (drug_idf in names(drug_identifiers_list)) {
    colnames(Drug_info) <- drug_identifiers_list[[drug_idf]]
    df_metadata$batch <- df_metadata[[drug_idf]]
    df_metadata[[drug_idf]] <- remove_drug_batch(df_metadata[[drug_idf]])
    req_col <- c(drug_idf, setdiff(colnames(df_metadata), colnames(Drug_info)))
    df_metadata <- merge(df_metadata[, req_col, with = FALSE], Drug_info, by = drug_idf, all.x = TRUE)
    df_metadata[[drug_idf]] <- df_metadata$batch
    df_metadata$batch <- NULL
  }
  stopifnot(nrows_df == nrow(df_metadata))
  df_metadata
}

#' Remove batch from Gnumber
#' 
#' @param drug drug name
#'
#' @examples
#' remove_drug_batch("DRUG.123")
#'
#' @return Gnumber without a batch
#' @export
remove_drug_batch <- function(drug) {
  gsub("\\.[0-9]+.*", "", drug)
}
