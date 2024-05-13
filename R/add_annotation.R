#' add_CellLine_annotation
#'
#' add cellline annotation to a data.table with metadata
#'
#' @param dt_metadata data.table with metadata
#' @param DB_cellid_header string with colnames with cell line identifier
#'                         in the annotation file
#' @param DB_cell_annotate character vector with mandatory colnames 
#'                         used in the annotation file
#' @param fname string with file name with annotation
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @param annotationPackage string indication name of the package containing cellline annotation
#' @param externalSource string with path to external file with annotation data; by default it checks 
#' 'GDR_CELLLINE_ANNOTATION' env var. This file should contain columns such as gnumber, drug_name and drug_moa
#' @keywords annotation
#' @details
#' The logic of adding celline annotation for dt_metadata based on 
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
    dt_metadata,
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
    annotationPackage = if ("gDRinternal" %in% .packages(all.available = TRUE)) {
      "gDRinternal"
    } else {
      "gDRtestData"
    },
    externalSource = Sys.getenv("GDR_CELLLINE_ANNOTATION")
) {
  
    # Assertions:
    checkmate::assert_data_table(dt_metadata)
    checkmate::assert_string(fill, null.ok = TRUE)
    checkmate::assert_string(externalSource, null.ok = TRUE)
    
    cellline <- gDRutils::get_env_identifiers("cellline")
    cellline_name <- gDRutils::get_env_identifiers("cellline_name")
    add_clid <- gDRutils::get_header("add_clid")
    
    if (all(c(cellline, cellline_name, add_clid) %in% names(dt_metadata))) {
      return(dt_metadata)
    }
    
    CLs_info <- if (nchar(externalSource) && file.exists(externalSource)) {
      data.table::fread(externalSource)
    } else if (annotationPackage == "gDRtestData") {
      data.table::fread(
        system.file("annotation_data", fname, package = annotationPackage), header = TRUE
      )
    } else {
      eval(parse(text = paste0(annotationPackage,
                               "::",
                               "get_cell_line_annotations")))(dt_metadata[[cellline]])
    }
  
    CLs_info <- CLs_info[, c(DB_cellid_header, DB_cell_annotate), with = FALSE]
  
    validatedCLs <- unique(dt_metadata[[cellline]]) %in% CLs_info[[DB_cellid_header]]
    missingTblCellLines <- NULL
    if (!is.null(fill) && !all(validatedCLs)) {
      unValidatedCellLine <- unique(dt_metadata[[cellline]])[!validatedCLs]
      missingTblCellLines <- data.table::data.table(
        parental_identifier = unValidatedCellLine,
        cell_line_name = unValidatedCellLine,
        cell_line_identifier = unValidatedCellLine,
        doubling_time = NA,
        primary_tissue = fill,
        subtype = fill
      )
    }
    
    if (any(!dt_metadata[[cellline]] %in% CLs_info[[DB_cellid_header]]) && 
        !is.null(missingTblCellLines)) {
      CLs_info <- rbind(CLs_info, missingTblCellLines)
    }
    cols_to_fill <- names(CLs_info)[!names(CLs_info) %in% c("doubling_time")]
    for (col in cols_to_fill) {
      CLs_info[[col]][is.na(CLs_info[[col]])] <- fill
    }
    cols_num <- "doubling_time"
    CLs_info[, (cols_num) := as.numeric(get(cols_num))]
    
    
    colnames(CLs_info) <- unlist(c(cellline, add_clid[c("cellline_name",
                                                        "cellline_tissue",
                                                        "cellline_ref_div_time",
                                                        "cellline_parental_identifier",
                                                        "cellline_subtype")]))
    
    futile.logger::flog.info("Merge with Cell line info")
    nrows_df <- nrow(dt_metadata)
    dt_metadata <- CLs_info[dt_metadata, on = cellline]
    stopifnot(nrows_df == nrow(dt_metadata))
    dt_metadata
  }


#' add_Drug_annotation
#'
#' add drug annotation to a data.table with metadata
#'
#' @param dt_metadata data.table with metadata
#' @param fname string with file name with annotation
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @param annotationPackage string indication name of the package containing drug annotation
#' @param externalSource string with path to external file with annotation data; by default it checks 
#' 'GDR_DRUG_ANNOTATION' env var. This file should contain columns such as gnumber, drug_name, and drug_moa
#' @keywords annotation
#' @details The logic of adding drug annotation for dt_metadata 
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
    dt_metadata,
    fname = "drugs.csv",
    fill = "unknown",
    annotationPackage = if ("gDRinternal" %in% .packages(all.available = TRUE)) {
      "gDRinternal"
    } else {
      "gDRtestData"
    },
    externalSource = Sys.getenv("GDR_DRUG_ANNOTATION")
) {
  
  # Assertions:
  checkmate::assert_data_table(dt_metadata)
  checkmate::assert_string(fill, null.ok = TRUE)
  checkmate::assert_string(externalSource, null.ok = TRUE)
  
  nrows_df <- nrow(dt_metadata)
  
  drug <- unlist(gDRutils::get_env_identifiers(c(
    "drug", "drug2", "drug3"), simplify = FALSE))
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  drug_name <- unlist(gDRutils::get_env_identifiers(c(
    "drug_name", "drug_name2", "drug_name3"), simplify = FALSE))
  drug_moa <- unlist(gDRutils::get_env_identifiers(c(
    "drug_moa", "drug_moa2", "drug_moa3"), simplify = FALSE))
  drug_idx <- which(drug %in% names(dt_metadata))
  drugsTreated <- unique(unlist(subset(dt_metadata, select = drug[drug_idx])))
  
  if (all(c(drug[["drug"]],
            drug_name[["drug_name"]],
            drug_moa[["drug_moa"]]) %in% names(dt_metadata))) {
    return(dt_metadata)
  } else {
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
    Drug_info <- if (nchar(externalSource) && file.exists(externalSource)) {
      data.table::fread(externalSource)
    } else if (annotationPackage == "gDRtestData") {
      data.table::fread(
        system.file("annotation_data", fname, package = annotationPackage), header = TRUE
      )
    } else {
      eval(parse(text = paste0(annotationPackage,
                               "::",
                               "get_drug_annotations")))()
    }
      
      
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
      drug_name <- intersect(drug_name, colnames(dt_metadata))
      drug <- intersect(drug, colnames(dt_metadata))
      dt_metadata[, (drug_name) := get(drug)]
      return(dt_metadata)
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
      dt_metadata$batch <- dt_metadata[[drug_idf]]
      dt_metadata[[drug_idf]] <- remove_drug_batch(dt_metadata[[drug_idf]])
      req_col <- c(drug_idf, setdiff(colnames(dt_metadata), colnames(Drug_info)))
      dt_metadata <- Drug_info[dt_metadata[, req_col, with = FALSE],  on = drug_idf]
      dt_metadata[[drug_idf]] <- dt_metadata$batch
      dt_metadata$batch <- NULL
    }
    stopifnot(nrows_df == nrow(dt_metadata))
    dt_metadata
  }
}

#' Remove batch from Gnumber
#' 
#' @param drug drug name
#' @keywords annotation
#'
#' @examples
#' remove_drug_batch("DRUG.123")
#'
#' @return Gnumber without a batch
#' @export
remove_drug_batch <- function(drug) {
  gsub("\\.[0-9]+.*", "", drug)
}


#' Retrieve the drug annotation from the annotated dt input
#'
#' @param dt annotated data.table
#'
#' @return data.table with drug annotation
#' @export
#'
#' @examples
#' dt <- data.table::data.table(Gnumber = "A",
#' DrugName = "drugA",
#' drug_moa = "drug_moa_A")
#' get_drug_annotation_from_dt(dt)
get_drug_annotation_from_dt <- function(dt) {
  checkmate::assert_data_table(dt)
  drug_cols <- intersect(gDRutils::get_env_identifiers()
                         [grep("drug",names(gDRutils::get_env_identifiers()))],
                         names(dt))
  dt_drug <- dt[, unlist(drug_cols), with = FALSE]
  dt_long <- data.table::melt(dt_drug,
                              measure.vars = patterns(paste0("^", unlist(drug_cols[c("drug",
                                                                                     "drug_name",
                                                                                     "drug_moa")]))),
                              value.name = unlist(drug_cols[c("drug",
                                                              "drug_name",
                                                              "drug_moa")]))
  dt_long[, variable := NULL]
  data.table::setnames(dt_long,
                       unlist(drug_cols[c("drug",
                                          "drug_name",
                                          "drug_moa")]),
                       c("gnumber", "drug_name", "drug_moa"))
  unique_dt <- unique(dt_long)
  unique_dt[!gnumber %in% gDRutils::get_env_identifiers("untreated_tag"), ]
}


#' Retrieve the cell line annotation from the annotated dt input
#'
#' @param dt annotated data.table
#'
#' @return data.table with cell line annotation
#' @export
#'
#' @examples
#' dt <- data.table::data.table(Gnumber = "A",
#' clid = "CL123",
#' CellLineName = "cl name",
#' Tissue = "Bone",
#' parental_identifier = "some cl",
#' subtype = "cortical",
#' ReferenceDivisionTime = 5)
#' get_cellline_annotation_from_dt(dt)
get_cellline_annotation_from_dt <- function(dt) {
  checkmate::assert_data_table(dt)
  cell_cols <- c(gDRutils::get_env_identifiers("cellline"),
                 gDRutils::get_header("add_clid"))
  cell_dt <- dt[, unlist(cell_cols), with = FALSE]
  data.table::setnames(cell_dt,
                       unlist(cell_cols),
                       c("cell_line_identifier",
                         "cell_line_name",
                         "primary_tissue",
                         "parental_identifier",
                         "subtype",
                         "doubling_time"))
  unique(cell_dt)
}
