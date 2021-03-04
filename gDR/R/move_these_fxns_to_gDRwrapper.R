#' add_CellLine_annotation
#'
#' add cellline annotation to a data.frame with metadata
#'
#' @param df_metadata a data.frame with metadata
#' @param fill_DB_with_unknown a logical indicating whether DB should be filled with unknown cell lines
#' 
#' @return a data.frame with metadata with annotated cell lines
#' @export
#'
# TODO: Move this to gDRwrapper.
add_CellLine_annotation <- function(df_metadata,
                                    fill_DB_with_unknown = FALSE) {
  
    # Assertions:
    stopifnot(inherits(df_metadata, "data.frame"))
    checkmate::assert_logical(fill_DB_with_unknown)
    
    DB_cellid_header <- "cell_line_identifier"
    DB_cell_annotate <- c("cell_line_name", "primary_tissue", "doubling_time", "parental_identifier", "subtype")
    # corresponds to columns gDRutils::get_header("add_clid"): name, tissue, doubling time
    
    # the logic of adding celline annotation for df_metadata is based on the function get_cell_lines from the gDRwrapper
    # we added additional parameter 'fill_DB_with_unknown' that allows to fill the DB with clid info for these cell lines
    # that are not present in the DB.
    # Other fields are set as "UNKNOWN". If the fill_DB_with_unknown is set as FALSE we add unknown cell lines
    # only to the tibble.
    # This approach will be corrected once we will implement final solution for adding cell lines.
    validateCLs <- gDRwrapper::validate_cell_lines(unique(df_metadata[[gDRutils::get_identifier("cellline")]]))
    if(!validateCLs){
      missingTblCellLines <- tibble::tibble(parental_identifier = "UNKNOWN",
                                            cell_line_name = "UNKNOWN",
                                            cell_line_identifier = unique(df_metadata[,gDRutils::get_identifier("cellline")]),
                                            doubling_time = "UNKNOWN",
                                            primary_tissue = "UNKNOWN",
                                            subtype = "UNKNOWN")
      
      if(fill_DB_with_unknown){
        addMissingCellLines <- gDRwrapper::add_drugs(missingTblCellLines)
      }
    }
    CLs_info <- tryCatch( {
        CLs_info <- gDRwrapper::get_cell_lines()
        CLs_info <- CLs_info[CLs_info$cell_line_identifier %in% unique(df_metadata[[gDRutils::get_identifier("cellline")]]), ]
        CLs_info <- CLs_info[, c(DB_cellid_header, DB_cell_annotate), with = FALSE]
        CLs_info
    }, error = function(e) {
      futile.logger::flog.error("Failed to load cell line info from DB: %s", e)
        data.frame()
    })

    if (nrow(CLs_info) == 0) return(df_metadata)

    colnames(CLs_info)[1:4] <- c(gDRutils::get_identifier("cellline"), gDRutils::get_header("add_clid"))
    CLIDs <- unique(df_metadata[[gDRutils::get_identifier("cellline")]])
    bad_CL <- CLs_info[gDRutils::get_identifier("cellline") %in% CLIDs][[gDRutils::get_identifier("cellline")]]
    if (any(bad_CL)) {
        futile.logger::flog.warn("Cell line ID %s not found in cell line database",
                     paste(CLIDs[bad_CL], collapse = " ; "))
        temp_CLIDs = data.frame(CLIDs[bad_CL], CLIDs[bad_CL])
        temp_CLIDs[, 1+(2:length(gDRutils::get_header("add_clid")))] = NA
        colnames(temp_CLIDs) = c(gDRutils::get_identifier("cellline"),
                      gDRutils::get_header("add_clid"))
        CLs_info <- rbind(CLs_info, temp_CLIDs)
        }

    futile.logger::flog.info("Merge with Cell line info")
    nrows_df <- nrow(df_metadata)
    df_metadata <- base::merge(df_metadata, CLs_info, by = gDRutils::get_identifier("cellline"), all.x = TRUE)
    stopifnot(nrows_df == nrow(df_metadata))
    return(df_metadata)

}


#' add_Drug_annotation
#'
#' add drug annotation to a data.frame with metadata
#'
#' @param df_metadata a data.frame with metadata
#' @param fill_DB_with_unknown a logical indicating whether DB should be filled with unknown drugs
#' 
#'
#' @return a data.frame with metadata with annotated drugs
#' @export
#'
# TODO: Move this to gDRwrapper.
add_Drug_annotation <- function(df_metadata,
                                fill_DB_with_unknown = FALSE) {
  
        # Assertions:
        stopifnot(inherits(df_metadata, "data.frame"))
        checkmate::assert_logical(fill_DB_with_unknown)
  
        nrows_df <- nrow(df_metadata)

        DB_drug_identifier <- c("gnumber", "drug_name", "drug_moa")
        # the logic of adding drug annotation for df_metadata is based on the function get_drugs from the gDRwrapper
        # we added additional parameter 'fill_DB_with_unknown' that allows to fill the DB with drug_name and gnumber, for these drugs,
        # that are not present in the DB
        # Other fields are set as "UNKNOWN". If the fill_DB_with_unknown is set as FALSE we add unkonown cell lines
        # only to the tibble.
        # This approach will be corrected once we will implement final solution for adding cell lines.

        drugsTreated <- unique(df_metadata[[gDRutils::get_identifier("drug")]])
        
        drugsTreated <- drugsTreated[!drugsTreated %in% gDRutils::get_identifier("untreated_tag")]
        validateDrugs <- gDRwrapper::validate_drugs(drugsTreated)
        if(!validateDrugs){
          missingTblDrugs <- tibble::tibble(drug_name = drugsTreated,
                                            drug_moa = "UNKNOWN",
                                            gnumber = drugsTreated)
          if(fill_DB_with_unknown){
            addMissingDrugs <- gDRwrapper::add_drugs(missingTblDrugs)
          }
          
        }
        Drug_info <- tryCatch({
          # TODO: refactor this part of code once we switch to DataFrameMatrix class
          gDrugs <- gDRwrapper::get_drugs()[, c(..DB_drug_identifier)]
          #gDrugs[, 1] <- gsub("\\..*", "", gDrugs$gnumber) # remove batch number from DB_drug_identifier
          gDrugs
        }, error = function(e) {
          futile.logger::flog.error("Failed to load drug info from DB: %s", e)
            data.frame()
        })

        if (nrow(Drug_info) == 0) {
            df_metadata[, gDRutils::get_identifier("drugname")] = df_metadata[, gDRutils::get_identifier("drug")]
            return(df_metadata)
        }
        
        # -----------------------

        colnames(Drug_info)[1:2] <- c("drug", "drug_name")
        Drug_info <-
          rbind(data.frame(
            drug = gDRutils::get_identifier("untreated_tag"),
            drug_name = gDRutils::get_identifier("untreated_tag"),
            drug_moa = gDRutils::get_identifier("untreated_tag")
          ),
          Drug_info)
        Drug_info <- Drug_info[!duplicated(Drug_info[["drug"]]),]
        DrIDs <- unique(unlist(df_metadata[,grep(gDRutils::get_identifier("drug"), colnames(df_metadata)), with = FALSE]))
        if(any(!gsub("\\..*", "", drugsTreated) %in% Drug_info$drug) & exists("missingTblDrugs")){
          Drug_info <- rbind(Drug_info, data.table::setnames(missingTblDrugs[!drugsTreated %in% Drug_info$drug, ], names(Drug_info)))
        }
        bad_DrID <- !(gsub("\\..*", "", DrIDs) %in% Drug_info$drug) & !is.na(DrIDs)
        if (any(bad_DrID)) {
            # G number, but not registered
            ok_DrID <- attr(regexpr("^G\\d*",DrIDs), "match.length")==9
            if (any(ok_DrID)) {
              futile.logger::flog.warn("cleanup_metadata: Drug %s  not found in gCSI database; use G# as DrugName",
                                       paste(DrIDs[ok_DrID & bad_DrID], collapse = " ; "))
              Drug_info <-
                rbind(Drug_info, data.frame(drug = DrIDs[ok_DrID & bad_DrID],
                                            DrugName = DrIDs[ok_DrID & bad_DrID]))
            } else {
              futile.logger::flog.error("Drug %s not in the correct format for database",
                  paste(DrIDs[!ok_DrID], collapse = ' ; '))
            }
        }
        colnames(Drug_info)[2] <- gDRutils::get_identifier("drugname")
        futile.logger::flog.info("Merge with Drug_info for Drug 1")
        df_metadata[[paste0(gDRutils::get_identifier("drug"), "_temp")]] <- gsub("\\..*", "", df_metadata[[gDRutils::get_identifier("drug")]])
        df_metadata <- base::merge(df_metadata, Drug_info, by.x = gsub("\\..*", "", paste0(gDRutils::get_identifier("drug"), "_temp")), by.y = "drug", all.x = TRUE)
        df_metadata <- df_metadata[, .SD, .SDcols = !endsWith(names(df_metadata), "temp")]
        # add info for columns Gnumber_*
        for (i in grep(paste0(gDRutils::get_identifier("drug"),"_\\d"), colnames(df_metadata))) {
            df_metadata[is.na(df_metadata[[i]]), i] = gDRutils::get_identifier("untreated_tag")[1] # set missing values to Untreated
            Drug_info_ <- Drug_info
            colnames(Drug_info_)[2:length(colnames(Drug_info_))] <- paste0(colnames(Drug_info_)[c(2,3)], substr(colnames(df_metadata)[i], 8, 12))
            futile.logger::flog.info("Merge with Drug_info for %s", colnames(df_metadata)[[i]])
            df_metadata[[paste0(colnames(df_metadata)[[i]], "_temp")]] <- gsub("\\..*", "", df_metadata[[i]])
            df_metadata <- base::merge(df_metadata, Drug_info_, by.x = gsub("\\..*", "", paste0(colnames(df_metadata)[[i]], "_temp")), by.y = "drug", all.x = TRUE)
            df_metadata <- df_metadata[, .SD, .SDcols = !endsWith(names(df_metadata), "temp")]
        }
    stopifnot(nrows_df == nrow(df_metadata))

    return(df_metadata)
}


#' update_metadata
#'
#' Update metadata in a SummarizedExperiment
#'
#' @param SE a SummarizedExperiment
#'
#' @return a SummarizedExperiment with additional metadata
#' @export
#'
update_experiment_metadata = function(metadata,
                                      metadataList = NULL) {
  # Assertions:
  checkmate::assert_list(metadata)
  if(length(metadata$experiment_metadata)>0 && length(metadata$experiment_metadata$name)>0) {
    description <- as.character(metadata$experiment_metadata$description)
    experiment_name <- as.character(metadata$experiment_metadata$name)
    expert_unixid <- as.character(metadata$experiment_metadata$experimentalist)
    qcs_id <- tryCatch(as.character(strsplit(as.character(metadata$experiment_metadata$name), " ")[[1]][1]),
                       error=function(cond) {
                         NULL })
  }
  checkmate::assert_true(all(names(metadataList) %in% c("expert_unixid", "unix_id", "description",
                                                    "assay_id", "date_experiment",
                                                    "source_id", "state_id",
                                                    "experiment_name", "experiment_project", "qcs_id",
                                                    "labhead_unixid")))
  if(!is.null(metadataList)) {
    for (element in names(metadataList)){
      assign(element, metadataList[[element]])
    }
  }
  
  if (!exists("source_id")) {
    source_id <-
      subset(gDRwrapper::get_source_keys(), source_type == "UNKNOWN")[["source_id"]]
  }
  
  if (!exists("assay_id")) {
    assay_id <-
      subset(gDRwrapper::get_assay_types(), assay_type == "UNKNOWN")[["assay_id"]]
  }
  
  if (!exists("state_id")) {
    state_id <-
      tryCatch({
        subset(gDRwrapper::get_state_types(), state_name == "active")[["state_id"]]
      },
      error=function(cond) {
        return(1)
      })
  }
  
  metadata$experiment_metadata <- tibble::tibble(
    expert_unixid = ifelse(exists("expert_unixid"), expert_unixid, Sys.getenv("USER")),
    unix_id = ifelse(exists("unix_id"), unix_id, Sys.getenv("USER")),
    description = ifelse(exists("description"), description, NA),
    assay_id = ifelse(exists("assay_id"), assay_id, NA),
    date_experiment = ifelse(exists("date_experiment"), date_experiment, as.character(Sys.Date())),
    source_id = ifelse(exists("source_id"), source_id, NA),
    state_id = ifelse(exists("state_id"), state_id, NA),
    experiment_name = ifelse(exists("experiment_name"), experiment_name, NA),
    experiment_project = ifelse(exists("experiment_project"), experiment_project, NA),
    qcs_id = ifelse(exists("qcs_id"), qcs_id, NA),
    labhead_unixid = ifelse(exists("labhead_unixid"), labhead_unixid, NA),
  )
  return(metadata)
}
