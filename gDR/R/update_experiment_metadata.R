#' update_metadata
#'
#' Update metadata in a SummarizedExperiment
#'
#' @param SE a SummarizedExperiment
#'
#' @return a SummarizedExperiment with additional metadata
#' @export
#'

update_experiment_metadata = function(
                           metadata,
                           metadataList) {
  # Assertions:
  checkmate::assert_list(metadata)
  if(length(metadata$experiment_metadata)>0 && exists('metadata$experiment_metadata$name')) {
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
                                                    "experiment_name", "qcs_id",
                                                    "labhead_unixid")))
  for (element in names(metadataList)){
    assign(element, metadataList[[element]])
  }
  
  metadata$experiment_metadata <- data.frame(
    expert_unixid = ifelse(exists("expert_unixid"), expert_unixid, Sys.getenv("USER")),
    unix_id = ifelse(exists("unix_id"), unix_id, Sys.getenv("USER")),
    description = ifelse(exists("description"), description, NA),
    assay_id = ifelse(exists("assay_id"), assay_id, NA),
    date_experiment = ifelse(exists("date_experiment"), date_experiment, NA),
    source_id = ifelse(exists("source_id"), source_id, NA),
    state_id = ifelse(exists("state_id"), state_id, 1),
    experiment_name = ifelse(exists("experiment_name"), experiment_name, NA),
    qcs_id = ifelse(exists("qcs_id"), qcs_id, NA),
    labhead_unixid = ifelse(exists("labhead_unixid"), labhead_unixid, NA),
    stringsAsFactors = FALSE
  )
  return(metadata)
}

