#' update_metadata
#'
#' Update metadata in a SummarizedExperiment
#'
#' @param SE a SummarizedExperiment
#'
#' @return a SummarizedExperiment with additional metadata
#' @export
#'

update_metadata = function(SE,
                           expert_unixid = Sys.getenv("USER"),
                           description = NA,
                           assay_id = NA,
                           date_processed = NA,
                           source_id = NA,
                           state_id = NA,
                           experiment_name = NA,
                           qcs_id = NA,
                           labhead_unixid = NA) {
  # Assertions:
  checkmate::assert_class(SE, "SummarizedExperiment")
  if(length(S4Vectors::metadata(SE)$experiment_metadata)>0) {
    description <- as.character(metadata(SE)$experiment_metadata$description)
    experiment_name <- as.character(metadata(SE)$experiment_metadata$name)
    expert_unixid <- as.character(metadata(SE)$experiment_metadata$experimentalist)
    qcs_id <- as.character(strsplit(as.character(metadata(SE)$experiment_metadata$name), " ")[[1]][1])
  }
  metadata(SE)$experiment_metadata <- data.frame(
    expert_unixid = expert_unixid,
    description = description,
    assay_id = assay_id,
    date_processed = date_processed,
    source_id = source_id,
    state_id = state_id,
    experiment_name = experiment_name,
    qcs_id = qcs_id,
    labhead_unixid = labhead_unixid
  )
  return(SE)
}

