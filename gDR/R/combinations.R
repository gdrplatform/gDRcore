#' Add codrug group
#'
#' @param se 
#'
#' @return
#' @export
#'
# TODO: Utilize the set_SE_metadata functions. 
add_codrug_group_SE <- function(se) {

  r_data <- SummarizedExperiment::rowData(se)
  if (!(paste0(gDRutils::get_identifier("drugname"), '_2') %in% colnames(r_data))) return(se)

  # find the pairs of drugs with relevant metadata
  drug_ids <- paste0(gDRutils::get_identifier("drugname"), c('', '_2'))
  other_metadata <- c(paste0(gDRutils::get_identifier("drug"), c('', '_2')),
            setdiff(colnames(r_data), c('Concentration_2', drug_ids,
                paste0(gDRutils::get_identifier("drug"), c('', '_2')),
                paste0(gDRutils::get_identifier("drug_moa"), c('', '_2')))))
  drug_pairs <- unique(r_data[, c(drug_ids, other_metadata)])
  drug_pairs <- drug_pairs[ !(drug_pairs[,drug_ids[2]] %in% gDRutils::get_identifier('untreated_tag')),]

  pair_list <- vector('list', nrow(drug_pairs))
  # loop through the pairs to assess the number of individual concentration pairs
  for (idp in 1:nrow(drug_pairs)) {
    row_idx <- r_data[,drug_ids[1]] %in% unlist(drug_pairs[idp, drug_ids]) &
            r_data[,drug_ids[2]] %in% c(unlist(drug_pairs[idp, drug_ids]),
                gDRutils::get_identifier('untreated_tag')) &
            apply(as.matrix(
                IRanges::LogicalList(c(
                  lapply(setdiff(other_metadata,
                      paste0(gDRutils::get_identifier("drug"), c('', '_2'))),
                    function(y) # matching the metadata
                    r_data[, y] == drug_pairs[idp, y])
                  ))), 2, all)

    # reverse engineer the type of combination experiment
    flat_data <- gDRutils::assay_to_dt(se[row_idx, ], 'Averaged')
    flat_data <- flat_data[flat_data$Concentration_2 > 0,]
    conc_1 <- table(flat_data$Concentration)
    conc_2 <- table(flat_data$Concentration_2)
    n_conc_pairs <- nrow(unique(flat_data[,c('Concentration', 'Concentration_2')]))
    conc_ratio <- table(round(log10(flat_data$Concentration / flat_data$Concentration_2), 2))
    conc_ratio <- conc_ratio[!names(conc_ratio) %in% c('Inf', '-Inf')]

    condition <- paste(paste(other_metadata, unlist(drug_pairs[idp,other_metadata]), sep = '='),
                  collapse=' ')
    if (length(conc_ratio) <= 2) {
      type <- 'co-dilution'
      print(sprintf('Found %s combination with %s and %s: ratio of %.2f, %i concentrations (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp, 2], 10**as.numeric(names(conc_ratio)),
            length(conc_2), condition))
      condition = c(as.list(drug_pairs[idp,]), list(conc_ratio = 10**as.numeric(names(conc_ratio))))
    } else if (n_conc_pairs == length(conc_1)*length(conc_2) & length(conc_2) >= 4) {
      type <- 'matrix'
      print(sprintf('Found %s combination with %s and %s: %i x %i concentrations (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], length(conc_1), length(conc_2), condition))
      condition = c(as.list(drug_pairs[idp,]), list(
            Concentration = as.numeric(names(conc_1)), Concentration_2 = as.numeric(names(conc_2))))
    } else if (length(conc_2)<4) {
      type <- 'fixed'
      print(sprintf('Found %s combination of %s with %s at %.3g uM (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], as.numeric(names(conc_2)), condition))
      condition = c(as.list(drug_pairs[idp,]), list( Concentration_2 = as.numeric(names(conc_2))))
    } else {
      type <- 'other'
      print(sprintf('Found %s combination with %s and %s: %i concentration pairs (%s)',
        type, drug_pairs[idp,1], drug_pairs[idp,2], n_conc_pairs, condition))
      condition = c(as.list(drug_pairs[idp,]), list(
            Concentration = as.numeric(names(conc_1)), Concentration_2 = as.numeric(names(conc_2))))
    }

    pair_list[[idp]] <- list(condition = condition,
                          rows = rownames(r_data)[row_idx],
                          type = type)
  }

  metadata(se)$drug_combinations <- pair_list
  return(se)
}
