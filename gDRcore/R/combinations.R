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
  if (!(paste0(gDRutils::get_env_identifiers("drugname"), "_2") %in% colnames(r_data)) ||
     all(r_data[[paste0(gDRutils::get_env_identifiers("drugname"), "_2")]] %in%
       gDRutils::get_env_identifiers("untreated_tag"))) return(se)

  # find the pairs of drugs with relevant metadata
  drug_ids <- paste0(gDRutils::get_env_identifiers("drugname"), c("", "_2"))
  other_metadata <- c(#paste0(gDRutils::get_env_identifiers("drug"), c("", "_2")),
            setdiff(colnames(r_data), c("Concentration_2", drug_ids,
                paste0(gDRutils::get_env_identifiers("drug"), c("", "_2")),
                paste0(gDRutils::get_env_identifiers("drug_moa"), c("", "_2")))))
  drug_pairs <- unique(r_data[, c(drug_ids, other_metadata)])
  drug_pairs <- drug_pairs[!(drug_pairs[, drug_ids[2]] %in% gDRutils::get_env_identifiers("untreated_tag")), ]

  pair_list <- vector("list", nrow(drug_pairs))
  # loop through the pairs to assess the number of individual concentration pairs
  for (idp in 1:nrow(drug_pairs)) {
    row_idx <- r_data[, drug_ids[1]] %in% unlist(drug_pairs[idp, drug_ids]) &
            r_data[, drug_ids[2]] %in% c(unlist(drug_pairs[idp, drug_ids]),
                gDRutils::get_env_identifiers("untreated_tag")) &
            apply(as.matrix(
                IRanges::LogicalList(c(
                  lapply(other_metadata,
                    function(y) # matching the metadata
                    r_data[, y] == drug_pairs[idp, y])
                  ))), 2, all)

    # reverse engineer the type of combination experiment
    flat_data <- gDRutils::convert_se_assay_to_dt(se[row_idx, ], "Averaged")
    flat_data <- flat_data[flat_data$Concentration_2 > 0, ]
    conc_1 <- table(flat_data$Concentration)
    conc_2 <- table(flat_data$Concentration_2)
    n_conc_pairs <- nrow(unique(flat_data[, c("Concentration", "Concentration_2")]))
    conc_ratio <- table(round(log10(flat_data$Concentration / flat_data$Concentration_2), 2))
    conc_ratio <- conc_ratio[!names(conc_ratio) %in% c("Inf", "-Inf")]

    condition_name <- paste(paste(other_metadata, unlist(drug_pairs[idp, other_metadata]), sep = "="),
                  collapse = " ")
    if (length(conc_ratio) <= 2) {
      type <- "co-dilution"
      print(sprintf("Found %s combination with %s and %s: ratio of %.2f, %i concentrations (%s)",
          type, drug_pairs[idp, 1], drug_pairs[idp, 2], 10**as.numeric(names(conc_ratio)),
            length(conc_2), condition_name))
      condition <- c(as.list(drug_pairs[idp, ]), list(conc_ratio = 10**as.numeric(names(conc_ratio))))
    } else if (n_conc_pairs == length(conc_1) * length(conc_2) & length(conc_2) >= 4) {
      type <- "matrix"
      print(sprintf("Found %s combination with %s and %s: %i x %i concentrations (%s)",
          type, drug_pairs[idp, 1], drug_pairs[idp, 2], length(conc_1), length(conc_2), condition_name))
      condition <- c(as.list(drug_pairs[idp, ]), list(
            Concentration = as.numeric(names(conc_1)), Concentration_2 = as.numeric(names(conc_2))))
    } else if (length(conc_2) < 4) {
      type <- "fixed"
      print(sprintf("Found %s combination of %s with %s at %.3g uM (%s)",
          type, drug_pairs[idp, 1], drug_pairs[idp, 2], as.numeric(names(conc_2)), condition_name))
      condition <- c(as.list(drug_pairs[idp, ]), list(Concentration_2 = as.numeric(names(conc_2))))
    } else {
      type <- "other"
      print(sprintf("Found %s combination with %s and %s: %i concentration pairs (%s)",
        type, drug_pairs[idp, 1], drug_pairs[idp, 2], n_conc_pairs, condition_name))
      condition <- c(as.list(drug_pairs[idp, ]), list(
            Concentration = as.numeric(names(conc_1)), Concentration_2 = as.numeric(names(conc_2))))
    }

    pair_list[[idp]] <- list(condition = condition,
                          rows = rownames(r_data)[row_idx],
                          type = type,
                          name = sprintf("%s x %s (%s)", drug_pairs[idp, 1], drug_pairs[idp, 2], 
                              condition_name))
  }

  S4Vectors::metadata(se)$drug_combinations <- pair_list
  return(se)
}




# Perform the calculations on all references.
# Then, take the mean to report the final reference normalized value.
.calculate_references <- function(ref_df) {
  RV_vec <- ref_df$RefReadout / ref_df$UntrtReadout
  GR_vec <- calculate_GR_value(rel_viability = RV_vec, 
                               corrected_readout = ref_df$RefReadout, 
                               day0_readout = ref_df$Day0Readout, 
                               untrt_readout = ref_df$UntrtReadout, 
                               ndigit_rounding = ndigit_rounding, 
                               duration = duration, 
                               ref_div_time = ref_div_time, 
                               cl_name = cl_name)

  ref_rel_viability[i, j] <- round(mean(RV_vec, na.rm = TRUE), ndigit_rounding)
  ref_GR_value[i, j] <- round(mean(GR_vec, na.rm = TRUE), ndigit_rounding)
  div_time[i, j] <- round(duration / log2(mean(ref_df$UntrtReadout / ref_df$Day0Readout,
                                               na.rm = TRUE)), ndigit_rounding)
}
