.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$",
                          gDRutils::get_env_identifiers("drug_name"),
                          gDRutils::get_env_identifiers("drug_name"))

.untreated_tag_patterns <- vapply(gDRutils::get_env_identifiers("untreated_tag"), sprintf, fmt = "^%s$", character(1))
.untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")
