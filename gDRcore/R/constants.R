.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$", 
                          gDRutils::get_identifier("drugname"), 
                          gDRutils::get_identifier("drugname"))

.untreated_tag_patterns <- vapply(gDRutils::get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
.untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")

DICTIONARY <- list(
  Untreated = "untreated",
  Vehicle = "vehicle",
  vehcle = "vehicle"
)
