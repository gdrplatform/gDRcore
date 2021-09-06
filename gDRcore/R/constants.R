.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$",
                          gDRutils::get_identifier("drugname"),
                          gDRutils::get_identifier("drugname"))
