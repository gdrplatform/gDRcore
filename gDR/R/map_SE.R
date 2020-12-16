#' map_SE
#'
#' Perfmorm mapping for normalization
#'
#' @param normSE a SummarizedExperiment with normalization assaay
#' @param ctrlSE a SummarizedExperiment object with information for controls
#' @param row_endpoint_value_filter an array with key values for end points
#' @param Keys a list of all identified keys
#' @param T0 a logical indicating if the mapping should be performer for Time=0 (FALSE by default)
#'
#' @return a list of mapping
#' @export

map_SE <- function(normSE, ctrlSE, row_endpoint_value_filter, Keys, T0 = FALSE){
    # Assertions:
    checkmate::assert_class(normSE, "SummarizedExperiment")
    checkmate::assert_class(ctrlSE, "SummarizedExperiment")
    checkmate::assert_array(row_endpoint_value_filter)
    checkmate::assert_list(Keys)
    checkmate::assert_logical(T0)

    mappingFactor <- ifelse(T0, "Day0", "untrt_Endpoint")

    keyValuesList <- list(key_values = row_endpoint_value_filter)
    if(T0){
      keyValuesList <- list(T0 = SummarizedExperiment::rowData(ctrlSE)[, gDRutils::get_identifier("duration")] == 0)
    }

    matchFactor <- ifelse(T0, "T0", gDRutils::get_identifier("duration"))

    lapply(rownames(normSE), function(x) {
    # define matix with matching metadata
      ctrl_metadata_idx = intersect(Keys[[mappingFactor]],
                                    names(SummarizedExperiment::rowData(ctrlSE)))
      names(ctrl_metadata_idx) = ctrl_metadata_idx
      match_mx <-
        IRanges::LogicalList(c(
          lapply(ctrl_metadata_idx, function(y) # matching the metadata
            SummarizedExperiment::rowData(ctrlSE)[,y] ==
              SummarizedExperiment::rowData(normSE)[x,y]),
          c(keyValuesList, list(conc = apply(cbind(array(0, nrow(ctrlSE)), # padding to avoid empty df;
                                  SummarizedExperiment::rowData(ctrlSE)[, agrep("Concentration",
                                                                                colnames(SummarizedExperiment::rowData(ctrlSE))), drop = FALSE]), 1,
                            function(x)
                              all(x == 0))

          ))))
      match_idx <- which(apply(as.matrix(match_mx), 2, all)) # test matching conditions
      if (length(match_idx) == 0) {
        # if not exact match, try to find best match (as many metadata fields as possible)
        futile.logger::flog.warn("Missing untreated controls %s for: %s",
            ifelse(T0, '(T=0)', '(endpoint)'), x)
        idx <-
          apply(as.matrix(match_mx), 2, function(y)
            sum(y, na.rm = TRUE)) *
          match_mx[[matchFactor]]
        if (any(idx > 0)) {
          match_idx <- which.max(idx)
          futile.logger::flog.warn("Found partial match:",
                                   rownames(ctrlSE)[match_idx])
        } else { # failed to find any potential match
          futile.logger::flog.warn("No partial match found")
        }
      }
      return(rownames(ctrlSE)[match_idx])
    })
}
