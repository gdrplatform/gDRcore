#' average_SE
#'
#' Avereage normalized SummarizedExperiment of DR data
#'
#' @param normSE a SummarizedExperiment with normalized DR data
#' @param TrtKeys a vector of keys used for averaging (NULL by default)
#'
#' @return a SummarizedExperiment with additional assay with averaged DR data
#' @export
#'

average_SE <- function(normSE, TrtKeys = NULL, include_masked = F) {

  # Assertions:
  checkmate::assert_class(normSE, "SummarizedExperiment")
  checkmate::assert_vector(TrtKeys, null.ok = TRUE)


    avgSE <- normSE
    if (is.null(TrtKeys)) {
        if ("Keys" %in% names(metadata(normSE))) {
          TrtKeys <- metadata(normSE)$Keys$Trt
          TrtKeys <- setdiff(TrtKeys, metadata(normSE)$Keys$discard_keys)
        } else {
          TrtKeys <- identify_keys(normSE)$Trt
        }
    }
    metadata(normSE)$Keys$Trt <- TrtKeys

    SummarizedExperiment::assay(avgSE, "Averaged") <- SummarizedExperiment::assay(avgSE, "Normalized")
    avgSE <- aapply(avgSE, function(x) {
        # bypass 'masked' filter
        x$masked = x$masked & !include_masked

        subKeys <- intersect(TrtKeys, colnames(x))
        if (sum(!x$masked) >= 1) {
            df_av <- aggregate(x[ !x$masked ,
                                  c("GRvalue", "RelativeViability","CorrectedReadout")],
                            by = as.list(x[ !x$masked , subKeys, drop = FALSE]),
                            FUN = function(y) mean(y, na.rm = TRUE))
            df_std <- aggregate(x[!x$masked, c("GRvalue", "RelativeViability")],
                                by = as.list(x[ !x$masked, subKeys, drop = FALSE]),
                                FUN = function(x) sd(x, na.rm = TRUE))
            colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")] =
                paste0("std_",
                    colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")])
            return( merge(df_av, df_std, by = subKeys) )
        } else { # case: (nrow(x) == 0 || all(x$masked))
            df_ = as.data.frame(matrix(0,0,length(subKeys)+5))
            colnames(df_) = c(subKeys,
                  c("GRvalue", "RelativeViability","CorrectedReadout"),
                  paste0("std_", c("GRvalue", "RelativeViability")))
            return(df_)
        } 
    }, "Averaged")

    SummarizedExperiment::assay(avgSE, "Avg_Controls") <- SummarizedExperiment::assay(avgSE, "Controls")
    avgSE <- aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys <- intersect(TrtKeys, colnames(x))
            df_av <- DataFrame(lapply(x[, c("Day0Readout", "UntrtReadout",
                    "RefGRvalue", "RefRelativeViability",
                    "RefReadout", "DivisionTime")], FUN = function(y) mean(y, na.rm = TRUE)))
            return( df_av )
        } else return(x)
    }, "Avg_Controls")

    return(avgSE)
}
<<<<<<< HEAD
=======

#' metrics_SE
#'
#' Calculate metrics for DR data
#'
#' @param avgSE a SummarizedExperiment with averaged and normalized assays
#' @param studyConcThresh a numeric with study concentration threshold (4 by default)
#'
#' @return a SummarizedExperiment with additional assay with metrics
#' @export
#'

metrics_SE = function(avgSE, studyConcThresh = 4) {

    # Assertions:
    checkmate::assert_class(avgSE, "SummarizedExperiment")
    checkmate::assert_number(studyConcThresh)

    stopifnot(is.numeric(studyConcThresh))
    # this is not used as we enforce the same conditions as the input SE; not collapsing allowed
    # if (is.null(DoseRespKeys)) {
    #     if ("Keys" %in% names(metadata(avgSE))) DoseResp = metadata(avgSE)$Keys$DoseResp
    #     else DoseRespKeys = identify_keys(avgSE)$DoseResp
    # } else {
    #     metadata(avgSE)$Keys$DoseResp = DoseRespKeys
    # }

    metricsSE <- avgSE
    SummarizedExperiment::assay(metricsSE, "Metrics") <- SummarizedExperiment::assay(metricsSE, "Averaged")

    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in a foor loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    mSE_m <- SummarizedExperiment::assay(metricsSE, "Metrics")
    a_SE = SummarizedExperiment::assay(metricsSE, "Averaged")
    aCtrl_SE = SummarizedExperiment::assay(metricsSE, "Avg_Controls")
    for (i in rownames(metricsSE)) {
        for (j in colnames(metricsSE)) {
            df_ <- a_SE[[i, j]]
            if (!is.null(df_) && all(dim(df_) > 0)) { # studyConcThresh is embeded in RVGRfits
                mSE_m[[i, j]] <- DataFrame(gDRutils::RVGRfits(df_,
                    e_0 = aCtrl_SE[[i, j]]$RefRelativeViability,
                    GR_0 = aCtrl_SE[[i, j]]$RefGRvalue,
                    n_point_cutoff = studyConcThresh))
            } else {
                out <- DataFrame(matrix(NA, 0, length(gDRutils::get_header("response_metrics"))+2))
                colnames(out) <- c(gDRutils::get_header("response_metrics"), "maxlog10Concentration", "N_conc")
                mSE_m[[i, j]] <- out
            }
        }
    }
    SummarizedExperiment::assay(metricsSE, "Metrics") <- mSE_m
    return(metricsSE)
}

#' identify_keys
#'
#' Identify keys in the DR data represented by dataframe or SummarizedExperiment or MultiAssayExperiment objects
#'
#' @param df_se_mae a dataframe or SummarizedExperiment or MultiassayExperiment with keys
#'
#' @return a list of keys
#' @export
#'

identify_keys <- function(df_se_mae) {

  # Assertions:
  stopifnot(inherits(df_se_mae, c("data.frame", "MultiAssayExperiment", "SummarizedExperiment")))


    if (any(class(df_se_mae) %in% c("MultiAssayExperiment", "SummarizedExperiment"))) {
        if ("MultiAssayExperiment" %in% class(df_se_mae)) {
            # if MAE, convert to SE based on the treated SE (could be optimized)
            df_se_mae <- df_se_mae[["treated"]]
            se_untrt <-  df_se_mae[["untreated"]]
        } else se_untrt <- NULL
        all_keys <- unique(c(
            colnames(SummarizedExperiment::rowData(df_se_mae)),
            colnames(SummarizedExperiment::colData(df_se_mae)),
            unlist(lapply(SummarizedExperiment::assay(df_se_mae), colnames))))
    } else { # case of a data frame
        all_keys <- colnames(df_se_mae)
    }

    keys <- list(Trt = setdiff(all_keys, "Barcode"),
            DoseResp = setdiff(all_keys,  "Barcode"),
            ref_Endpoint = setdiff(all_keys, c("Concentration",
                                            gDRutils::get_identifier("drug"),
                                            gDRutils::get_identifier("drugname"))),
            untrt_Endpoint = all_keys[ c(-agrep("Concentration", all_keys),
                                            -agrep(gDRutils::get_identifier("drug"), all_keys),
                                            -agrep(gDRutils::get_identifier("drugname"), all_keys))])
    keys[["Day0"]] <- setdiff(keys[["untrt_Endpoint"]], gDRutils::get_identifier("duration"))
    keys <- lapply(keys, function(x) setdiff(x, c(gDRutils::get_header("raw_data"),
        gDRutils::get_header("normalized_results"), "Template", gDRutils::get_identifier("WellPosition"), gDRutils::get_header("averaged_results"),
            gDRutils::get_header("metrics_results"), "ReferenceDivisionTime"
    )))
    keys <- lapply(keys, sort)

    # check if all values of a key is NA
    for (k in keys[["untrt_Endpoint"]]) {

        if ("SummarizedExperiment" %in% class(df_se_mae)) {
            # check the metadata fields for NA
            if (k %in% colnames(SummarizedExperiment::rowData(df_se_mae))) df_ <- SummarizedExperiment::rowData(df_se_mae)
            else if (k %in% colnames(SummarizedExperiment::colData(df_se_mae))) df_ <- SummarizedExperiment::colData(df_se_mae)
            else next # not a metadata

            if (all(is.na(df_[,k]))) keys <- lapply(keys, function(x) setdiff(x, k))

            if (!is.null(se_untrt) && k %in% colnames(SummarizedExperiment::rowData(se_untrt))) {
                df_ <- SummarizedExperiment::rowData(se_untrt)
                if (all(is.na(df_[df_[,gDRutils::get_identifier("duration")] == 0, k]))) {
                    keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
                }
            }
        } else { # case of a data frame
            if (all(is.na(df_se_mae[, k]))) {
                keys <- lapply(keys, function(x) setdiff(x, k))
            }
            if (all(is.na(df_se_mae[df_se_mae[,gDRutils::get_identifier("duration")] == 0, k]))) {
                keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
            }
        }
    }
  return(keys)
}

#' cleanup_metadata
#'
#' Cleanup a dataframe with metadata
#'
#' @param df_metadata a dataframe with metadata
#'
#' @return a dataframe with cleaned metadata
#'
#' @export

cleanup_metadata <- function(df_metadata) {

  # Assertions:
  stopifnot(inherits(df_metadata, "data.frame"))

  # clean up numberic fields
  df_metadata[, gDRutils::get_identifier("duration")] <-
    round(as.numeric(df_metadata[, gDRutils::get_identifier("duration")]), 6)
  # identify potential numeric fields and replace NA by 0 - convert strings in factors
  for (c in setdiff(1:dim(df_metadata)[2], c(
    agrep(gDRutils::get_identifier("drug"), colnames(df_metadata)),
    agrep("Concentration", colnames(df_metadata)),
    grep(paste(
      c(
        gDRutils::get_identifier("cellline"),
        gDRutils::get_header("manifest"),
        gDRutils::get_identifier("WellPosition")
      ),
      collapse = "|"
    ), colnames(df_metadata))
  ))) {
    vals <- unique(df_metadata[, c])

    if (is.character(vals)) {
      num_vals <- as.numeric(vals)
      if (sum(is.na(num_vals)) > 2 || all(is.na(num_vals))) {
        df_metadata[, c] <- factor(df_metadata[, c])
        futile.logger::flog.warn("Metadata field %s converted to factors",
                colnames(df_metadata)[c])
      } else {
        is.na(df_metadata[, c]) <- 0
        df_metadata[, c] <- as.numeric(df_metadata[, c])
        futile.logger::flog.warn("Metadata field %s converted to numeric values",
                colnames(df_metadata)[c])
      }
    }
  }
    # TODO: specific to GNE database --> need to be replaced by a function
    df_metadata <- add_CellLine_annotation(df_metadata)

    # check that Gnumber_* are in the format 'G####' and add common name (or Vehicle or Untreated)

    for (i in agrep(gDRutils::get_identifier("drug"), colnames(df_metadata))) { # correct case issues
        for (w in gDRutils::get_identifier("untreated_tag")) {
            df_metadata[grep(w, df_metadata[,i], ignore.case = T),i] <- w
        }
    }
    # -----------------------

    df_metadata <- add_Drug_annotation(df_metadata)

    # clean up concentration fields
    for (i in agrep("Concentration", colnames(df_metadata))) {
        trt_n <- ifelse(regexpr("_\\d", colnames(df_metadata)[i]) > 0,
                            substr(colnames(df_metadata)[i], 15, 20), 1)
        DrugID_col <- ifelse(trt_n == 1, gDRutils::get_identifier("drug"), paste0(gDRutils::get_identifier("drug"), "_", trt_n))
        df_metadata[df_metadata[,DrugID_col] %in% gDRutils::get_identifier("untreated_tag"), i] <- 0 # set all untreated to 0

        DrugID_0 <- setdiff(unique(df_metadata[ df_metadata[,i] == 0, DrugID_col]), gDRutils::get_identifier("untreated_tag"))
        DrugID_0 <- DrugID_0[!is.na(DrugID_0)]
        if (length(DrugID_0) > 0) {
          futile.logger::flog.warn("Some concentration for %s are 0: %s",
                                   DrugID_col,
                                   paste(DrugID_0, collapse = " ; "))

        }
        df_metadata[,i] <- 10 ** round(log10(as.numeric(df_metadata[, i])), 6)
        # df_metadata[,i] <- round(as.numeric(df_metadata[, i]), 10) # avoid mismatch due to string truncation
    }
  return(df_metadata)
}



#' Order_result_df
#'
#' Order a dataframe with results
#'
#' @param df_ a dataframe with results
#'
#' @return a ordered dataframe with results
#' @export

Order_result_df <- function (df_) {

  # Assertions:
  stopifnot(inherits(df_, "data.frame"))

  cols <- c(gDRutils::get_header("ordered_1"),
            setdiff(colnames(df_),
                    c(
                      gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
                    )),
            gDRutils::get_header("ordered_2"))
  cols <- intersect(cols, colnames(df_))

  row_order_col <-
    intersect(
      c(
        gDRutils::get_header("add_clid")[1],
        gDRutils::get_identifier("duration"),
        gDRutils::get_identifier("drugname"),
        "Concentration",
        paste0(c(
          paste0(gDRutils::get_identifier("drugname"), "_"), "Concentration_"
        ),
        sort(c(2:10, 2:10))),
        setdiff(colnames(df_), c(
          gDRutils::get_header("ordered_1"), gDRutils::get_header("ordered_2")
        ))
      ),
      cols
    )

  df_ <- df_[do.call(order, df_[, row_order_col]), cols]

  return(df_)
}


#' Add codrug group
#'
#' @param SE 
#'
#' @return
#' @export
#'
add_codrug_group = function(SE) {

  r_data = SummarizedExperiment::rowData(SE)
  if (!(paste0(gDRutils::get_identifier()$drugname, '_2') %in% colnames(r_data))) return(SE)

  # find the pairs of drugs with relevant metadata
  drug_ids = paste0(gDRutils::get_identifier()$drugname, c('', '_2'))
  other_metadata = c(paste0(gDRutils::get_identifier()$drug, c('', '_2')),
            setdiff(colnames(r_data), c('Concentration_2', drug_ids,
                paste0(gDRutils::get_identifier()$drug, c('', '_2')))))
  drug_pairs = unique(r_data[, c(drug_ids, other_metadata)])
  drug_pairs = drug_pairs[ !(drug_pairs[,drug_ids[2]] %in% gDRutils::get_identifier('untreated_tag')),]

  pair_list = vector('list', nrow(drug_pairs))
  # loop through the pairs to assess the number of individual concentration pairs
  for (idp in 1:nrow(drug_pairs)) {
    row_idx = r_data[,drug_ids[1]] %in% unlist(drug_pairs[idp, drug_ids]) &
            r_data[,drug_ids[2]] %in% c(unlist(drug_pairs[idp, drug_ids]),
                gDRutils::get_identifier('untreated_tag')) &
            apply(as.matrix(
                IRanges::LogicalList(c(
                  lapply(setdiff(other_metadata,
                      paste0(gDRutils::get_identifier()$drug, c('', '_2'))),
                    function(y) # matching the metadata
                    r_data[,y] == drug_pairs[idp,y])
                  ))), 2, all)

    # reverse engineer the type of combination experiment
    flat_data = gDRutils::assay_to_df(SE[row_idx, ], 'Averaged')
    flat_data = flat_data[flat_data$Concentration_2 > 0,]
    conc_1 = table(flat_data$Concentration)
    conc_2 = table(flat_data$Concentration_2)
    n_conc_pairs = nrow(unique(flat_data[,c('Concentration', 'Concentration_2')]))
    conc_ratio = table(round(log10(flat_data$Concentration / flat_data$Concentration_2),2))
    conc_ratio = conc_ratio[!names(conc_ratio) %in% c('Inf', '-Inf')]

    condition = paste(paste(other_metadata, unlist(drug_pairs[idp,other_metadata]), sep = '='),
                  collapse=' ')
    if (length(conc_ratio) <= 2) {
      type = 'co-dilution'
      print(sprintf('Found %s combination with %s and %s: ratio of %.2f, %i concentrations (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], 10**as.numeric(names(conc_ratio)),
            length(conc_2), condition))
    } else if (n_conc_pairs == length(conc_1)*length(conc_2) & length(conc_2) >= 4) {
      type = 'matrix'
      print(sprintf('Found %s combination with %s and %s: %i x %i concentrations (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], length(conc_1), length(conc_2), condition))
    } else if (length(conc_2)<4) {
      type = 'fixed'
      print(sprintf('Found %s combination of %s with %s at %.3g uM (%s)',
          type, drug_pairs[idp,1], drug_pairs[idp,2], as.numeric(names(conc_2)), condition))
    } else {
      type = 'other'
      print(sprintf('Found %s combination with %s and %s: %i concentration pairs (%s)',
        type, drug_pairs[idp,1], drug_pairs[idp,2], n_conc_pairs, condition))
    }

    pair_list[[idp]] = list(condition = unlist(drug_pairs[idp,]),
                          rows = rownames(r_data)[row_idx],
                          type = type)
  }

  metadata(SE)$drug_combinations = pair_list
  return(SE)
}
