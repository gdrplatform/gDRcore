

# openxlsx skip the first emprty rows and cannot be overridden --> use readxl
#' @import readxl
#' @import readr
#' @import stringr


#' @export
get_identifier <- function(x = NULL) {
  identifiersList <- list(
    duration = "Duration",
    
    cellline = "clid",
    DB_cell = "clid",
    
    drug = "Gnumber",
    DB_drug = "drug",
    drugname = "DrugName",
    # corresponds to the fieLd  "gcsi_drug_name" from gCellGenomics::getDrugs()
    
    untreated_tag = c("untreated", "vehicle"),
    # flag to identify control treatments
    
    WellPosition = c("WellRow", "WellColumn")
  )
  if (!is.null(x) &&
      x %in% names(identifiersList))
    return(identifiersList[[x]])
  else
    return(identifiersList)
}

#######-------------------------------------------------------
# these should not be changed and are protected field names
#' @export
get_header <- function(x = NULL) {
  headersList <- list(
    manifest = c("Barcode", "Template", get_identifier("duration")),
    raw_data = c(
      "ReadoutValue",
      "BackgroundValue",
      "UntrtReadout",
      "Day0Readout"
    ),
    normalized_results = c(
      "CorrectedReadout",
      "GRvalue",
      "RelativeViability",
      "DivisionTime"
    ),
    averaged_results = c("std_GRvalue", "std_RelativeViability"),
    metrics_results = c(
      "maxlog10Concentration",
      "N_conc",
      "mean_viability",
      "ic50",
      "e_max",
      "ec50",
      "e_inf",
      "e_0",
      "h_ic",
      "ic_r2",
      "flat_fit_ic",
      "GR_AOC",
      "GR50",
      "GRmax",
      "GEC50",
      "GRinf",
      "GR_0",
      "h_GR",
      "GR_r2",
      "flat_fit_GR"
    ),
    add_clid = c("CellLineName", "Tissue", "ReferenceDivisionTime")
    # corresponds to the fieLd  "celllinename", "primarytissue", "doublingtime" from gneDB CLIDs
  )
  headersList[["controlled"]] <- c(
    get_identifier("cellline"),
    headersList[["manifest"]],
    get_identifier("drug"),
    "Concentration",
    paste0(get_identifier("drug"), "_", 2:10),
    paste0("Concentration_", 2:10)
  )
  headersList[["reserved"]] <-
    c(
      headersList[["add_clid"]],
      get_identifier("drugname"),
      paste0(get_identifier("drugname"), "_", 2:10),
      headersList[["raw_data"]],
      headersList[["normalized_results"]],
      headersList[["averaged_results"]],
      headersList[["metrics_results"]],
      "WellRow",
      "WellColumn"
    )
  
  headersList[["ordered_1"]] <- c(
    headersList[["add_clid"]][1:2],
    get_identifier("duration"),
    get_identifier("drugname"),
    "Concentration",
    paste0(c(
      paste0(get_identifier("drugname"), "_"), "Concentration_"
    ),
    sort(c(2:10, 2:10)))
  )
  headersList[["ordered_2"]] <- c(
    headersList[["normalized_results"]],
    headersList[["averaged_results"]],
    headersList[["metrics_results"]],
    headersList[["raw_data"]],
    headersList[["add_clid"]][-2:-1],
    get_identifier("cellline"),
    get_identifier("drug"),
    paste0(get_identifier("drug"), "_", 2:10),
    headersList[["manifest"]],
    "WellRow",
    "WellColumn"
  )
  
  if (!is.null(x) &&
      x %in% names(headersList))
    return(headersList[[x]])
  else
    return(headersList)
}

#' Load data
#' 
#' This functions loads and checks the data file(s)
#' 
#' @param manifest_file character, file path(s) to manifest(s)
#' @param df_template_files data.frame, with datapaths and names of results file(s)
#' or character with file path of templates file(s)
#' @param results_file  data.frame, with datapaths and names of results file(s)
#' or character with file path of results file(s)
#' @param log_str character, file path to logs
#' @param instrument character
#' @export
load_data <-
  function(manifest_file,
           df_template_files,
           results_file,
           log_str,
           instrument = "EnVision") {
    log_str <- c(log_str, "", "load_merge_data")
    
    if (is.data.frame(df_template_files)) {
      # for the shiny app
      template_file <- df_template_files$datapath
      template_filename <- df_template_files$name
    } else {
      template_filename <- df_template_files
    }
    
    manifest <- load_manifest(manifest_file, log_str)
    treatments <- load_templates(df_template_files, log_str)
    data <- load_results(results_file, log_str, instrument)
    
    # check the all template files are available
    if (!all(unique(manifest$Template[manifest$Barcode %in% data$Barcode])
             %in% basename(template_filename))) {
      ErrorMsg <- paste("Some template files are missing:",
                        paste(setdiff(
                          unique(manifest$Template[manifest$Barcode %in% data$Barcode]),
                          basename(template_filename)
                        ), collapse = " ; "))
      stop(ErrorMsg)
    }
    return(list(
      manifest = manifest,
      treatments = treatments,
      data = data
    ))
  }



#' Load manifest
#' 
#' This functions loads and checks the manifest file(s)
#' 
#' @param manifest_file character, file path(s) to manifest(s)
#' @param log_str character, file path to logs
#' @export
load_manifest <- function (manifest_file, log_str) {
  # manifest_file is a string or a vector of strings
  
  log_str <- c(log_str, "", "load_manifest")
  available_formats <- c("text/tsv",
                         "text/tab-separated-values",
                         "xlsx", "xls", "tsv")
  
  # read files
  manifest_data <- lapply(manifest_file, function(x) {
    manifest_ext <- tools::file_ext(x)
    if (manifest_ext %in% c("xlsx", "xls")) {
      df <- tryCatch({
        readxl::read_excel(x, col_names = TRUE)
      }, error = function(e) {
        stop(sprintf("Error reading the Manifest file. Please see the logs:\n%s", e))
      }) 
    } else if (manifest_ext %in% c("text/tsv",
                                   "text/tab-separated-values",
                                   "tsv")) {
      df <- tryCatch({
        readr::read_tsv(x, col_names = TRUE, skip_empty_rows = TRUE)
      }, error = function(e) {
        stop(sprintf("Error reading the Manifest file. Please see the logs:\n%s", e))
      })
    } else {
        stop(sprintf(
          "%s file format is not supported.
          Please convert your file to one of the follwoing: %s",
          manifest_ext,
          stringi::stri_flatten(available_formats, collapse = ", ")
        )
      )
    }
  })
  
  # replace Time by Duration for backwards compatibility
  manifest_data <- lapply(manifest_data, function(x) {
    if ("Time" %in% colnames(x)) {
      colnames(x)[colnames(x) == "Time"] <-
        get_identifier("duration")
    }
    return(x)
  })
  
  # check default headers are in each df
  dump <- sapply(1:length(manifest_file),
                 function(i)
                   check_metadata_names(
                     colnames(manifest_data[[i]]),
                     log_str,
                     df_name = manifest_file[[i]],
                     df_type = "manifest"
                   ))
  
  cat_manifest_data <- dplyr::bind_rows(manifest_data)
  colnames(cat_manifest_data) <-
    check_metadata_names(colnames(cat_manifest_data),
                         log_str, "manifest")
  
  # check that barcodes are unique
  if (dim(cat_manifest_data)[1] != length(unique(cat_manifest_data$Barcode)))
    stop("Barcodes in Manifest must be unique!")
  
  cat_manifest_data$Template <- basename(cat_manifest_data$Template)
  
  print("Manifest loaded successfully")
  return(cat_manifest_data)
}


#' Load templates from
#' 
#' This functions loads and checks the template file(s)
#' 
#' @param df_template_files data.frame, with datapaths and names of results file(s)
#' or character with file path of templates file(s)
#' @param log_str character, file path to logs
#' @export
load_templates <- function (df_template_files, log_str) {
  # template_file is a string or a vector of strings
  log_str <- c(log_str, "", "load_templates")
  
  if (is.data.frame(df_template_files)) {
    # for the shiny app
    template_file <- df_template_files$datapath
    template_filename <- df_template_files$name
  } else {
    template_file <- df_template_files
    template_filename <- basename(template_file)
  }
  
  all_templates <- data.frame()
  if (any(grepl("\\.xlsx?$", template_filename))) {
    idx <- grepl("\\.xlsx?$", template_filename)
    all_templates_1 <- load_templates_xlsx(template_file[idx], log_str)
    all_templates <- rbind(all_templates, all_templates_1)
  }
  if (any(grepl("\\.[ct]sv$", template_filename))) {
    idx <- grepl("\\.[ct]sv$", template_filename)
    print(paste("Reading", template_filename[idx], "with load_templates_tsv"))
    all_templates_2 <- load_templates_tsv(template_file[idx], log_str)
    all_templates <- rbind(all_templates, all_templates_2)
  }
  
  return(all_templates)
  
}

#' Load results
#' 
#' This functions loads and checks the results file(s)
#' 
#' @param df_results_files  data.frame, with datapaths and names of results file(s)
#' or character with file path of results file(s)
#' @param log_str character, file path to logs
#' @param intrument character
#' @export
#' @export
load_results <-
  function(df_results_files, log_str, instrument = "EnVision") {
    if (is.data.frame(df_results_files)) {
      # for the shiny app
      results_file <- df_results_files$datapath
      results_filename <- df_results_files$name
    } else {
      results_file <- df_results_files
      results_filename <- basename(results_file)
    }
    stopifnot(sapply(results_file, file.exists))
    
    if (instrument == "EnVision") {
      all_results <-
        load_results_EnVision(results_file, log_str)
    } else if (instrument == "long_tsv") {
      all_results <-
        load_results_tsv(results_file, log_str)
    }
    return(all_results)
  }



# individual functions

#' Load templates from tsv
#' 
#' This functions loads and checks the template file(s)
#' 
#' @param template_file character, file path(s) to template(s)
#' @param log_str character, file path to logs
load_templates_tsv <-
  function(template_file,
           log_str) {
    template_filename <- basename(template_file)
    
    # read columns in files
    templates <- lapply(template_file, function(x)
      readr::read_tsv(x, col_names = TRUE, skip_empty_rows = TRUE))
    names(templates) <- template_filename
    # check WellRow/WellColumn is present in each df
    dump <- sapply(1:length(template_file),
                   function(i)
                     if (!(all(
                       get_identifier("WellPosition") %in% colnames(templates[[i]])
                     ))) {
                       print(paste(
                         template_filename[[i]],
                         "missing",
                         get_identifier("WellPosition"),
                         "as header"
                       ))
                     })
    # check drug_identifier is present in each df
    dump <- sapply(1:length(template_file),
                   function(i)
                     check_metadata_names(
                       setdiff(colnames(templates[[i]]), get_identifier("WellPosition")),
                       log_str,
                       df_name = template_filename[[i]],
                       df_type = "template"
                     ))
    
    metadata_fields <- NULL
    all_templates <- data.frame()
    for (iF in 1:length(template_file)) {
      print(paste("Loading", template_filename[iF]))
      # first check that the sheet names are ok
      # identify drug_identifier sheet (case insensitive)
      Gnumber_idx <- grep(paste0(get_identifier("drug"), "$"),
                          colnames(templates[[iF]]),
                          ignore.case = TRUE)
      Conc_idx <-
        grepl("Concentration", colnames(templates[[iF]]), ignore.case = TRUE)
      # case of untreated plate
      if (sum(Conc_idx) == 0) {
        if (length(Gnumber_idx) == 0) {
          ErrorMsg <- sprintf(
            "In untreated template file %s, sheet name must be %",
            template_file[[iF]],
            get_identifier("drug")
          )
          stop(ErrorMsg)
        }
        df <- templates[[iF]][, get_identifier("drug")]
        if (!(all(toupper(df)[!is.na(df)]) %in% toupper(get_identifier("untreated_tag")))) {
          ErrorMsg <- sprintf(
            "In untreated template file %s, entries must be %s",
            template_file[[iF]],
            paste(get_identifier("untreated_tag"), collapse = " or ")
          )
          stop(ErrorMsg)
        }
      } else {
        # normal case
        check_metadata_names(colnames(templates[[iF]]),
                             log_str,
                             df_name = template_filename[iF],
                             df_type = "template_treatment")
      }
      
      df_template <- templates[[iF]]
      for (iS in colnames(df_template)) {
        # check if metadata field already exist and correct capitalization if needed
        if (!(iS %in% metadata_fields)) {
          if (!is.null(metadata_fields) &&
              toupper(iS) %in% toupper(metadata_fields)) {
            oldiS <- iS
            iS <-
              metadata_fields[toupper(iS) == toupper(metadata_fields)]
            print(paste(oldiS, "corrected to match case with ", iS))
            colnames(df_template)[colnames(df_template) == oldiS] <-
              iS
          } else {
            metadata_fields <- c(metadata_fields, iS)
          }
        }
      }
      df_template$Template <- template_filename[iF]
      colnames(df_template) <-
        check_metadata_names(colnames(df_template), log_str,
                             df_name = template_filename[iF])
      all_templates <- dplyr::bind_rows(all_templates, df_template)
      
    }
    print("Templates loaded successfully!")
    return(all_templates)
  }

#' Load templates from xlsx
#' 
#' This functions loads and checks the template file(s)
#' 
#' @param template_file character, file path(s) to template(s)
#' @param log_str character, file path to logs
load_templates_xlsx <-
  function(template_file,
           log_str) {
    template_filename <- basename(template_file)
    # read sheets in files
    template_sheets <- lapply(template_file, readxl::excel_sheets)
    # check drug_identifier is present in each df
    dump <- sapply(1:length(template_file),
                   function(i)
                     check_metadata_names(
                       template_sheets[[i]],
                       log_str,
                       df_name = template_file[[i]],
                       df_type = "template"
                     ))
    
    metadata_fields <- NULL
    all_templates <- data.frame()
    for (iF in 1:length(template_file)) {
      print(paste("Loading", template_filename[iF]))
      # first check that the sheet names are ok
      # identify drug_identifier sheet (case insensitive)
      Gnumber_idx <- grep(paste0(get_identifier("drug"), "$"),
                          template_sheets[[iF]],
                          ignore.case = TRUE)
      Conc_idx <- grepl("Concentration", template_sheets[[iF]], ignore.case = TRUE)
      # case of untreated plate
      if (sum(Conc_idx) == 0) {
        if (length(Gnumber_idx) == 0) {
          ErrorMsg <- sprintf(
            "In untreated template file %s, sheet name must be %",
            template_file[[iF]],
            get_identifier("drug")
          )
          stop(ErrorMsg)
        }
        tryCatch({
          df <-
            readxl::read_excel(
              template_file[[iF]],
              sheet = Gnumber_idx,
              col_names = paste0("x", 1:48),
              range = "A1:AV32"
            )
        }, error = function(e) {
          stop(sprintf("Error loading template. See logs: %s", e))
        })
        if (!(all(toupper(unlist(df)[!is.na(unlist(df))]) %in% 
                  toupper(get_identifier("untreated_tag"))))) {
          
          ErrorMsg <- sprintf(
            "In untreated template file %s, entries must be %s",
            template_file[[iF]],
            paste(get_identifier("untreated_tag"), collapse = " or ")
          )
          stop(ErrorMsg)
        }
      } else {
        # normal case
        dump <- check_metadata_names(template_sheets[[iF]],
                                     log_str,
                                     df_name = template_filename[iF],
                                     df_type = "template_treatment")
      }
      # read the different sheets and check for plate size
      # enforce range to avoid skipping empty rows at the beginning
      tryCatch({
        df <-
          readxl::read_excel(
            template_file[[iF]],
            sheet = template_sheets[[iF]][Gnumber_idx],
            col_names = paste0("x", 1:48),
            range = "A1:AV32",
            col_types = "text"
          )
      }, error = function(e) {
        stop(sprintf("Error loading template. See logs: %s", e))
      })
      # get the plate size
      n_row <-
        2 ** ceiling(log2(max(which(
          apply(!is.na(df), 1, any)
        ))))
      n_col <-
        1.5 * 2 ** ceiling(log2(max(which(
          apply(!is.na(df), 2, any)
        )) / 1.5))
      n_row <- max(n_row, n_col / 1.5)
      n_col <- max(1.5 * n_row, n_col)
      plate_range <-
        ifelse(n_col < 26, paste0("A1:", LETTERS[n_col], n_row), "A1:AV32")
      
      # need to adapt for 1536 well plates
      df_template <-
        base::expand.grid(WellRow = LETTERS[1:n_row], WellColumn = 1:n_col)
      
      for (iS in template_sheets[[iF]]) {
        tryCatch({
          df <- as.data.frame(readxl::read_excel(
            template_file[[iF]],
            sheet = iS,
            col_names = paste0("x", 1:n_col),
            range = plate_range
          ))
        }, error = function(e) {
          stop(sprintf("Error loading %s. Please check logs: %s", template_file[[iF]], e))
        })
        df$WellRow <- LETTERS[1:n_row]
        df_melted <- reshape2::melt(df, id.vars = "WellRow")
        
        # check if metadata field already exist and correct capitalization if needed
        if (!(iS %in% metadata_fields)) {
          if (!is.null(metadata_fields) &&
              toupper(iS) %in% toupper(metadata_fields)) {
            oldj <- iS
            iS <-
              metadata_fields[toupper(iS) == toupper(metadata_fields)]
            print(paste(oldj, "corrected to match case with ", iS))
          } else {
            metadata_fields <- c(metadata_fields, iS)
          }
        }
        colnames(df_melted)[3] <- iS
        colnames(df_melted)[colnames(df_melted) == "variable"] <-
          "WellColumn"
        df_melted$WellColumn <-
          gsub("x", "", df_melted$WellColumn)
        df_template <-
          base::merge(df_template, df_melted, by = c("WellRow", "WellColumn"))
      }
      df_template$Template <- template_filename[iF]
      colnames(df_template) <-
        check_metadata_names(colnames(df_template), log_str,
                             df_name = template_filename[iF])
      all_templates <- dplyr::bind_rows(all_templates, df_template)
      
    }
    print("Templates loaded successfully!")
    return(all_templates)
  }

#' Load results from tsv
#' 
#' This functions loads and checks the results file(s)
#' 
#' @param results_file character, file path(s) to template(s)
#' @param log_str character, file path to logs
load_results_tsv <-
  function(results_file, log_str) {
    # results_file is a string or a vector of strings
    log_str <- c(log_str, "", "load_results")
    
    results_filename <- basename(results_file)
    
    # read all files
    all_results <- data.frame()
    for (iF in 1:length(results_file)) {
      print(paste("Reading file", results_file[iF]))
      tryCatch({
        df <-
          readr::read_tsv(results_file[iF],
                          col_names = TRUE,
                          skip_empty_rows = TRUE)
      }, error = function(e) {
        stop(sprintf("Error reading %s", results_file[[iF]]))
      })
      # skip_empty_rows flag needs to be TRUE even if it ends up not skipping empty rows
      if (dim(df)[2] == 1) {
        tryCatch({
          # likely a csv file
          df <-
            read_csv(results_file[iF],
                     col_names = TRUE,
                     skip_empty_rows = TRUE)
        }, error = function(e) {
          stop(sprintf("Error reading %s", results_file[[iF]]))
        })
      }
      
      for (coln in c("Barcode",
                     get_identifier("WellPosition"),
                     "ReadoutValue")) {
        if (!(coln %in% colnames(df))) {
          ErrorMsg(coln, "needs to be a column of", results_filename[iF])
        }
      }
      if (dim(unique(df[, c("Barcode", get_identifier("WellPosition"))]))[1] !=
          dim(df[, c("Barcode", get_identifier("WellPosition"))])[1]) {
        ErrorMsg("Multiple rows with the same Barcode and Well in",
                 results_filename[iF])
      }
      if (!("BackgroundValue" %in% colnames(df)))
        df$BackgroundValue <- 0
      
      print(paste("File", results_filename[iF],
                  "read;", dim(df)[1], "wells"))
      all_results <- rbind(all_results, df)
      
      print("File done")
    }
    
    if (dim(unique(df[, c("Barcode", get_identifier("WellPosition"))]))[1] !=
        dim(df[, c("Barcode", get_identifier("WellPosition"))])[1]) {
      ErrorMsg("Multiple rows with the same Barcode and Well across all files")
    }
    
    return(all_results)
  }


#' Load results from xlsx
#' 
#' This functions loads and checks the results file(s)
#' 
#' @param results_file character, file path(s) to template(s)
#' @param log_str character, file path to logs
load_results_EnVision <-
  function(results_file, log_str) {
    results_filename <- basename(results_file)
    # results_file is a string or a vector of strings
    log_str <- c(log_str, "", "load_results")
    
    # test if the result files are .tsv or .xls(x) files
    isExcel <- sapply(results_file, function(x) {
      return(tools::file_ext(x) %in% c("xlsx", "xls"))
    })
    
    # read sheets in files; warning if more than one sheet (unexpected but can be handled)
    results_sheets <- vector("list", length(results_file))
    results_sheets[!isExcel] <- 0
    results_sheets[isExcel] <-
      lapply(results_file[isExcel], excel_sheets)
    if (any(lapply(results_sheets, length) > 1)) {
      WarnMsg <- paste("multiple sheets in result file:",
                       results_file[lapply(results_sheets, length) > 1])
      warning(WarnMsg)
    }
    
    # read all files and sheets
    all_results <- data.frame()
    for (iF in 1:length(results_file)) {
      for (iS in results_sheets[[iF]]) {
        print(paste("Reading file", results_file[[iF]], "; sheet", iS))
        if (iS == 0) {
          tryCatch({
            df <-
              readr::read_tsv(results_file[[iF]],
                              col_names = FALSE,
                              skip_empty_rows = TRUE)
          }, error = function(e) {
            stop(sprintf("Error reading %s, sheet %s", results_file[[iF]], iS))
          })
          # skip_empty_rows flag needs to be TRUE even if it ends up not skipping empty rows
          if (dim(df)[2] == 1) {
            tryCatch({
              # likely a csv file
              df <-
                readr::read_csv(results_file[[iF]],
                                col_names = FALSE,
                                skip_empty_rows = TRUE)
            }, error = function(e) {
              stop(sprintf("Error reading %s, sheet %s", results_file[[iF]], iS))
            })
          }
        } else {
          # expect an Excel spreadsheet
          if (length(results_sheets[[iF]]) > 1) {
            # if multiple sheets, assume 1 plate per sheet
            tryCatch({
              df <-
                readxl::read_excel(
                  results_file[[iF]],
                  sheet = iS,
                  col_names = paste0("x", 1:48),
                  range = "A1:AV32"
                )
            }, error = function(e) {
              stop(sprintf("Error reading %s, sheet %s", results_file[[iF]], iS))
            })
          } else {
            tryCatch({
              df <- readxl::read_excel(results_file[[iF]],
                                       sheet = iS,
                                       col_names = FALSE)
            }, error = function(e) {
              stop(sprintf("Error reading %s, sheet %s", results_file[[iF]], iS))
            })
            colnames(df) <-
              col_names <- paste0("x", 1:dim(df)[2])
          }
          df <-
            df[, !apply(df[1:48, ], 2, function(x) all(is.na(x)))] 
          # remove extra columns
          # limit to first 48 rows in case Protocol information is
          # exported which generate craps at the end of the file
        }
        full_rows <- !apply(df[,-6:-1], 1, function(x) all(is.na(x)))
        # not empty rows
        # before discarding the rows; move ''Background information'' in the next row
        Bckd_info_idx <- which(as.data.frame(df)[, 1] %in% 'Background information')
        if (length(Bckd_info_idx) > 0) {
          df[Bckd_info_idx + 1, 1] = df[Bckd_info_idx, 1]
          df[Bckd_info_idx, 1] = ''
        }
        
        # don't consider the first columns as these may be metadata
        # if big gap, delete what is at the bottom (Protocol information)
        gaps <-
          min(which(full_rows)[(diff(which(full_rows)) > 20)] + 1, dim(df)[1])
        df <-
          df[which(full_rows)[which(full_rows) <= gaps], ] # remove extra rows
        df <-
          df[, !apply(df, 2, function(x) all(is.na(x)))] 
        # remove empty columns
        
        # get the plate size
        n_col <- 1.5 * 2 ** ceiling(log2(dim(df)[2] / 1.5))
        n_row <- n_col / 1.5
        
        # get the barcode(s) in the sheet; expected in column C (third one)
        Barcode_idx <-
          which(as.data.frame(df)[, 3] %in% "Barcode")
        # run through all plates
        for (iB in Barcode_idx) {
          # two type of format depending on where Background information is placed
          if (any(as.data.frame(df)[iB + (1:4), 1] %in% "Background information")) {
            ref_bckgrd <-
              which(as.data.frame(df)[iB + (1:4), 1] %in% "Background information")
            readout_offset <- 1 + ref_bckgrd
            stopifnot(as.character(df[iB + ref_bckgrd, 4]) %in% 'Signal')
            BackgroundValue <-
              as.numeric(df[iB + ref_bckgrd + 1, 4])
          } else {
            # export without background information
            # case of " Exported with EnVision Workstation version 1.13.3009.1409 "
            readout_offset <- 1
            BackgroundValue <- 0
          }
          
          # check the structure of file is ok
          check_values <-
            as.matrix(df[iB + readout_offset + c(0, 1, n_row, n_row + 1), n_col])
          if (any(c(is.na(check_values[2:3]),!is.na(check_values[c(1, 4)])))) {
            ErrorMsg <-
              paste(
                "In result file",
                results_filename[[iF]],
                "(sheet",
                iS,
                ") readout values are misplaced for plate",
                as.character(df[iB + 1, 3])
              )
            stop(ErrorMsg)
          }
          
          readout <-
            as.matrix(df[iB + readout_offset + (1:n_row), 1:n_col])
          
          # check that the plate size is consistent and contains values
          if (any(is.na(readout))) {
            ErrorMsg <-
              paste(
                "In result file",
                results_filename[[iF]],
                "(sheet",
                iS,
                ") readout values are missing for plate",
                as.character(df[iB + 1, 3])
              )
            stop(ErrorMsg)
          }
          
          df_results <- data.frame(
            Barcode = as.character(df[iB + 1, 3]),
            WellRow = LETTERS[1:n_row],
            WellColumn = as.vector(t(matrix(
              1:n_col, n_col, n_row
            ))),
            ReadoutValue = as.numeric(as.vector(readout)),
            BackgroundValue = BackgroundValue
          )
          print(paste(
            "Plate",
            as.character(df[iB + 1, 3]),
            "read;",
            dim(df_results)[1],
            "wells"
          ))
          all_results <- rbind(all_results, df_results)
        }
        print("File done")
      }
    }
    return(all_results)
  }



#' @export
check_metadata_names <-
  function(col_df,
           log_str,
           df_name = "",
           df_type = NULL) {
    log_str <- c(log_str, "   check_metadata_names")
    # first check for required column names
    if (!is.null(df_type)) {
      if (df_type == "manifest") {
        expected_headers <- get_header("manifest")
      } else if (df_type == "template") {
        expected_headers <- get_identifier("drug")
      } else if (df_type == "template_treatment") {
        expected_headers <- c(get_identifier("drug"), "Concentration")
      }
      
      headersOK <- toupper(expected_headers) %in% toupper(col_df)
      if (any(!headersOK)) {
        ErrorMsg <- paste(
          df_name,
          "does not contains all expected headers for a",
          df_type,
          "; ",
          paste(expected_headers[!(expected_headers %in% col_df)], collpase = " ; "),
          " required"
        )
        log_str <- c(log_str, "Error in check_metadata_names:")
        log_str <- c(log_str, ErrorMsg)
        writeLines(log_str)
        stop(ErrorMsg)
      }
      if (df_type == "template_treatment") {
        # assess if multiple drugs and proper pairing
        n_drug <-
          agrep(get_identifier("drug"), col_df, ignore.case = TRUE)
        n_conc <-
          agrep("Concentration", col_df, ignore.case = TRUE)
        if (length(n_drug) != length(n_conc)) {
          ErrorMsg <- paste(
            "Treatment template",
            df_name,
            "does not contains the same number of Gnumber_* and Concentration_* sheets"
          )
          log_str <-
            c(log_str, "Error in check_metadata_names:")
          log_str <- c(log_str, ErrorMsg)
          writeLines(log_str)
          stop(ErrorMsg)
        }
        if (length(n_drug) > 1) {
          trt_sheets <- c(
            paste0(get_identifier("drug"), "_",
                   2:length(n_drug)),
            paste0("Concentration_", 2:length(n_conc))
          )
          if (!(all(toupper(trt_sheets) %in% toupper(col_df)))) {
            ErrorMsg <- paste(
              "Treatment template",
              df_name,
              "does not contains: ",
              paste(trt_sheets[!(toupper(trt_sheets) %in% toupper(col_df))],
                    collapse = " ; ")
            )
            log_str <-
              c(log_str, "Error in check_metadata_names:")
            log_str <- c(log_str, ErrorMsg)
            writeLines(log_str)
            stop(ErrorMsg)
          }
        }
      }
    }
    check_headers <-
      setdiff(get_header("reserved"), get_identifier("WellPosition"))
    
    
    corrected_names <- col_df
    
    # remove spaces and convert to WordUppercase
    names_spaces <- regexpr("\\s", corrected_names) > 0
    if (any(names_spaces)) {
      for (i in which(names_spaces)) {
        s <- strsplit(corrected_names[i], " ")[[1]]
        corrected_names[i] <-
          paste(toupper(substring(s, 1, 1)),
                substring(s, 2),
                sep = "",
                collapse = "")
      }
      
      WarnMsg <- paste(
        "Metadata field names for",
        df_name,
        "cannot contain spaces --> corrected to: ",
        paste(corrected_names[names_spaces], collapse = " ; ")
      )
      log_str <- c(log_str, "Warning in check_metadata_names:")
      log_str <- c(log_str, WarnMsg)
      
      warning(WarnMsg)
    }
    
    # check for wrong metadata field names (including dash, starting with number, ... )
    bad_names <-
      regexpr("\\W", corrected_names) > 0 |
      regexpr("\\d", corrected_names) == 1
    if (any(bad_names)) {
      ErrorMsg <- paste(
        "Metadata field names for",
        df_name,
        "cannot contain special characters or start with a number: ",
        paste(corrected_names[bad_names], collapse = " ; ")
      )
      log_str <- c(log_str, "Error in check_metadata_names:")
      log_str <- c(log_str, ErrorMsg)
      writeLines(log_str)
      stop(ErrorMsg)
    }
    
    # common headers that are written in a specific way
    # throw warning if close match and correct upper/lower case for consistency
    for (i in 1:length(get_header("controlled"))) {
      case_match <- setdiff(
        grep(paste0(get_header("controlled")[i], "$"), corrected_names, ignore.case = TRUE),
        grep(paste0(get_header("controlled")[i], "$"), corrected_names)
      )
      if (length(case_match) > 0) {
        WarnMsg <-
          paste(
            "Header",
            corrected_names[case_match],
            "in",
            df_name,
            "corrected to",
            get_header("controlled")[i]
          )
        corrected_names[case_match] <-
          get_header("controlled")[i]
        log_str <-
          c(log_str, "Warning in check_metadata_names:")
        log_str <- c(log_str, WarnMsg)
        warning(WarnMsg)
      }
    }
    
    # check for headers that are reserved for downstream analyses
    if (any(corrected_names %in% check_headers)) {
      ErrorMsg <- paste(
        "Metadata field name: ",
        paste(intersect(check_headers, corrected_names), collapse = " ; "),
        " in",
        df_name,
        "is not valid (reserved for output)"
      )
      log_str <- c(log_str, "Error in check_metadata_names:")
      log_str <- c(log_str, ErrorMsg)
      writeLines(log_str)
      stop(ErrorMsg)
    }
    
    return(corrected_names)
  }
