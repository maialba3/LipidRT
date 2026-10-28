# checksource
#' check source from data (lipidms, msdial, other)
#'
#' Check the source of the data (LipidMS, MS-DIAL, or other) and identify the
#' row where relevant data begins
#'
#' @param features Path to a file (txt or csv) containing the annotated peak list
#' extracted using external software (e.g., MS-DIAL)
#' @param sep Character used to separate columns in text (.txt) or CSV files. By default ",".
#' @param sheet Sheet to read. Either a string (the name of a sheet), or an integer
#' (the position of the sheet). Ignored if the sheet is specified via range.
#' If neither argument specifies the sheet, defaults to the first sheet.
#'
#' @return header_row and source
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checksource <- function(features, sep = ",", sheet = NULL) {
  # If the input is a file path
  if (is.character(features)) {
    extension <- tools::file_ext(features)

    # Read a preview of the file (first 10 rows).
    preview <- if (extension == "csv") {
      read.csv(features, header = FALSE, sep = sep, nrows = 20, check.names = FALSE)
    } else if (extension == "xlsx") {
      # Read Excel file without header
      readxl::read_excel(features, sheet = sheet, col_names = FALSE, n_max = 20)
    } else if (extension == "txt") {
      read.table(features, header = FALSE, sep = sep, nrows = 20, check.names = FALSE)
    }

    # Search for the header row and detect the source of the data
    header_info <- apply(preview, 1, function(row) {
      row <- as.character(row)
      if (all(c("Alignment ID", "Average Rt(min)", "Average Mz",
                "Metabolite name", "Adduct type") %in% row)) {
        return("msdial")
      } else if (all(c("mz", "RT", "LipidMSid", "Adduct") %in% row)) {
        return("lipidms")
      } else if (all(c("lipidname") %in% row) &&
                 any(c("rtmin", "rtsec") %in% tolower(row))) {
        return("other")
      } else {
        return(NA)
      }
    })

    # Retrieve the first valid row and identify the data source
    header_row <- which(!is.na(header_info))[1]
    source_data <- header_info[header_row]

    if (is.na(source_data) || is.null(header_row)) {
      stop("Data does not have the expected structure: 'mz', 'RTmin' or 'RTsec', 'lipidname', 'adduct'.")
    }

    return(list(source_data = source_data,
                header_row = header_row))

  } else if (is.data.frame(features)) {
    # If the input is already a data frame, detect the source based on the column names
    colnames_data <- colnames(features)

    if (all(c("Alignment ID", "Average Rt(min)", "Average Mz",
              "Metabolite name", "Adduct type") %in% colnames_data)) {
      source_data <- "msdial"
    } else if (all(c("mz", "RT", "LipidMSid", "Adduct") %in% colnames_data)) {
      source_data <- "lipidms"
    } else if (all(c("lipidname") %in% colnames_data) &&
               any(c("rtmin", "rtsec") %in% tolower(colnames_data))) {
      source_data <- "other"
    } else {
      stop("Data does not have the expected structure: 'mz', 'RTmin' or 'RTsec', 'lipidname', 'adduct'.")
    }
    return(source_data)

  } else {
    stop("Input 'data' must be a data.frame or a valid path to a .txt/.csv file.")
  }
}


# standardizedata
#' For each feature, separate all possible annotated lipids within each group.
#'
#' For each feature, separate all possible annotated lipids within each group, as well
#' as obtain the lipid class, the number of carbons and double bonds, and the
#' adduct of each one.
#'
#' @param features original peak list.
#' @param source source of data (lipidms, msdial o other). If other, the file
#' must have at least 5 columns: mz, RT, lipidname, adduct and confidenceLevel.
#'
#' @return processedData
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
standardizedata <- function(features, source) {

  # Keep only the relevant columns based on the detected data source
  data <- preparedata(features = features, source = source)

  # Unknown data
  data_isNA <- as.data.frame(data[is.na(data$lipidname) | data$lipidname %in% c("Unknown", ""), , drop = FALSE])

  processedData_isNA <- data_isNA[, c("mz", "RTsec", "LipidRTid", "lipidname")]
  colnames(processedData_isNA) <- c("mz", "RTsec", "LipidRTid", "raw_lipidname")


  processedData_isNA <- cbind(
    processedData_isNA,
    data.frame(
      lipidname       = character(nrow(processedData_isNA)),
      adduct          = character(nrow(processedData_isNA)),
      IS              = character(nrow(processedData_isNA)),
      confidenceLevel = character(nrow(processedData_isNA)),
      LipidClass      = character(nrow(processedData_isNA)),
      LipidSubclass   = character(nrow(processedData_isNA)),
      C               = numeric(nrow(processedData_isNA)),
      DB              = numeric(nrow(processedData_isNA)),
      Filter          = character(nrow(processedData_isNA)),
      roworder        = character(nrow(processedData_isNA)),
      stringsAsFactors = FALSE
    )
  )

  # annotated data
  processedData <- data.frame()
  separator <- "\\|"
  LipidRTids <- setdiff(unique(data$LipidRTid), data_isNA$LipidRTid)

  for (i in LipidRTids) {
    ss <- data[!is.na(data$LipidRTid) & data$LipidRTid == i, , drop = FALSE]
    raw_lipidname <- unlist(strsplit(ss$lipidname, split = separator))

    if (grepl("no MS2: |low score: ", ss$lipidname)) {

      # Handle cases with "no MS2: " or "low score: " annotations (MS-DIAL)
      if (grepl("no MS2: ", ss$lipidname)) {
        raw_lipidname <- sapply(raw_lipidname, function(x) {
          if (!startsWith(x, "no MS2: ")) paste0("no MS2: ", x) else x
        })
      } else if (grepl("low score: ", ss$lipidname)) {
        raw_lipidname <- sapply(raw_lipidname, function(x) {
          if (!startsWith(x, "low score: ")) paste0("low score: ", x) else x
        })
      }

      processed_data <- data.frame(
        mz = as.numeric(ss$mz),
        RTsec = as.numeric(ss$RTsec),
        LipidRTid = ss$LipidRTid,
        raw_lipidname = as.vector(raw_lipidname),
        lipidname = "",
        adduct = "",
        IS = ss$IS,
        confidenceLevel = "Not pass QC",
        LipidClass = "",
        LipidSubclass = "",
        C = NA,
        DB = NA,
        Filter = "Not pass QC",
        roworder = seq_along(raw_lipidname),
        stringsAsFactors = FALSE
      )


    } else {
      # Annotations (lipidname) ................................................
      lipidname <- unlist(strsplit(ss$lipidname, split = separator))

      # Identify lipid class and subclass ......................................
      class_info <- checkclass(ids = lipidname, general = FALSE)
      LipidClass <- as.character(class_info$lipidclass)
      LipidSubclass <- as.character(class_info$lipidsubclass)

      # Adduct .................................................................
      if (ss$adduct == "" | is.na(ss$adduct)) {
        adduct <- ""
      } else {
        adduct <- unlist(strsplit(ss$adduct, split = separator))
        adduct <- gsub("\\[|\\]2-|\\]-|\\]2\\+|\\]\\+", "", adduct)
      }

      # If only one adduct but multiple lipids, replicate adduct
      if (length(adduct) == 1 && length(lipidname) > 1) {
        adduct <- rep(adduct, length(lipidname))
      }

      # IS .....................................................................
      if (ss$IS == "" | is.na(ss$adduct)) {
        IS <- ""
      } else {
        IS <- unlist(strsplit(ss$IS, split = separator))
      }

      if (length(IS) == 1 && length(lipidname) > 1) {
        IS <- rep(IS, length(lipidname))
      }

      IS[is.na(IS) | IS == "NA"] <- ""

      # confidenceLevel ........................................................
      confidenceLevel <- unlist(strsplit(ss$confidenceLevel, split = separator))

      if (length(confidenceLevel) == 0) {
        confidenceLevel <- c()
        for (l in seq_along(lipidname)) {
          if (LipidClass[l] != "") {
            lipid_confidenceLevel <- determineconfidenceLevel(lipidname[l])
          } else {
            lipid_confidenceLevel <- "Unsupported class"
          }
          confidenceLevel <- c(confidenceLevel, lipid_confidenceLevel)
        }

      } else if ((length(confidenceLevel) == 1) && length(lipidname) > 1) {
        if (all(confidenceLevel %in% LipidRT::confLevels$level)) {
          confidenceLevel <- rep(confidenceLevel, length(lipidname))
        } else {
          confidenceLevel <- c()
          for (l in seq_along(lipidname)) {
            if (LipidClass[l] != "") {
              lipid_confidenceLevel <- determineconfidenceLevel(lipidname[l])
            } else {
              lipid_confidenceLevel <- "Unsupported class"
            }
            confidenceLevel <- c(confidenceLevel, lipid_confidenceLevel)
          }
        }
      }

      if (length(confidenceLevel) == 1 && length(lipidname) > 1) {
        confidenceLevel <- rep(confidenceLevel, length(lipidname))
      }

      # C, DB, Filter ..........................................................
      C <- c()
      DB <- c()
      Filter <- c()

      lipid <- lipidname

      for (l in seq_along(lipidname)) {
        if (LipidClass[l] != "") {

          lipid[l] <- lipidname[l]
          lipid[l] <- convertlipidannotation(lipid = lipid[l],
                                             lipidclass = LipidClass[l],
                                             lipidsubclass = LipidSubclass[l])
          Filter[l] <- "Kept"


          # Extract carbons and double bonds
          c_db <- gsub("no MS2: |low score: ", "", lipid[l])
          c_db <- gsub("^[A-Za-z0-9_\\-]+\\(", "", c_db)
          c_db <- gsub("\\)$", "", c_db)

          # Remove d/t only when they are part of shorthand sphingo notation
          c_db <- gsub(
            "(^|[/_\\-])(d|t|O-|P-)(?=[0-9]+:[0-9]+)",
            "\\1",
            c_db,
            perl = TRUE
          )

          c_db <- unlist(strsplit(c_db, "\\(|\\)|/|:|_|-"))
          c_db <- c_db[c_db != ""]

          if (length(c_db) >= 2) {
            c_part <- suppressWarnings(as.numeric(c_db[c(1:length(c_db)) %% 2 != 0]))
            db_part <- suppressWarnings(as.numeric(c_db[c(1:length(c_db)) %% 2 == 0]))

            C[l] <- if (all(is.na(c_part))) NA else sum(c_part, na.rm = TRUE)
            DB[l] <- if (all(is.na(db_part))) NA else sum(db_part, na.rm = TRUE)
          } else {
            C[l] <- NA
            DB[l] <- NA
          }

        } else {
          LipidClass[l] <- ""
          LipidSubclass[l] <- ""
          C[l] <- NA
          DB[l] <- NA
          Filter[l] <- "Unsupported class"
          confidenceLevel[l] <- "Unsupported class"
        }
      }

      # Si hay mas de una subclase y son PX o LPX

      if (length(unique(LipidSubclass)) > 1) {

      }

      # ........................................................................
      processed_data <- data.frame(
        mz = as.numeric(ss$mz),
        RTsec = as.numeric(ss$RTsec),
        LipidRTid = ss$LipidRTid,
        raw_lipidname = as.vector(lipidname),
        lipidname = as.vector(lipid),
        adduct = as.vector(adduct),
        IS = as.vector(IS),
        confidenceLevel = as.vector(confidenceLevel),
        LipidClass = as.vector(LipidClass),
        LipidSubclass = as.vector(LipidSubclass),
        C = as.numeric(C),
        DB = as.numeric(DB),
        Filter = as.vector(Filter),
        roworder = seq_along(lipidname),
        stringsAsFactors = FALSE
      )

      # Mark duplicated annotations within each group
      confLevels <- LipidRT::confLevels
      confLevels$isobaric_order <- confLevels$order
      confLevels$isobaric_order[confLevels$isobaric_order == 4] <- 3

      if (nrow(processed_data) > 0) {
        split_data <- split(
          processed_data,
          paste(processed_data$LipidSubclass, processed_data$adduct, sep = "_")
        )

        split_data <- lapply(split_data, function(df) {
          if (nrow(df) > 0 && all(!(df$Filter %in% c("Unsupported class", "Not pass QC")))) {
            df$conflev <- confLevels$order[match(df$confidenceLevel, confLevels$level)]
            df <- df[order(-df$conflev), , drop = FALSE]
            if (nrow(df) > 1) {
              df$Filter[2:nrow(df)] <- "Ignore duplicated"
            }
            df$conflev <- NULL
            df <- df[order(df$roworder), , drop = FALSE]
          }
          return(df)
        })

        processed_data <- do.call(rbind, split_data)
        processed_data <- processed_data[order(processed_data$roworder), , drop = FALSE]

      }

    }

    processedData <- rbind(processedData, processed_data)

  }

  processedTotalData <- rbind(processedData, processedData_isNA)
  processedTotalData <- processedTotalData[order(processedTotalData$LipidRTid), ]

  # Add index column
  processedTotalData$index <- seq_len(nrow(processedTotalData))

  # Keep only valid lipid classes
  processedTotalData <- processedTotalData[!is.na(processedTotalData$LipidSubclass), , drop = FALSE]


  processedTotalData$confidenceLevel <- factor(
    processedTotalData$confidenceLevel,
    levels = c("FA position", "FA", "Subclass", "MS-only", "MSMS", "Predicted",
               "Other", "Unsupported class", "Not pass QC", ""),
    ordered = TRUE)

  return(processedTotalData)
}



# preparedata
#' prepare data table based on source format
#'
#' Cleans the feature table depending on the input source format
#' (e.g., LipidMS, MSDIAL, or other). Filters out invalid or unknown entries,
#' selects relevant columns, and unifies column names for downstream processing.
#'
#' @param features original peak list.
#' @param source source of data (lipidms, msdial o other). If other, the file
#' must have at least 5 columns: mz, RT, lipidname, adduct and confidenceLevel.
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
preparedata <- function(features, source) {

  if (!("LipidRTid" %in% colnames(features))) {
    features$LipidRTid <- seq_len(nrow(features))
  }

  if (!"IS" %in% colnames(features)) {
    features$IS <- rep("", nrow(features))
  } else {
    features$IS[is.na(features$IS)] <- ""
  }

  if (source == "lipidms") { # LipidMS .........................................

    if (!("confidenceLevel" %in% colnames(features))) {
      features$confidenceLevel <- ""
    }

    data <- features[(features$LipidMSid != "" &
                             !is.na(features$LipidMSid)) |
                            features$isotope == "[M+0]",
                          c("mz", "RT", "LipidRTid", "LipidMSid",
                            "Adduct", "confidenceLevel", "IS"),
                          drop = FALSE]


  } else if (source == "msdial") { # MSDIAL ....................................

    if (!("confidenceLevel" %in% colnames(features))) {
      if ("Annotation tag (VS1.0)" %in% colnames(features)) {

        annotation <- features$`Annotation tag (VS1.0)`

        features$confidenceLevel <- ""

        features$confidenceLevel[annotation %in% c("Lipid assigned", "Chain resolved")] <- "FA position"
        features$confidenceLevel[annotation == "Chain semi-resolved"] <- "FA"
        features$confidenceLevel[annotation == "Class resolved"] <- "Subclass"
        features$confidenceLevel[annotation %in% c(
          "Annotated by experimental MS/MS",
          "annotated by user-defined text library",
          "Unknown"
        )] <- "MSMS"

      } else {
        features$confidenceLevel <- ""
      }
    }


    # Convert RT to seconds and ensure it is numeric
    features$RTsec <- as.numeric(gsub(",", ".", features$`Average Rt(min)`)) * 60
    # Convert mz values to numeric format.
    features$mz <- as.numeric(gsub(",", ".", features$`Average Mz`))

    data <- features[((features$`Metabolite name` != "Unknown" &
           !is.na(features$`Metabolite name`))) |
             (!grepl("^\\s*([1-9]\\d*|0*\\.\\d*[1-9]\\d*)\\s*$",
                     features$`Isotope tracking weight number`)),
           c("mz", "RTsec", "LipidRTid", "Metabolite name",
             "Adduct type", "confidenceLevel", "IS"),
           drop = FALSE]


  } else if (source == "other") { # Other ......................................

    if ("RTmin" %in% colnames(features) &&
        !("RTsec" %in% colnames(features))) {
      features$RTsec <- features$RTmin*60
    }

    if (!"mz" %in% colnames(features)) {
      features$mz <- rep(NA, nrow(features))
    }

    if (!"adduct" %in% colnames(features)) {
      features$adduct <- rep("", nrow(features))
    }

    if (!("confidenceLevel" %in% colnames(features))) {
      features$confidenceLevel <- ""
    }


    data <- features[, c("mz", "RTsec", "LipidRTid", "lipidname",
                            "adduct", "confidenceLevel", "IS"), drop = FALSE]
  } # ..........................................................................

  # Standardize column names.
  colnames(data) <- c("mz", "RTsec", "LipidRTid", "lipidname", "adduct",
                      "confidenceLevel", "IS")

  # convert to numeric: mz, RT
  data$mz <- as.numeric(data$mz)
  data$RTsec <- as.numeric(data$RTsec)

  return(data)
}





# determineconfidenceLevel
#' Identify the confidence level of a lipid annotation.
#'
#' The confidence level of a lipid annotation is determined based on how the fatty
#' acid chains are separated within the structure.
#' If the chains are separated by a slash (/), it indicates FA position level.
#' If they are separated by an underscore (_), it corresponds to FA level.
#' If there is no separator, the annotation is considered to be at the Subclass level.
#'
#' @param lipid character containing the converted annotation.
#'
#' @return confidenceLevel
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
determineconfidenceLevel <- function(lipid) {
    # extract content within brackets
    bracket_content <- gsub(".*\\((.*)\\).*", "\\1", lipid)

    # assign annotation level based on lipid name
    if (grepl("/", bracket_content)) {
      confidenceLevel <- "FA position"
    } else if (grepl("_", bracket_content)) {
      confidenceLevel <- "FA"
    } else {
      confidenceLevel <- "Subclass"
    }
    return(confidenceLevel)
}

# norm
#' Normalize lipid annotations
#'
#' This helper standardizes lipid annotations before class recognition and
#' conversion. It harmonizes aliases, spacing, oxygen notation and common
#' shorthand names.
#'
#' @param x Character string representing the lipid annotation to be normalized.
#'
#' @return A character string with the normalized lipid annotation.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
norm <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)

  # Normalize oxygen notation: ;2O -> ;O2
  x <- gsub(";([0-9]+)O\\b", ";O\\1", x, perl = TRUE)

  # Normalize separators
  x <- gsub("\u2212", "-", x, fixed = TRUE)  # minus sign
  x <- gsub("\u2013", "-", x, fixed = TRUE)  # en dash

  # Remove score prefixes
  x <- gsub("^no MS2: ", "", x, perl = TRUE)
  x <- gsub("^low score: ", "", x, perl = TRUE)

  # Common aliases
  x <- gsub("^CER\\b", "Cer", x, perl = TRUE)
  x <- gsub("^SPB\\b", "Sph", x, perl = TRUE)
  x <- gsub("^SPH\\b", "Sph", x, perl = TRUE)

  # Ceramide-related aliases
  x <- gsub("^Cer1P\\b", "CerP", x, perl = TRUE)
  x <- gsub("^GlcCer\\b|^GalCer\\b|^Hex-Cer\\b|^HexCer\\b", "HexCer", x, perl = TRUE)

  # Only convert ST when it is clearly used as a class alias
  x <- gsub("^Sulfatide(?=\\s|\\()|^ST(?=\\s|\\()", "SHexCer", x, perl = TRUE)

  x <- gsub("^AHexCer\\b", "AHexCer", x, perl = TRUE)
  x <- gsub("^ACER\\b|^ACer\\b", "AcylCer", x, perl = TRUE)

  # Carnitines
  x <- gsub("^CAR\\b|^ACAR\\b|^ACar\\b|^Car\\b", "Carnitine", x, perl = TRUE)
  x <- gsub("^(\\d+:\\d+)\\s+Carnitine$", "Carnitine \\1", x, perl = TRUE)

  # Ether shorthand
  x <- gsub("PC\\bO-", "PC-O\\b", x, perl = TRUE)
  x <- gsub("PC\\bP-", "PC-P\\b", x, perl = TRUE)
  x <- gsub("PE\\bO-", "PE-O\\b", x, perl = TRUE)
  x <- gsub("PE\\bP-", "PE-P\\b", x, perl = TRUE)
  x <- gsub("PG\\bO-", "PG-O\\b", x, perl = TRUE)
  x <- gsub("PG\\bP-", "PG-P\\b", x, perl = TRUE)
  x <- gsub("PI\\bO-", "PI-O\\b", x, perl = TRUE)
  x <- gsub("PI\\bP-", "PI-P\\b", x, perl = TRUE)
  x <- gsub("PS\\bO-", "PS-O\\b", x, perl = TRUE)
  x <- gsub("PS\\bP-", "PS-P\\b", x, perl = TRUE)
  x <- gsub("PA\\bO-", "PA-O\\b", x, perl = TRUE)
  x <- gsub("PA\\bP-", "PA-P\\b", x, perl = TRUE)
  x <- gsub("LPC\\bO-", "LPC-O\\b", x, perl = TRUE)
  x <- gsub("LPC\\bP-", "LPC-P\\b", x, perl = TRUE)
  x <- gsub("LPE\\bO-", "LPE-O\\b", x, perl = TRUE)
  x <- gsub("LPE\\bP-", "LPE-P\\b", x, perl = TRUE)

  x <- gsub("PC\\(O-", "PC-O(", x, perl = TRUE)
  x <- gsub("PC\\(P-", "PC-P(", x, perl = TRUE)
  x <- gsub("PE\\(O-", "PE-O(", x, perl = TRUE)
  x <- gsub("PE\\(P-", "PE-P(", x, perl = TRUE)
  x <- gsub("PG\\(O-", "PG-O(", x, perl = TRUE)
  x <- gsub("PG\\(P-", "PG-P(", x, perl = TRUE)
  x <- gsub("PI\\(O-", "PI-O(", x, perl = TRUE)
  x <- gsub("PI\\(P-", "PI-P(", x, perl = TRUE)
  x <- gsub("PS\\(O-", "PS-O(", x, perl = TRUE)
  x <- gsub("PS\\(P-", "PS-P(", x, perl = TRUE)
  x <- gsub("PA\\(O-", "PA-O(", x, perl = TRUE)
  x <- gsub("PA\\(P-", "PA-P(", x, perl = TRUE)
  x <- gsub("LPC\\(O-", "LPC-O(", x, perl = TRUE)
  x <- gsub("LPC\\(P-", "LPC-P(", x, perl = TRUE)
  x <- gsub("LPE\\(O-", "LPE-O(", x, perl = TRUE)
  x <- gsub("LPE\\(P-", "LPE-P(", x, perl = TRUE)

  x <- gsub("PCo\\(", "PC-O(", x, perl = TRUE)
  x <- gsub("PCp\\(", "PC-P(", x, perl = TRUE)
  x <- gsub("PEo\\(", "PE-O(", x, perl = TRUE)
  x <- gsub("PEp\\(", "PE-P(", x, perl = TRUE)
  x <- gsub("PGo\\(", "PG-O(", x, perl = TRUE)
  x <- gsub("PGp\\(", "PG-P(", x, perl = TRUE)
  x <- gsub("PIo\\(", "PI-O(", x, perl = TRUE)
  x <- gsub("PIp\\(", "PI-P(", x, perl = TRUE)
  x <- gsub("PSo\\(", "PS-O(", x, perl = TRUE)
  x <- gsub("PSp\\(", "PS-P(", x, perl = TRUE)
  x <- gsub("PAo\\(", "PA-O(", x, perl = TRUE)
  x <- gsub("PAp\\(", "PA-P(", x, perl = TRUE)
  x <- gsub("LPCo\\(", "LPC-O(", x, perl = TRUE)
  x <- gsub("LPCp\\(", "LPC-P(", x, perl = TRUE)
  x <- gsub("LPEo\\(", "LPE-O(", x, perl = TRUE)
  x <- gsub("LPEp\\(", "LPE-P(", x, perl = TRUE)

  # x <- gsub("^LPC\\(O-", "LPC-O(", x, perl = TRUE)
  # x <- gsub("^LPC\\(P-", "LPC-P(", x, perl = TRUE)
  # x <- gsub("^LPE\\(O-", "LPE-O(", x, perl = TRUE)
  # x <- gsub("^LPE\\(P-", "LPE-P(", x, perl = TRUE)
  # x <- gsub("^LPG\\(O-", "LPG-O(", x, perl = TRUE)
  # x <- gsub("^LPG\\(P-", "LPG-P(", x, perl = TRUE)
  # x <- gsub("^LPI\\(O-", "LPI-O(", x, perl = TRUE)
  # x <- gsub("^LPI\\(P-", "LPI-P(", x, perl = TRUE)
  # x <- gsub("^LPS\\(O-", "LPS-O(", x, perl = TRUE)
  # x <- gsub("^LPS\\(P-", "LPS-P(", x, perl = TRUE)
  # x <- gsub("^LPA\\(O-", "LPA-O(", x, perl = TRUE)
  # x <- gsub("^LPA\\(P-", "LPA-P(", x, perl = TRUE)
  #
  # x <- gsub("^LPCo\\(", "LPC-O(", x, perl = TRUE)
  # x <- gsub("^LPCp\\(", "LPC-P(", x, perl = TRUE)
  # x <- gsub("^LPEo\\(", "LPE-O(", x, perl = TRUE)
  # x <- gsub("^LPEp\\(", "LPE-P(", x, perl = TRUE)
  # x <- gsub("^LPGo\\(", "LPG-O(", x, perl = TRUE)
  # x <- gsub("^LPGp\\(", "LPG-P(", x, perl = TRUE)
  # x <- gsub("^LPIo\\(", "LPI-O(", x, perl = TRUE)
  # x <- gsub("^LPIp\\(", "LPI-P(", x, perl = TRUE)
  # x <- gsub("^LPSo\\(", "LPS-O(", x, perl = TRUE)
  # x <- gsub("^LPSp\\(", "LPS-P(", x, perl = TRUE)
  # x <- gsub("^LPAo\\(", "LPA-O(", x, perl = TRUE)
  # x <- gsub("^LPAp\\(", "LPA-P(", x, perl = TRUE)

  x <- gsub("TG\\(O-", "TG-O\\b", x, perl = TRUE)
  x <- gsub("TG\\(P-", "TG-P\\b", x, perl = TRUE)
  x <- gsub("DG\\(O-", "DG-O\\b", x, perl = TRUE)
  x <- gsub("DG\\(P-", "DG-P\\b", x, perl = TRUE)

  x <- gsub("TG\\bO-", "TG-O\\b", x, perl = TRUE)
  x <- gsub("TG\\bP-", "TG-P\\b", x, perl = TRUE)
  x <- gsub("DG\\bO-", "DG-O\\b", x, perl = TRUE)
  x <- gsub("DG\\bP-", "DG-P\\b", x, perl = TRUE)

  return(x)
}



# is_deuterated_annotation
#' Detect isotopically labelled annotations
#'
#' This helper detects deuterated annotations that should be excluded from
#' class recognition and conversion.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return Logical value.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
is_deuterated_annotation <- function(x) {
  x <- norm(x)
  r <- grepl("\\(d\\d+\\)", x, perl = TRUE) ||
    grepl("\\[[^]]*d\\d+[^]]*\\]", x, perl = TRUE)
  return(r)
}


# is_oxidized_non_sphingo_pl
#' Detect oxidized non-sphingolipid phospholipids
#'
#' This helper identifies oxidized phospholipids and lysophospholipids that
#' should be excluded.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return Logical value.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
is_oxidized_non_sphingo_pl <- function(x) {
  x <- norm(x)

  is_pl <- grepl("^(PC|PE|PG|PI|PS|PA|LPC|LPE|LPG|LPI|LPS|PC-O|PC-P|PE-O|PE-P|LPC-O|LPC-P|LPE-O|LPE-P)(\\(|\\s)", x, perl = TRUE) ||
    grepl("^(PC|PE|LPC|LPE)\\s+[OP]-", x, perl = TRUE)

  is_gly <- grepl("^(MGDG|DGDG|MG|DG|TG)(\\(|\\s)", x, perl = TRUE)

  has_ox <- grepl(";O\\d*\\b", x, perl = TRUE) ||
    grepl("\\([^)]*(OH|OOH|oxo)[^)]*\\)", x, perl = TRUE)

  r <- (is_pl || is_gly) && has_ox

  return(r)
}


# has_modified_fa_sphingo
#' Detect unsupported modified fatty acyl chains in sphingolipids
#'
#' This helper identifies oxidized or hydroxylated fatty acyl chains in
#' sphingolipid annotations that should be excluded.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return Logical value.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
has_modified_fa_sphingo <- function(x) {
  x <- norm(x)

  r <- grepl("/[^/]*;O\\d*\\b", x, perl = TRUE) ||
    grepl("/[^/]*\\([^)]*(OH|OOH|oxo|keto)[^)]*\\)", x, perl = TRUE) ||
    grepl("\\(FA[^)]*\\)$", x, perl = TRUE) # "...(FA C:DB)"

  return(r)
}


# invalid_sphingo_annotation
#' Detect unsupported sphingolipid annotations
#'
#' This helper excludes unsupported sphingolipid annotations, such as classes
#' with invalid total oxygen counts.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return Logical value.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
invalid_sphingo_annotation <- function(x) {
  x <- norm(x)

  invalid_total_o <- grepl("^(Cer|HexCer|SM|CerP|Sph|AHexCer|SHexCer|SphP|AcylCer)\\s+\\d+:\\d+;O[4-9]\\d*(?:$|[ )])", x, perl = TRUE) ||
    grepl("^(Cer|HexCer|SM|CerP|Sph)\\(\\d+:\\d+;O[4-9]\\d*\\)$", x, perl = TRUE)

  return(invalid_total_o)
}


# base_lipidclass
#' Get base lipid class from lipid subclass
#'
#' This helper maps sphingolipid subclasses carrying d/t information to their
#' base lipid class.
#'
#' @param lipidsubclass Character string representing the lipid subclass.
#'
#' @return Character string with the base lipid class.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
base_lipidclass <- function(lipidsubclass) {

  if (lipidsubclass %in% c("Cer_d", "Cer_t")) {
    return("Cer")
  } else if (lipidsubclass %in% c("HexCer_d", "HexCer_t")) {
    return("HexCer")
  } else if (lipidsubclass %in% c("AHexCer_d", "AHexCer_t")) {
    return("AHexCer")
  } else if (lipidsubclass %in% c("SHexCer_d", "SHexCer_t")) {
    return("SHexCer")
  } else if (lipidsubclass %in% c("SM_d", "SM_t")) {
    return("SM")
  } else if (lipidsubclass %in% c("Sph_d", "Sph_t")) {
    return("Sph")
  } else if (lipidsubclass %in% c("SphP_d", "SphP_t")) {
    return("SphP")
  } else if (lipidsubclass %in% c("CerP_d", "CerP_t")) {
    return("CerP")
  } else if (lipidsubclass %in% c("AcylCer_d", "AcylCer_t")) {
    return("AcylCer")
  } else {
    return(lipidsubclass)
  }
}


# sphingo_base_type
#' Infer sphingoid base type
#'
#' This helper infers whether a sphingolipid belongs to a d- or t-type
#' sphingoid base.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return A character string with the inferred base type: \code{"d"},
#' \code{"t"}, or \code{""}.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
sphingo_base_type <- function(x) {

  x <- norm(x)

  # Explicit d/t always wins
  if (grepl("(^|[\\( /_-])t\\d+:\\d+\\b", x, perl = TRUE)) {
    return("t")
  }
  if (grepl("(^|[\\( /_-])d\\d+:\\d+\\b", x, perl = TRUE)) {
    return("d")
  }

  # MS-DIAL class tags
  if (grepl("-(NP|AP|BDS|BP)\\b", x, perl = TRUE)) {
    return("t")
  }
  if (grepl("-(NS|NDS|AS|ADS|HDS|EOS)\\b", x, perl = TRUE)) {
    return("d")
  }

  # Do not infer d/t if FA is oxidized or hydroxylated
  if (has_modified_fa_sphingo(x)) {
    return("")
  }

  # Molecular species: base oxygenation
  if (grepl("\\d+:\\d+;O3/\\d+:\\d+(?:$|[) ])", x, perl = TRUE)) {
    return("t")
  }
  if (grepl("\\d+:\\d+;O2/\\d+:\\d+(?:$|[) ])", x, perl = TRUE)) {
    return("d")
  }

  # Total composition
  if (grepl("\\d+:\\d+;O3(?:$|[) ])", x, perl = TRUE)) {
    return("t")
  }
  if (grepl("\\d+:\\d+;O2(?:$|[) ])", x, perl = TRUE)) {
    return("d")
  }

  return("")
}


# guess_sphingo_subclass
#' Guess sphingolipid subclass from an annotation
#'
#' This helper identifies sphingolipid subclasses, including d/t subclasses for
#' ceramides and glycosphingolipids.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return A character string with the guessed sphingolipid subclass.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
guess_sphingo_subclass <- function(x) {

  x <- norm(x)

  if (invalid_sphingo_annotation(x)) {
    return("")
  }

  base_type <- sphingo_base_type(x)

  if (grepl("^SphP(\\(|\\s)", x, perl = TRUE)) {
    if (base_type == "t") return("SphP_t")
    return("SphP_d")
  }

  if (grepl("^Sph(\\(|\\s)", x, perl = TRUE)) {
    if (base_type == "t") return("Sph_t")
    return("Sph_d")
  }

  if (grepl("^CerP(\\(|\\s)", x, perl = TRUE)) {
    if (base_type == "t") return("CerP_t")
    return("CerP_d")
  }

  if (grepl("^Cer\\s+P-", x, perl = TRUE)) {
    if (base_type == "t") return("CerP_t")
    return("CerP_d")
  }

  if (grepl("^SM(\\(|\\s)", x, perl = TRUE)) {

    if (has_modified_fa_sphingo(x)) {
      return("")
    }

    if (base_type == "t") return("SM_t")
    return("SM_d")
  }

  if (grepl("^AcylCer(\\(|\\s)", x, perl = TRUE)) {
    if (base_type == "t") return("AcylCer_t")
    return("AcylCer_d")
  }

  if (grepl("^SHexCer(\\(|\\s|-)", x, perl = TRUE)) {

    if (has_modified_fa_sphingo(x)) {
      return("")
    }

    valid_shexcer <- (
      grepl("^SHexCer\\([dt]\\d+:\\d+\\)$", x, perl = TRUE) ||
        grepl("^SHexCer\\([dt]\\d+:\\d+[/_\\-]\\d+:\\d+\\)$",
              x, perl = TRUE) ||
        grepl("^SHexCer\\([0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^SHexCer\\([0-9]+:[0-9]+[/_\\-][0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^SHexCer\\s+[dt]?\\d+:\\d+$",
              x, perl = TRUE) ||
        grepl("^SHexCer\\s+[0-9]+:[0-9]+;O[23]$",
              x, perl = TRUE) ||
        grepl("^SHexCer\\s+[0-9]+:[0-9]+;O[23]/[0-9]+:[0-9]+$",
              x, perl = TRUE) ||
        grepl("^SHexCer-[A-Z]+\\s+",
              x, perl = TRUE)
    )

    if (!valid_shexcer) {
      return("")
    }

    if (base_type == "t") return("SHexCer_t")
    return("SHexCer_d")
  }


  if (grepl("^AHexCer(\\(|\\s|-)", x, perl = TRUE)) {

    if (has_modified_fa_sphingo(x)) {
      return("")
    }

    valid_ahexcer <- (
      grepl("^AHexCer\\([dt]\\d+:\\d+\\)$",
            x, perl = TRUE) ||
        grepl("^AHexCer\\([dt]\\d+:\\d+[/_\\-]\\d+:\\d+\\)$",
              x, perl = TRUE) ||
        grepl("^AHexCer\\([0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^AHexCer\\([0-9]+:[0-9]+[/_\\-][0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^AHexCer\\s+[dt]?\\d+:\\d+$",
              x, perl = TRUE) ||
        grepl("^AHexCer\\s+[0-9]+:[0-9]+;O[23]$",
              x, perl = TRUE) ||
        grepl("^AHexCer\\s+[0-9]+:[0-9]+;O[23]/[0-9]+:[0-9]+$",
              x, perl = TRUE) ||
        grepl("^AHexCer-[A-Z]+\\s+",
              x, perl = TRUE)
    )

    if (!valid_ahexcer) {
      return("")
    }

    if (base_type == "t") return("AHexCer_t")
    return("AHexCer_d")
  }


  if (grepl("^HexCer(\\(|\\s|-)", x, perl = TRUE)) {

    if (has_modified_fa_sphingo(x)) {
      return("")
    }

    valid_hexcer <- (
      grepl("^HexCer\\([dt]\\d+:\\d+\\)$",
            x, perl = TRUE) ||
        grepl("^HexCer\\([dt]\\d+:\\d+[/_\\-]\\d+:\\d+\\)$",
              x, perl = TRUE) ||
        grepl("^HexCer\\([0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^HexCer\\([0-9]+:[0-9]+[/_\\-][0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^HexCer\\s+[dt]?\\d+:\\d+$",
              x, perl = TRUE) ||
        grepl("^HexCer\\s+[0-9]+:[0-9]+;O[23]$",
              x, perl = TRUE) ||
        grepl("^HexCer\\s+[0-9]+:[0-9]+;O[23]/[0-9]+:[0-9]+$",
              x, perl = TRUE) ||
        grepl("^HexCer-[A-Z]+\\s+",
              x, perl = TRUE)
    )

    if (!valid_hexcer) {
      return("")
    }

    if (base_type == "t") return("HexCer_t")
    return("HexCer_d")
  }


  if (grepl("^Cer(\\(|\\s|-)", x, perl = TRUE)) {

    if (has_modified_fa_sphingo(x)) {
      return("")
    }

    valid_cer <- (
      grepl("^Cer\\([dt]\\d+:\\d+\\)$",
            x, perl = TRUE) ||
        grepl("^Cer\\([dt]\\d+:\\d+[/_\\-]\\d+:\\d+\\)$",
              x, perl = TRUE) ||
        grepl("^Cer\\([0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^Cer\\([0-9]+:[0-9]+[/_\\-][0-9]+:[0-9]+\\)$",
              x, perl = TRUE) ||
        grepl("^Cer\\s+[dt]?\\d+:\\d+$",
              x, perl = TRUE) ||
        grepl("^Cer\\s+[0-9]+:[0-9]+;O[23]$",
              x, perl = TRUE) ||
        grepl("^Cer\\s+[0-9]+:[0-9]+;O[23]/[0-9]+:[0-9]+$",
              x, perl = TRUE) ||
        grepl("^Cer-[A-Z]+\\s+",
              x, perl = TRUE)
    )

    if (!valid_cer) {
      return("")
    }

    if (base_type == "t") return("Cer_t")
    return("Cer_d")
  }


  return("")
}

# guess_sphingo_class
#' Guess sphingolipid base class from an annotation
#'
#' This helper returns the base sphingolipid class, collapsing d/t subclasses
#' into their parent class.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return A character string with the guessed sphingolipid class.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
guess_sphingo_class <- function(x) {
  sph_subclass <- guess_sphingo_subclass(x)

  if (sph_subclass == "") {
    return("")
  }

  return(base_lipidclass(sph_subclass))
}


# guess_lipid_class
#' Guess lipid class from an annotation
#'
#' This helper identifies the lipid base class from an annotation.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return A character string with the guessed lipid class.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
guess_lipid_class <- function(x) {

  x <- norm(x)

  if (x == "" || is.na(x)) {
    return("")
  }

  # Reject isotopically labelled compounds
  if (is_deuterated_annotation(x)) {
    return("")
  }

  # Reject oxidized PL
  if (is_oxidized_non_sphingo_pl(x)) {
    return("")
  }

  # Reject unsupported ether glycoglycerolipids and neutral glycerolipids
  if (grepl("^(MGDG|DGDG|MG|DG|TG)(\\(|\\s)", x, perl = TRUE) &&
      grepl("(O-|P-|O_|P_)", x, perl = TRUE)) {
    return("")
  }

  # Sphingolipids
  sph_class <- guess_sphingo_class(x)
  if (sph_class != "") {
    return(sph_class)
  }

  # Fatty acyls and related
  if (grepl("^FAHFA(\\(|\\s)", x, perl = TRUE)) return("FAHFA")
  if (grepl("^HFA(\\(|\\s)", x, perl = TRUE)) return("HFA")
  if (grepl("^FA(\\(|\\s)", x, perl = TRUE)) return("FA")
  if (grepl("^NAGlySer(\\(|\\s)", x, perl = TRUE)) return("NAGlySer")
  if (grepl("^NAGly(\\(|\\s)", x, perl = TRUE)) return("NAGly")
  if (grepl("^NAOrn(\\(|\\s)", x, perl = TRUE)) return("NAOrn")
  if (grepl("^NAE(\\(|\\s)", x, perl = TRUE)) return("NAE")

  # Lysophospholipid ethers first
  if (grepl("^LPC\\s*\\(?O-|^LPC-O|^LPCo\\(", x, perl = TRUE)) return("LPC-O")
  if (grepl("^LPC\\s*\\(?P-|^LPC-P|^LPCp\\(", x, perl = TRUE)) return("LPC-P")
  if (grepl("^LPE\\s*\\(?O-|^LPE-O|^LPEo\\(", x, perl = TRUE)) return("LPE-O")
  if (grepl("^LPE\\s*\\(?P-|^LPE-P|^LPEp\\(", x, perl = TRUE)) return("LPE-P")
  if (grepl("^LPG\\s*\\(?O-|^LPG-O|^LPGo\\(", x, perl = TRUE)) return("LPG-O")
  if (grepl("^LPG\\s*\\(?P-|^LPG-P|^LPGp\\(", x, perl = TRUE)) return("LPG-P")
  if (grepl("^LPI\\s*\\(?O-|^LPI-O|^LPIo\\(", x, perl = TRUE)) return("LPI-O")
  if (grepl("^LPI\\s*\\(?P-|^LPI-P|^LPIp\\(", x, perl = TRUE)) return("LPI-P")
  if (grepl("^LPS\\s*\\(?O-|^LPS-O|^LPSo\\(", x, perl = TRUE)) return("LPS-O")
  if (grepl("^LPS\\s*\\(?P-|^LPS-P|^LPSp\\(", x, perl = TRUE)) return("LPS-P")
  if (grepl("^LPA\\s*\\(?O-|^LPA-O|^LPAo\\(", x, perl = TRUE)) return("LPA-O")
  if (grepl("^LPA\\s*\\(?P-|^LPA-P|^LPAp\\(", x, perl = TRUE)) return("LPA-P")
  # if (grepl("^LPC\\s*\\(?O-|^LPC-O|^LPCo\\(", x, perl = TRUE)) return("LPC")
  # if (grepl("^LPC\\s*\\(?P-|^LPC-P|^LPCp\\(", x, perl = TRUE)) return("LPC")
  # if (grepl("^LPE\\s*\\(?O-|^LPE-O|^LPEo\\(", x, perl = TRUE)) return("LPE")
  # if (grepl("^LPE\\s*\\(?P-|^LPE-P|^LPEp\\(", x, perl = TRUE)) return("LPE")
  # if (grepl("^LPG\\s*\\(?O-|^LPG-O|^LPGo\\(", x, perl = TRUE)) return("LPG")
  # if (grepl("^LPG\\s*\\(?P-|^LPG-P|^LPGp\\(", x, perl = TRUE)) return("LPG")
  # if (grepl("^LPI\\s*\\(?O-|^LPI-O|^LPIo\\(", x, perl = TRUE)) return("LPI")
  # if (grepl("^LPI\\s*\\(?P-|^LPI-P|^LPIp\\(", x, perl = TRUE)) return("LPI")
  # if (grepl("^LPS\\s*\\(?O-|^LPS-O|^LPSo\\(", x, perl = TRUE)) return("LPS")
  # if (grepl("^LPS\\s*\\(?P-|^LPS-P|^LPSp\\(", x, perl = TRUE)) return("LPS")
  # if (grepl("^LPA\\s*\\(?O-|^LPA-O|^LPAo\\(", x, perl = TRUE)) return("LPA")
  # if (grepl("^LPA\\s*\\(?P-|^LPA-P|^LPAp\\(", x, perl = TRUE)) return("LPA")

  # Standard lysophospholipids
  if (grepl("^LPA(\\(|\\s)", x, perl = TRUE)) return("LPA")
  if (grepl("^LPC(\\(|\\s)", x, perl = TRUE)) return("LPC")
  if (grepl("^LPE(\\(|\\s)", x, perl = TRUE)) return("LPE")
  if (grepl("^LPG(\\(|\\s)", x, perl = TRUE)) return("LPG")
  if (grepl("^LPI(\\(|\\s)", x, perl = TRUE)) return("LPI")
  if (grepl("^LPS(\\(|\\s)", x, perl = TRUE)) return("LPS")

  # Phospholipid ethers first
  if (grepl("^PC\\s*\\(?O-|^PC-O|^PCo\\(", x, perl = TRUE)) return("PC-O")
  if (grepl("^PC\\s*\\(?P-|^PC-P|^PCp\\(", x, perl = TRUE)) return("PC-P")
  if (grepl("^PE\\s*\\(?O-|^PE-O|^PEo\\(", x, perl = TRUE)) return("PE-O")
  if (grepl("^PE\\s*\\(?P-|^PE-P|^PEp\\(", x, perl = TRUE)) return("PE-P")
  if (grepl("^PG\\s*\\(?O-|^PG-O|^PGo\\(", x, perl = TRUE)) return("PG-O")
  if (grepl("^PG\\s*\\(?P-|^PG-P|^PGp\\(", x, perl = TRUE)) return("PG-P")
  if (grepl("^PI\\s*\\(?O-|^PI-O|^PIo\\(", x, perl = TRUE)) return("PI-O")
  if (grepl("^PI\\s*\\(?P-|^PI-P|^PIp\\(", x, perl = TRUE)) return("PI-P")
  if (grepl("^PS\\s*\\(?O-|^PS-O|^PSo\\(", x, perl = TRUE)) return("PS-O")
  if (grepl("^PS\\s*\\(?P-|^PS-P|^PSp\\(", x, perl = TRUE)) return("PS-P")
  if (grepl("^PA\\s*\\(?O-|^PA-O|^PAo\\(", x, perl = TRUE)) return("PA-O")
  if (grepl("^PA\\s*\\(?P-|^PA-P|^PAp\\(", x, perl = TRUE)) return("PA-P")
  # if (grepl("^PC\\s*\\(?O-|^PC-O|^PCo\\(", x, perl = TRUE)) return("PC")
  # if (grepl("^PC\\s*\\(?P-|^PC-P|^PCp\\(", x, perl = TRUE)) return("PC")
  # if (grepl("^PE\\s*\\(?O-|^PE-O|^PEo\\(", x, perl = TRUE)) return("PE")
  # if (grepl("^PE\\s*\\(?P-|^PE-P|^PEp\\(", x, perl = TRUE)) return("PE")
  # if (grepl("^PG\\s*\\(?O-|^PG-O|^PGo\\(", x, perl = TRUE)) return("PG")
  # if (grepl("^PG\\s*\\(?P-|^PG-P|^PGp\\(", x, perl = TRUE)) return("PG")
  # if (grepl("^PI\\s*\\(?O-|^PI-O|^PIo\\(", x, perl = TRUE)) return("PI")
  # if (grepl("^PI\\s*\\(?P-|^PI-P|^PIp\\(", x, perl = TRUE)) return("PI")
  # if (grepl("^PS\\s*\\(?O-|^PS-O|^PSo\\(", x, perl = TRUE)) return("PS")
  # if (grepl("^PS\\s*\\(?P-|^PS-P|^PSp\\(", x, perl = TRUE)) return("PS")
  # if (grepl("^PA\\s*\\(?O-|^PA-O|^PAo\\(", x, perl = TRUE)) return("PA")
  # if (grepl("^PA\\s*\\(?P-|^PA-P|^PAp\\(", x, perl = TRUE)) return("PA")

  # Standard phospholipids
  if (grepl("^PC(\\(|\\s)", x, perl = TRUE)) return("PC")
  if (grepl("^PE(\\(|\\s)", x, perl = TRUE)) return("PE")
  if (grepl("^PG(\\(|\\s)", x, perl = TRUE)) return("PG")
  if (grepl("^PI(\\(|\\s)", x, perl = TRUE)) return("PI")
  if (grepl("^PS(\\(|\\s)", x, perl = TRUE)) return("PS")
  if (grepl("^PA(\\(|\\s)", x, perl = TRUE)) return("PA")
  if (grepl("^BMP(\\(|\\s)", x, perl = TRUE)) return("BMP")
  if (grepl("^LNAPE(\\(|\\s)|^LNAOE(\\(|\\s)", x, perl = TRUE)) return("LNAPE")

  # Neutral lipids and glycoglycerolipids
  if (grepl("^MGDG(\\(|\\s)", x, perl = TRUE)) return("MGDG")
  if (grepl("^DGDG(\\(|\\s)", x, perl = TRUE)) return("DGDG")
  if (grepl("^MG(\\(|\\s)", x, perl = TRUE)) return("MG")
  if (grepl("^DG(\\(|\\s)", x, perl = TRUE)) return("DG")
  if (grepl("^TG(\\(|\\s)", x, perl = TRUE)) return("TG")
  if (grepl("^CL(\\(|\\s)", x, perl = TRUE)) return("CL")

  # Others
  if (grepl("^Carnitine(\\(|\\s)", x, perl = TRUE)) return("Carnitine")
  if (grepl("^CE(\\(|\\s)|[Cc]holesterol ester$", x, perl = TRUE)) return("CE")

  return("")
}

# guess_lipid_info
#' Guess lipid class and subclass from an annotation
#'
#' This helper returns both the base lipid class and the lipid subclass. For
#' most classes both values are identical, while for selected sphingolipids the
#' subclass retains d/t information.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return Named character vector with \code{lipidclass} and
#' \code{lipidsubclass}.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
guess_lipid_info <- function(x) {

  x <- norm(x)

  if (x == "" || is.na(x)) {
    return(c(lipidclass = "", lipidsubclass = ""))
  }

  if (is_deuterated_annotation(x)) {
    return(c(lipidclass = "", lipidsubclass = ""))
  }

  if (is_oxidized_non_sphingo_pl(x)) {
    return(c(lipidclass = "", lipidsubclass = ""))
  }

  # Reject unsupported ether glycoglycerolipids and neutral glycerolipids
  if (grepl("^(MGDG|DGDG|MG|DG|TG)(\\(|\\s)", x, perl = TRUE) &&
      grepl("(O-|P-|O_|P_)", x, perl = TRUE)) {
    return(c(lipidclass = "", lipidsubclass = ""))
  }

  sph_subclass <- guess_sphingo_subclass(x)
  if (sph_subclass != "") {
    return(c(lipidclass = base_lipidclass(sph_subclass),
             lipidsubclass = sph_subclass))
  }

  lipidclass <- guess_lipid_class(x)

  if (lipidclass == "") {
    return(c(lipidclass = "", lipidsubclass = ""))
  }

  # Handle ether phospholipid subclasses:
  # PX-O, PX-P, LPX-O, LPX-P, here X is any uppercase letter
  if (grepl("^(LP[A-Z]|P[A-Z])-(O|P)$", lipidclass, perl = TRUE)) {
    lipidsubclass <- lipidclass
    lipidclass <- sub("-(O|P)$", "", lipidclass)

    return(c(
      lipidclass = lipidclass,
      lipidsubclass = lipidsubclass
    ))
  }

  return(c(lipidclass = lipidclass, lipidsubclass = lipidclass))
}

# checkclass
#' Identify lipid classes from annotations
#'
#' This function checks lipid annotations and identifies whether they contain
#' known lipid classes. It can either return a general indication of presence
#' or a data frame with lipid class and lipid subclass for each annotation.
#'
#' @param ids Vector containing the original annotations to be converted into a
#' standardized format.
#' @param general Logical. If \code{TRUE}, the function checks whether any of
#' the lipid classes are present in the input annotations and returns
#' \code{"Yes"} if at least one match is found. If \code{FALSE}, the function
#' returns a data frame with \code{lipidclass} and \code{lipidsubclass}.
#'
#' @return \code{"Yes"} or \code{""} if \code{general = TRUE}. Otherwise, a
#' data frame with \code{lipidclass} and \code{lipidsubclass}.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkclass <- function(ids, general = TRUE) {

  ids <- gsub("no MS2: |low score: ", "", ids)
  lipids <- unlist(strsplit(ids, "\\||\\);"))
  lipids <- norm(lipids)

  info <- t(vapply(lipids, guess_lipid_info, character(2)))
  info <- as.data.frame(info, stringsAsFactors = FALSE)

  if (general) {
    if (any(info$lipidclass != "")) {
      return("Yes")
    } else {
      return("")
    }
  } else {
    return(info)
  }
}


# sphingo_to_lipidrt
#' Convert sphingolipid annotations to LipidRT format
#'
#' This helper converts sphingolipid annotations from external nomenclatures
#' into a LipidRT-compatible format.
#'
#' @param x Character string representing the lipid annotation to be converted.
#'
#' @return A character string with the converted lipid annotation.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
sphingo_to_lipidrt <- function(x) {

  x <- norm(x)

  x <- gsub("^CER\\b", "Cer", x, perl = TRUE)
  x <- gsub("^SPB\\b", "Sph", x, perl = TRUE)
  x <- gsub("^SPH\\b", "Sph", x, perl = TRUE)
  x <- gsub("^GlcCer\\b|^GalCer\\b|^Hex-Cer\\b", "HexCer", x, perl = TRUE)
  x <- gsub("^Sulfatide(?=\\s|\\()|^ST(?=\\s|\\()", "SHexCer", x, perl = TRUE)
  x <- gsub("^Cer1P\\b", "CerP", x, perl = TRUE)

  if (grepl("^(Cer|HexCer|AHexCer|SHexCer|CerP|SM)\\b", x, perl = TRUE) &&
      has_modified_fa_sphingo(x)) {
    return("")
  }

  if (invalid_sphingo_annotation(x)) {
    return("")
  }

  # Keep already valid shorthand
  if (grepl("^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|Sph)\\([dt]\\d+:\\d+(?:[/_\\-]\\d+:\\d+)?\\)$",
            x, perl = TRUE)) {
    return(x)
  }

  # MS-DIAL class tags
  x <- gsub("^(Cer)-(NS|NDS|AS|ADS|HDS|EOS)\\s+", "Cer d", x, perl = TRUE)
  x <- gsub("^(Cer)-(NP|AP|BDS|BP)\\s+", "Cer t", x, perl = TRUE)

  x <- gsub("^(HexCer)-(NS|NDS|AS|ADS|HDS|EOS)\\s+", "HexCer d", x, perl = TRUE)
  x <- gsub("^(HexCer)-(NP|AP|BDS|BP)\\s+", "HexCer t", x, perl = TRUE)

  x <- gsub("^(AHexCer)-(NS|NDS|AS|ADS|HDS|EOS)\\s+", "AHexCer d", x, perl = TRUE)
  x <- gsub("^(AHexCer)-(NP|AP|BDS|BP)\\s+", "AHexCer t", x, perl = TRUE)

  x <- gsub("^(SHexCer)-(NS|NDS|AS|ADS|HDS|EOS)\\s+", "SHexCer d", x, perl = TRUE)
  x <- gsub("^(SHexCer)-(NP|AP|BDS|BP)\\s+", "SHexCer t", x, perl = TRUE)

  # Total composition O2 / O3
  x <- gsub(
    "^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|Sph|SphP|AcylCer)\\s+(\\d+:\\d+);(?:O2|2O)$",
    "\\1(d\\2)", x, perl = TRUE)

  x <- gsub(
    "^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|Sph|SphP|AcylCer)\\s+(\\d+:\\d+);(?:O3|3O)$",
    "\\1(t\\2)", x, perl = TRUE)

  # Molecular species O2 / O3
  x <- gsub(
    "^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|AcylCer)\\s+(\\d+:\\d+);(?:O2|2O)/(\\d+:\\d+)$",
    "\\1(d\\2/\\3)", x, perl = TRUE)

  x <- gsub(
    "^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|AcylCer)\\s+(\\d+:\\d+);(?:O3|3O)/(\\d+:\\d+)$",
    "\\1(t\\2/\\3)", x, perl = TRUE)

  # Shorthand without parentheses
  x <- gsub("^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|AcylCer)\\s+([dt]\\d+:\\d+[/_\\-]\\d+:\\d+(?:;O\\d+)?)$",
            "\\1(\\2)", x, perl = TRUE)
  x <- gsub("^Sph\\s+([dt]\\d+:\\d+)$", "Sph(\\1)", x, perl = TRUE)

  # Convert sphingolipid shorthand notation to LipidRT format
  x <- gsub(
    "^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|Sph|SphP|AcylCer|)\\s+([dt]?)(\\d+:\\d+)(.*)$",
    "\\1(\\2\\3\\4)",
    x,
    perl = TRUE
  )

  # Add default d sphingoid base when missing
  x <- gsub(
    "^(Cer|HexCer|AHexCer|SHexCer|SM|CerP|Sph)\\((\\d+:\\d+)(.*)\\)$",
    "\\1(d\\2\\3)",
    x,
    perl = TRUE
  )

  return(x)
}

# sum_total_chains
#' Sum total carbons and double bonds from chain annotations
#'
#' This helper extracts all chain annotations in the form \code{C:DB} and
#' returns the total composition.
#'
#' @param x Character string representing the lipid annotation.
#'
#' @return Character string with total composition \code{C:DB}.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
sum_total_chains <- function(x) {

  x <- norm(x)

  x <- gsub("^[A-Za-z0-9_\\-]+\\(", "", x, perl = TRUE)
  x <- gsub("\\)$", "", x, perl = TRUE)

  x <- gsub("N-", "", x, fixed = TRUE)
  x <- gsub("^O-", "", x, perl = TRUE)
  x <- gsub("^P-", "", x, perl = TRUE)

  chains <- regmatches(x, gregexpr("\\d+:\\d+", x, perl = TRUE))[[1]]

  if (length(chains) == 0) {
    return("")
  }

  carbons <- sum(as.numeric(sub(":.*", "", chains)))
  db <- sum(as.numeric(sub(".*:", "", chains)))

  return(paste0(carbons, ":", db))
}

# convertlipidannotation
#' Convert lipid annotations between formats
#'
#' This function transforms lipid annotation strings from one format to another
#' (e.g., from MS-DIAL or external sources to LipidRT format), applying
#' specific rules based on lipid class and source format.
#'
#' @param lipid Character string representing the lipid annotation to be
#' converted.
#' @param lipidclass Character string indicating the lipid base class.
#' @param lipidsubclass Character string indicating the lipid subclass.
#'
#' @return A character string with the converted lipid annotation.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
convertlipidannotation <- function(lipid, lipidclass, lipidsubclass = lipidclass) {

  lipid <- norm(lipid)

  if (is_deuterated_annotation(lipid)) {
    return("")
  }

  if (is_oxidized_non_sphingo_pl(lipid)) {
    return("")
  }

  # Reject unsupported ether glycoglycerolipids and neutral glycerolipids
  if (grepl("^(MGDG|DGDG|MG|DG|TG)(\\(|\\s)", lipid, perl = TRUE) &&
      grepl("(O-|P-|O_|P_)", lipid, perl = TRUE)) {
    return("")
  }

  lipid <- gsub("^CAR\\s+|^ACAR\\s+|^ACar\\s+|^Car\\s+", "Carnitine ", lipid, perl = TRUE)
  lipid <- gsub("^(\\d+:\\d+)\\s+Carnitine$", "Carnitine \\1", lipid, perl = TRUE)

  lipid <- gsub("^GlcCer\\s+|^GalCer\\s+|^Hex-Cer\\s+", "HexCer ", lipid, perl = TRUE)
  lipid <- gsub("^Sulfatide\\s+|^ST\\s+", "SHexCer ", lipid, perl = TRUE)
  lipid <- gsub("^Cer1P\\s+", "CerP ", lipid, perl = TRUE)
  lipid <- gsub("^ACER\\s+|^ACer\\s+", "AcylCer ", lipid, perl = TRUE)
  lipid <- gsub("^LNAOE\\s+", "LNAPE ", lipid, perl = TRUE)

  if (lipidclass %in% c("Cer", "HexCer", "AHexCer", "SHexCer",
                        "CerP", "SM", "Sph", "SphP", "AcylCer")) {

    lipid <- gsub("^CER\\s+P-\\s*", "CerP ", lipid, perl = TRUE)
    lipid <- gsub("^Cer_(?:NS|NDS|NP|AS|ADS|AP|HDS|EOS|BDS|BP)\\s+", "Cer ", lipid, perl = TRUE)
    lipid <- gsub("^(CER|Cer)\\s+", "Cer ", lipid, perl = TRUE)
    lipid <- gsub("^(SPB|SPH)\\s+", "Sph ", lipid, perl = TRUE)

    lipid <- sphingo_to_lipidrt(lipid)

    if (lipid == "") {
      return("")
    }
  }

  lipid <- gsub("^FA\\s+", "FA(", lipid, perl = TRUE)
  lipid <- gsub("^HFA\\s+", "HFA(", lipid, perl = TRUE)
  lipid <- gsub("^FAHFA\\s+", "FAHFA(", lipid, perl = TRUE)
  lipid <- gsub("^NAGlySer\\s+", "NAGlySer(", lipid, perl = TRUE)
  lipid <- gsub("^NAGly\\s+", "NAGly(", lipid, perl = TRUE)
  lipid <- gsub("^NAOrn\\s+", "NAOrn(", lipid, perl = TRUE)
  lipid <- gsub("^NAE\\s+", "NAE(", lipid, perl = TRUE)
  lipid <- gsub("^AcylCer\\s+", "AcylCer(", lipid, perl = TRUE)

  lipid <- gsub("^PC(?:\\s+O-|-O(?:\\s+|\\())", "PC(O-", lipid, perl = TRUE)
  lipid <- gsub("^PC(?:\\s+P-|-P(?:\\s+|\\())", "PC(P-", lipid, perl = TRUE)
  lipid <- gsub("^PE(?:\\s+O-|-O(?:\\s+|\\())", "PE(O-", lipid, perl = TRUE)
  lipid <- gsub("^PE(?:\\s+P-|-P(?:\\s+|\\())", "PE(P-", lipid, perl = TRUE)
  lipid <- gsub("^PG(?:\\s+O-|-O(?:\\s+|\\())", "PG(O-", lipid, perl = TRUE)
  lipid <- gsub("^PG(?:\\s+P-|-P(?:\\s+|\\())", "PG(P-", lipid, perl = TRUE)
  lipid <- gsub("^PI(?:\\s+O-|-O(?:\\s+|\\())", "PI(O-", lipid, perl = TRUE)
  lipid <- gsub("^PI(?:\\s+P-|-P(?:\\s+|\\())", "PI(P-", lipid, perl = TRUE)
  lipid <- gsub("^PS(?:\\s+O-|-O(?:\\s+|\\())", "PS(O-", lipid, perl = TRUE)
  lipid <- gsub("^PS(?:\\s+P-|-P(?:\\s+|\\())", "PS(P-", lipid, perl = TRUE)
  lipid <- gsub("^PA(?:\\s+O-|-O(?:\\s+|\\())", "PA(O-", lipid, perl = TRUE)
  lipid <- gsub("^PA(?:\\s+P-|-P(?:\\s+|\\())", "PA(P-", lipid, perl = TRUE)

  lipid <- gsub("^LPC(?:\\s+O-|-O(?:\\s+|\\())", "LPC(O-", lipid, perl = TRUE)
  lipid <- gsub("^LPC(?:\\s+P-|-P(?:\\s+|\\())", "LPC(P-", lipid, perl = TRUE)
  lipid <- gsub("^LPE(?:\\s+O-|-O(?:\\s+|\\())", "LPE(O-", lipid, perl = TRUE)
  lipid <- gsub("^LPE(?:\\s+P-|-P(?:\\s+|\\())", "LPE(P-", lipid, perl = TRUE)
  lipid <- gsub("^LPG(?:\\s+O-|-O(?:\\s+|\\())", "LPG(O-", lipid, perl = TRUE)
  lipid <- gsub("^LPG(?:\\s+P-|-P(?:\\s+|\\())", "LPG(P-", lipid, perl = TRUE)
  lipid <- gsub("^LPI(?:\\s+O-|-O(?:\\s+|\\())", "LPI(O-", lipid, perl = TRUE)
  lipid <- gsub("^LPI(?:\\s+P-|-P(?:\\s+|\\())", "LPI(P-", lipid, perl = TRUE)
  lipid <- gsub("^LPS(?:\\s+O-|-O(?:\\s+|\\())", "LPS(O-", lipid, perl = TRUE)
  lipid <- gsub("^LPS(?:\\s+P-|-P(?:\\s+|\\())", "LPS(P-", lipid, perl = TRUE)
  lipid <- gsub("^LPA(?:\\s+O-|-O(?:\\s+|\\())", "LPA(O-", lipid, perl = TRUE)
  lipid <- gsub("^LPA(?:\\s+P-|-P(?:\\s+|\\())", "LPA(P-", lipid, perl = TRUE)

  # lipid <- gsub("^LPC\\s+O-", "LPC-O(", lipid, perl = TRUE)
  # lipid <- gsub("^LPC\\s+P-", "LPC-P(", lipid, perl = TRUE)
  # lipid <- gsub("^LPE\\s+O-", "LPE-O(", lipid, perl = TRUE)
  # lipid <- gsub("^LPE\\s+P-", "LPE-P(", lipid, perl = TRUE)
  # lipid <- gsub("^LPG\\s+O-", "LPG-O(", lipid, perl = TRUE)
  # lipid <- gsub("^LPG\\s+P-", "LPG-P(", lipid, perl = TRUE)
  # lipid <- gsub("^LPI\\s+O-", "LPI-O(", lipid, perl = TRUE)
  # lipid <- gsub("^LPI\\s+P-", "LPI-P(", lipid, perl = TRUE)
  # lipid <- gsub("^LPS\\s+O-", "LPS-O(", lipid, perl = TRUE)
  # lipid <- gsub("^LPS\\s+P-", "LPS-P(", lipid, perl = TRUE)
  # lipid <- gsub("^LPA\\s+O-", "LPA-O(", lipid, perl = TRUE)
  # lipid <- gsub("^LPA\\s+P-", "LPA-P(", lipid, perl = TRUE)

  lipid <- gsub("^PC\\s+", "PC(", lipid, perl = TRUE)
  lipid <- gsub("^PE\\s+", "PE(", lipid, perl = TRUE)
  lipid <- gsub("^PG\\s+", "PG(", lipid, perl = TRUE)
  lipid <- gsub("^PI\\s+", "PI(", lipid, perl = TRUE)
  lipid <- gsub("^PS\\s+", "PS(", lipid, perl = TRUE)
  lipid <- gsub("^PA\\s+", "PA(", lipid, perl = TRUE)

  lipid <- gsub("^LPA\\s+", "LPA(", lipid, perl = TRUE)
  lipid <- gsub("^LPC\\s+", "LPC(", lipid, perl = TRUE)
  lipid <- gsub("^LPE\\s+", "LPE(", lipid, perl = TRUE)
  lipid <- gsub("^LPG\\s+", "LPG(", lipid, perl = TRUE)
  lipid <- gsub("^LPI\\s+", "LPI(", lipid, perl = TRUE)
  lipid <- gsub("^LPS\\s+", "LPS(", lipid, perl = TRUE)

  lipid <- gsub("^MGDG\\s+", "MGDG(", lipid, perl = TRUE)
  lipid <- gsub("^DGDG\\s+", "DGDG(", lipid, perl = TRUE)
  lipid <- gsub("^BMP\\s+", "BMP(", lipid, perl = TRUE)
  lipid <- gsub("^LNAPE\\s+", "LNAPE(", lipid, perl = TRUE)

  lipid <- gsub("^CL\\s+", "CL(", lipid, perl = TRUE)
  lipid <- gsub("^MG\\s+", "MG(", lipid, perl = TRUE)
  lipid <- gsub("^DG\\s+", "DG(", lipid, perl = TRUE)
  lipid <- gsub("^TG\\s+", "TG(", lipid, perl = TRUE)
  lipid <- gsub("^Carnitine\\s+", "Carnitine(", lipid, perl = TRUE)
  lipid <- gsub("^CE\\s+", "CE(", lipid, perl = TRUE)

  if (lipidclass == "LNAPE") {
    lipid <- gsub("^LNAPE\\s+", "LNAPE(", lipid, perl = TRUE)
    if (!grepl("\\)$", lipid)) {
      lipid <- paste0(lipid, ")")
    }
    return(lipid)
  }

  if (!(lipidclass %in% c("Cer", "HexCer", "AHexCer", "SHexCer",
                          "CerP", "SM", "Sph", "SphP", "AcylCer"))) {
    lipid <- gsub(";O\\d+|;O", "", lipid, perl = TRUE)
  }

  if (!grepl("\\)$", lipid)) {
    lipid <- paste0(lipid, ")")
  }

  return(lipid)
}

# organizertdata
#' Update annotated peaklist with the filtered data
#'
#' Update annotated peaklist with the filtered data
#'
#' @param rtdata rtdata object.
#' @param typeanalysis Character. Specifies whether the function is being used for
#' filtering or RT predictions. Possible values are "filtering" or "prediction".
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
organizertdata <- function(rtdata, typeanalysis) {

  typeanalysis <- match.arg(typeanalysis, choices = c("filtering", "prediction"))

  if (typeanalysis == "filtering") {
    # Modify the filter tag generated during the filtering process
    rtdata$filtering$data <- checkfiltervalue(data = rtdata$filtering$data)

    rtdata <- groupidfilter(rtdata = rtdata)

    # Info indicating how many lipids have been removed during filtering
    rtdata$filtering$information <- generateinfo(rtdata = rtdata,
                                                 typeanalysis = typeanalysis)

  } else if (typeanalysis == "prediction") {
    # Update annotated peaklist with the predicted data
    rtdata <- groupidrtpred(rtdata = rtdata)

    # Print a message indicating how many lipids have been predicted
    rtdata$prediction$information <- generateinfo(rtdata = rtdata,
                                                  typeanalysis = typeanalysis)
  }

  return(rtdata)
}

# removepreviousanalysis
#' If the data had already been filtered and/or predicted, remove the previous
#' results
#'
#' # If the data had already been filtered and/or predicted, remove the previous
#' results. his is necessary because both for 'msobject' and 'msbatch', there
#' would be problems if filtering or predictions were performed again on the same
#' data when previous results already exist.
#'
#' @param rtdata rtdata object.
#' @param typeanalysis Character. Specifies whether the function is being used for
#' filtering or RT predictions. Possible values are "filtering" or "prediction".
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
removepreviousanalysis <- function(rtdata, typeanalysis) {

  if (typeanalysis == "filtering") {
    # If filtering results already exist
    if ("lipidname_filter" %in% colnames(rtdata$features)) {
      # If prediction results doesn't exist only remove filtering-related columns
      if (!("lipidname_rtpred" %in% colnames(rtdata$features))) {
        rtdata$features <- rtdata$features[, !colnames(rtdata$features) %in%
                                                       c("filter_label","lipidname_filter",
                                                         "LipidClass_filter", "LipidSubclass_filter",
                                                         "C_filter", "DB_filter", "adduct_filter",
                                                         "confidenceLevel_filter", "Score_filter",
                                                         "GroupID")]

      } else {
        # If prediction results exist remove both filtering and prediction columns
        rtdata$features <- rtdata$features[, !colnames(rtdata$features) %in%
                                                       c("filter_label","lipidname_filter",
                                                         "LipidClass_filter", "LipidSubclass_filter",
                                                         "C_filter", "DB_filter", "adduct_filter",
                                                         "confidenceLevel_filter", "Score_filter",
                                                         "lipidname_rtpred", "LipidClass_rtpred",
                                                         "LipidSubclass_rtpred", "C_rtpred", "DB_rtpred",
                                                         "adduct_rtpred", "confidenceLevel_rtpred",
                                                         "RT_rtpred", "Score_rtpred",
                                                         "lipidname_final", "LipidClass_final",
                                                         "LipidSubclass_final", "C_final", "DB_final",
                                                         "adduct_final", "confidenceLevel_final",
                                                         "Score_final", "GroupID")]
      }
      # Remove the filtering object entirely
      rtdata$filtering <- list()
      rtdata$prediction <- list()
    }

  } else if (typeanalysis == "prediction") {
    # If prediction results already exist
    if ("lipidname_rtpred" %in% colnames(rtdata$features)) {
      # Reset prediction to an empty list
      rtdata$prediction <- list()

      # Remove prediction related columns
      rtdata$features <- rtdata$features[, !colnames(rtdata$features) %in%
                                                     c("lipidname_rtpred", "LipidClass_rtpred",
                                                       "LipidSubclass_rtpred", "C_rtpred",
                                                       "DB_rtpred", "adduct_rtpred",
                                                       "confidenceLevel_rtpred",
                                                       "RT_rtpred", "Score_rtpred",
                                                       "lipidname_final", "LipidClass_final",
                                                       "LipidSubclass_final", "C_final", "DB_final",
                                                       "adduct_final", "confidenceLevel_final",
                                                       "Score_final", "GroupID")]

      # If filtering results still exist, regroup features accordingly
      if ("lipidname_filter" %in% colnames(rtdata$features)) {
        rtdata <- groupidfilter(rtdata = rtdata)
      }
    }
  }

  return(rtdata)
}

# checkdatastructure
#' Checks the structure of the provided data to determine whether it is an
#' msbatch or an msobject.
#'
#' Checks the structure of the provided data to determine whether it is an
#' msbatch or an msobject. If the format is incorrect, it provides a message
#' explaining why it does not meet the required structure.
#'
#' @param rtdata rtdata object
#'
#' @return rtdata_structure list (structure, message)
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkdatastructure <- function(rtdata) {
  rtdata_structure <- "unknown"
  error_message <- NULL

  # Check if it is a list
  if (!is.list(rtdata)) {
    return(list(structure = rtdata_structure, message = "Error: The object is not a list."))
  }

  if (all(c("metadata", "features", "data","filtering", "prediction") %in% names(rtdata)) &&
      is.list(rtdata$metadata) &&
      is.data.frame(rtdata$features) &&
      is.data.frame(rtdata$data) &&
      is.list(rtdata$filtering) &&
      is.list(rtdata$prediction)) {
    rtdata_structure <- "rtdata"
    error_message <- NULL
  } else {
    error_message <- "Error: The imported object does not match the expected structure.
    Please make sure you are providing a valid 'rtdata' object."
  }

  return(list(structure = rtdata_structure, message = error_message))
}

# generateinfo
#' Provide information about lipid filtering and predictions.
#'
#' Provide information about lipid filtering and predictions. When filtering data
#' is provided <filtering$data>, the function provides detailed
#' information about the lipids processed in each class. Specifically, it outlines
#' how many lipids have been kept, how many have been recovered by the model, and
#' how many have been eliminated from each class. Additionally, it includes the
#' minimum, maximum, and average values of ScoreRT_pred for the lipids in each
#' class after performing the filtering.
#' When, instead of filtering data, prediction data <candidates> is provided, the
#' function generates an output similar to the filtering data, but with additional
#' information about the predictions. This includes the minimum, maximum, and
#' average values of ScoreRT_pred for the lipids in each class after making
#' predictions for the candidates, as well as the number of new annotations
#' proposed for features that had no previous annotation, based on the score
#' provided by the user.
#'
#' @param rtdata rtdata
#' @param typeanalysis character. Specifies whether the function is being used for
#' filtering or RT predictions. Possible values are "filtering" or "prediction".
#'
#' @return updated_info
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
generateinfo <- function(rtdata, typeanalysis = c("filtering", "prediction")) {

  typeanalysis <- match.arg(typeanalysis)

  ## ======================================================
  ## FILTERING
  ## ======================================================
  if (typeanalysis == "filtering") {

    data <- rtdata$filtering$data

    if (is.null(data) || nrow(data) == 0) {
      return(data.frame())
    }

    classes <- unique(data$LipidSubclass[!is.na(data$LipidSubclass) & data$LipidSubclass != ""])

    keep_summary <- tapply(
      data$checkfilter %in% c("Kept", "Kept manually by the user"),
      factor(data$LipidSubclass, levels = classes),
      sum, na.rm = TRUE
    )

    recover_summary <- tapply(
      data$checkfilter %in% c("Recovered with RT modelling"),
      factor(data$LipidSubclass, levels = classes),
      sum, na.rm = TRUE
    )

    remove_summary <- tapply(
      data$checkfilter %in% c(
        "In-source fragments", "Removed distant data",
        "Confidence level among isobars", "Removed by RT rules",
        "Removed manually by the user", "Outlier",
        "RT coherence between subclasses"
      ),
      factor(data$LipidSubclass, levels = classes),
      sum, na.rm = TRUE
    )

    ignore_summary <- tapply(
      data$checkfilter %in% c("Not pass QC"),
      factor(data$LipidSubclass, levels = classes),
      sum, na.rm = TRUE
    )

    data_filtermodel <- data[data$Filter %in% "Kept", , drop = FALSE]
    data_keep <- data[data$Filter %in% c(
      "Kept", "Kept manually by the user", "Recovered with RT modelling"
    ), , drop = FALSE]
    data_remove <- data[data$Filter %in% c(
      "In-source fragments", "Removed distant data",
      "Confidence level among isobars", "Removed by RT rules",
      "Removed manually by the user", "Outlier",
      "RT coherence between subclasses"
    ), , drop = FALSE]

    # Scores
    if ("Score_rtpred" %in% names(data)) {

      num_score <- as.numeric(data$Score_rtpred)

      score_filtermodel <- tapply(
        as.numeric(data_filtermodel$Score_rtpred),
        factor(data_filtermodel$LipidSubclass, levels = classes),
        mean, na.rm = TRUE
      )

      min_filtermodel <- tapply(
        as.numeric(data_filtermodel$Score_rtpred),
        factor(data_filtermodel$LipidSubclass, levels = classes),
        min, na.rm = TRUE
      )

      max_filtermodel <- tapply(
        as.numeric(data_filtermodel$Score_rtpred),
        factor(data_filtermodel$LipidSubclass, levels = classes),
        max, na.rm = TRUE
      )

      min_score_keep <- tapply(
        as.numeric(data_keep$Score_rtpred),
        factor(data_keep$LipidSubclass, levels = classes),
        min, na.rm = TRUE
      )

      mean_score_keep <- tapply(
        as.numeric(data_keep$Score_rtpred),
        factor(data_keep$LipidSubclass, levels = classes),
        mean, na.rm = TRUE
      )

      max_score_keep <- tapply(
        as.numeric(data_keep$Score_rtpred),
        factor(data_keep$LipidSubclass, levels = classes),
        max, na.rm = TRUE
      )

      min_score_remove <- tapply(
        as.numeric(data_remove$Score_rtpred),
        factor(data_remove$LipidSubclass, levels = classes),
        min, na.rm = TRUE
      )

      mean_score_remove <- tapply(
        as.numeric(data_remove$Score_rtpred),
        factor(data_remove$LipidSubclass, levels = classes),
        mean, na.rm = TRUE
      )

      max_score_remove <- tapply(
        as.numeric(data_remove$Score_rtpred),
        factor(data_remove$LipidSubclass, levels = classes),
        max, na.rm = TRUE
      )

    } else {
      score_filtermodel <- min_filtermodel <- max_filtermodel <- NA
      min_score_keep <- mean_score_keep <- max_score_keep <- NA
      min_score_remove <- mean_score_remove <- max_score_remove <- NA
    }

    updated_info <- data.frame(
      LipidSubclass = classes,
      keep = keep_summary,
      recover = recover_summary,
      remove = remove_summary,
      ignored = ignore_summary,

      min_filtermodel = round(min_filtermodel, 3),
      score_filtermodel = round(score_filtermodel, 3),
      max_filtermodel = round(max_filtermodel, 3),

      min_score_keep = round(min_score_keep, 3),
      mean_score_keep = round(mean_score_keep, 3),
      max_score_keep = round(max_score_keep, 3),

      min_score_remove = round(min_score_remove, 3),
      mean_score_remove = round(mean_score_remove, 3),
      max_score_remove = round(max_score_remove, 3),

      row.names = NULL
    )

    return(updated_info)
  }

  ## ======================================================
  ## PREDICTION
  ## ======================================================
  if (typeanalysis == "prediction") {

    classes <- names(rtdata$prediction$models)

    if (length(classes) == 0) {
      return(data.frame())
    }

    data <- rtdata$filtering$data
    candidates <- rtdata$prediction$candidates

    data_predictmodel <- tapply(
      data$checkfilter %in% c("Kept", "Recovered with RT modelling", "Kept manually by the user"),
      factor(data$LipidSubclass, levels = classes),
      sum, na.rm = TRUE
    )

    if (is.null(candidates) || nrow(candidates) == 0) {

      return(data.frame(
        LipidSubclass = NA,
        new_rtannotations = 0,
        discarded_rtannotations = 0,
        data_predictmodel = NA,

        min_predictmodel = NA,
        score_predictmodel = NA,
        max_predictmodel = NA,

        min_score_keep = NA,
        mean_score_keep = NA,
        max_score_keep = NA,
        min_score_discarded = NA,
        mean_score_discarded = NA,
        max_score_discarded = NA,
        min_score_remove = NA,
        mean_score_remove = NA,
        max_score_remove = NA,

        row.names = NULL
      ))
    }

    score_num <- as.numeric(candidates$Score_rtpred)

    score_predictmodel <- tapply(
      score_num,
      factor(candidates$LipidSubclass, levels = classes),
      mean, na.rm = TRUE
    )

    min_predictmodel <- tapply(
      score_num,
      factor(candidates$LipidSubclass, levels = classes),
      min, na.rm = TRUE
    )

    max_predictmodel <- tapply(
      score_num,
      factor(candidates$LipidSubclass, levels = classes),
      max, na.rm = TRUE
    )

    new_rtannotations <- tapply(
      candidates$checkcandidate %in% c("Kept", "Kept manually by the user"),
      factor(candidates$LipidSubclass, levels = classes),
      sum, na.rm = TRUE
    )

    discarded_rtannotations <- tapply(
      candidates$checkcandidate %in% c("Discarded", "Removed manually by the user"),
      factor(candidates$LipidSubclass, levels = classes),
      sum, na.rm = TRUE
    )

    candidates_keep <- candidates[candidates$checkcandidate %in%
                                    c("Kept", "Kept manually by the user"), , drop = FALSE]
    candidates_discarded <- candidates[candidates$checkcandidate %in%
                                         "Discarded", , drop = FALSE]
    candidates_remove <- candidates[candidates$checkcandidate %in%
                                      "Removed", , drop = FALSE]

    agg <- function(df, fun) {
      if (nrow(df) == 0) return(rep(NA, length(classes)))
      tapply(as.numeric(df$Score_rtpred),
             factor(df$LipidSubclass, levels = classes),
             fun, na.rm = TRUE)
    }

    updated_info <- data.frame(
      LipidSubclass = classes,
      new_rtannotations = new_rtannotations,
      discarded_rtannotations = discarded_rtannotations,
      data_predictmodel = data_predictmodel,

      min_predictmodel = round(min_predictmodel, 3),
      score_predictmodel = round(score_predictmodel, 3),
      max_predictmodel = round(max_predictmodel, 3),

      min_score_keep = round(agg(candidates_keep, min), 3),
      mean_score_keep = round(agg(candidates_keep, mean), 3),
      max_score_keep = round(agg(candidates_keep, max), 3),

      min_score_discarded = round(agg(candidates_discarded, min), 3),
      mean_score_discarded = round(agg(candidates_discarded, mean), 3),
      max_score_discarded = round(agg(candidates_discarded, max), 3),

      min_score_remove = round(agg(candidates_remove, min), 3),
      mean_score_remove = round(agg(candidates_remove, mean), 3),
      max_score_remove = round(agg(candidates_remove, max), 3),

      row.names = NULL
    )

    return(updated_info)
  }
}


# determinepolarity
#' This function automatically determines the polarity ("positive", "negative").
#'
#' This function automatically determines the polarity ("positive", "negative").
#' It checks the adducts against a reference list of known positive and negative
#' adducts for lipid classes.
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#'
#' @return polarity
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
determinepolarity <- function(data) {

  adductsTable <- LipidRT::adductsTable

    # Load default adducts
    adducts_pos <- adductsTable$adduct[adductsTable$charge > 0]
    adducts_neg <- adductsTable$adduct[adductsTable$charge < 0]

    adducts_pos <- adducts_pos[!is.na(adducts_pos) & adducts_pos != ""]
    adducts_neg <- adducts_neg[!is.na(adducts_neg) & adducts_neg != ""]


    # Get all adducts from the features and split multiple adducts separated by "|"
    adducts <- unique(data$adduct)

    adducts_valid <- adducts[adducts %in% c(adducts_pos, adducts_neg)]

    is_pos <- any(adducts_valid %in% adducts_pos)
    is_neg <- any(adducts_valid %in% adducts_neg)

    if (is_pos && !is_neg) {
      polarity <- "positive"
    } else if (!is_pos & is_neg) {
      polarity <- "negative"
    } else {
      if (length(adducts) == 0 && adducts != "") { # Target
        polarity <- NULL
        warning(
          paste("Adducts must correspond to either positive or negative polarity only.\n",
                "Please ensure all adducts belong to a single polarity group.\n")
        )
      } else {
        polarity <- "both" # both
        warning(
          paste("Adducts must correspond to either positive or negative polarity only.\n",
                "Please ensure all adducts belong to a single polarity group.\n")
        )
      }
    }

  return(polarity)
}


# predictionlayergraph
#' Generate RT prediction lines for plotting
#'
#' Generate RT prediction lines for plotting
#'
#' @param rtdata rtdata object.
#' @param lipidsubclass Character. Lipid class to analyze.
#' @param min_c Numeric. Minimum carbon number.
#' @param max_c Numeric. Maximum carbon number.
#' @param min_rt Numeric. Minimum retention time.
#' @param max_rt Numeric. Maximum retention time.
#' @param alpha Numeric. Transparency level.
#' @param typeanalysis Character. Specifies whether the function is being used for
#' filtering or RT predictions. Possible values are "filtering" or "prediction".
#'
#' @return prediction_layer
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
predictionlayergraph <- function(rtdata, lipidsubclass, min_c, max_c, min_rt, max_rt, alpha, typeanalysis) {

  typeanalysis <- match.arg(typeanalysis, choices = c("filtering", "prediction"))

  if (typeanalysis == "filtering") {
    predictions <- rtdata$filtering$predictions
    rtpred_min <- rtdata$filtering$restrictrtpred$min
    rtpred_max <- rtdata$filtering$restrictrtpred$max
  } else if (typeanalysis == "prediction") {
    predictions <- rtdata$prediction$predictions
    rtpred_min <- rtdata$prediction$restrictrtpred$min
    rtpred_max <- rtdata$prediction$restrictrtpred$max
  }

  if (!is.null(predictions) && nrow(predictions) > 0) {
    predict_line <- predictions[predictions$LipidSubclass == lipidsubclass &
                                  predictions$C >= min_c &
                                  predictions$C <= max_c &
                                  predictions$RTsec_rtpred >= min_rt &
                                  predictions$RTsec_rtpred <= max_rt, c(1:5), drop = FALSE]

    # Filter prediction lines
    counts <- table(predict_line$DB)
    keep <- names(counts[counts >= 2])
    predict_line <- predict_line[predict_line$DB %in% keep, ]

    if (nrow(predict_line) > 0) {
      # -------------------------
      # Local minimum region
      # -------------------------
      if (lipidsubclass %in% unique(rtpred_min$class)) {
        c_min <- unique(rtpred_min$C_min[
          rtpred_min$class == lipidsubclass
        ])
        predict_line_min <- predict_line[predict_line$C <= c_min, ]
        predict_line <- predict_line[predict_line$C >= c_min, ]

        # Retain only DBs with at least two points in predict_line_min
        counts_min <- table(predict_line_min$DB)
        keep_min <- names(counts_min[counts_min >= 2])
        predict_line_min <- predict_line_min[predict_line_min$DB %in% keep_min, ]
      } else {
        predict_line_min <- data.frame()
      }

      # -------------------------
      # Local maximum region
      # -------------------------
      if (lipidsubclass %in% unique(rtpred_max$class)) {
        c_max <- unique(rtpred_max$C_max[
          rtpred_max$class == lipidsubclass
        ])
        predict_line_max <- predict_line[predict_line$C >= c_max, ]
        predict_line <- predict_line[predict_line$C <= c_max, ]

        # Retain only DBs with at least two points in predict_line_max
        counts_max <- table(predict_line_max$DB)
        keep_max <- names(counts_max[counts_max >= 2])
        predict_line_max <- predict_line_max[predict_line_max$DB %in% keep_max, ]
      } else {
        predict_line_max <- data.frame()
      }
    } else {
      predict_line <- data.frame()
      predict_line_min <- data.frame()
      predict_line_max <- data.frame()
    }

  } else {
    predict_line <- data.frame()
    predict_line_min <- data.frame()
    predict_line_max <- data.frame()
  }

  C <- DB <- RTsec_rtpred <- c()
  # check if prediction layer is required
  prediction_layer <- list()

  if (nrow(predict_line) >= 2) {
    prediction_layer <- list(
      geom_line(
        data = predict_line,
        aes(x = C, y = RTsec_rtpred, group = factor(DB), color = factor(DB)),
        linewidth = 0.7, linetype = "dashed", alpha = alpha, inherit.aes = FALSE
      )
    )

    # Add the minimum segment (gray) when at least two points are available
    if (exists("predict_line_min") && !is.null(predict_line_min) && nrow(predict_line_min) >= 2) {
      prediction_layer <- c(
        prediction_layer,
        list(
          geom_line(
            data = predict_line_min,
            aes(x = C, y = RTsec_rtpred, group = factor(DB)),
            color = "grey60", linewidth = 0.8, linetype = "dotted", alpha = alpha,
            inherit.aes = FALSE
          )
        )
      )
    }

    # Add the maximum segment (gray) when at least two points are available
    if (exists("predict_line_max") && !is.null(predict_line_max) && nrow(predict_line_max) >= 2) {
      prediction_layer <- c(
        prediction_layer,
        list(
          geom_line(
            data = predict_line_max,
            aes(x = C, y = RTsec_rtpred, group = factor(DB)),
            color = "grey60", linewidth = 0.8, linetype = "dotted", alpha = alpha,
            inherit.aes = FALSE
          )
        )
      )
    }
  }

  return(prediction_layer)
}




# collapse_or_empty
#' Collapse unique, non-empty values into a single character string.
#'
#' This function removes \code{NA} values and empty strings from the input
#' vector. If no valid values remain, it returns an empty string.
#' If one unique value remains, it returns that value.
#' If multiple unique values remain, it collapses them into a single string
#' separated by \code{sep}.
#'
#' @param x A character vector.
#' @param sep Character string used to separate values when collapsing.
#'
#' @return A character string: empty, single value, or collapsed unique values.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
collapse_or_empty <- function(x, sep) {
  x <- x[!is.na(x) & x != ""]
  ux <- unique(x)

  if (length(ux) == 0) {
    ""
  } else if (length(ux) > 1) {
    paste(ux, collapse = sep)
  } else {
    ux
  }
}

# collapse_conditional
#' Collapse conditional, non-empty values into a single character string.
#'
#' This function removes \code{NA} values and empty strings from the input
#' vector. If no valid values remain, it returns an empty string.
#' If one unique value remains, it returns that value.
#' If multiple unique values remain, it collapses them into a single string
#' separated by \code{sep}.
#'
#' @param x A character vector.
#' @param sep Character string used to separate values when collapsing.
#' @param force_paste logical. TRUE unique value, FALSE all values.
#'
#' @return A character string: empty, single value, or collapsed unique values.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
collapse_conditional <- function(x, sep, force_paste = FALSE) {
  x <- x[x != "" & !is.na(x)]
  if (!length(x)) return("")
  if (force_paste) {
    paste(x, collapse = sep)
  } else {
    collapse_or_empty(x, sep)
  }
}


# groupuniqueids
#' Assign group identifiers to lipid features based on chemical identity
#' and proximity in retention time.
#'
#' This function groups lipid observations into unique GroupIDs by comparing
#' their chemical identity (LipidClass, LipidSubclass, C, DB) and retention
#' time. Two observations are assigned to the same group if they share the same
#' chemical identity and their retention times differ by no more than
#' \code{rttol_id}.
#'
#' The function operates in two modes:
#' \itemize{
#'   \item \strong{filtering}: Assigns GroupIDs sequentially within
#'   \code{rtdata$filtering$data}.
#'   \item \strong{prediction}: Assigns GroupIDs to predicted candidates by
#'   matching them to existing filtering groups when possible.
#' }
#'
#' Rows with missing values in any of the required columns are always assigned
#' new, unique GroupIDs.
#'
#' @param rtdata A list containing filtering and/or prediction results.
#' @param rt_col Character. Name of the retention time column.
#' @param rttol_id Numeric. Maximum allowed retention-time difference (in
#'   seconds) to consider two observations part of the same group.
#' @param typeanalysis Character. Either "filtering" or "prediction".
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
# groupuniqueids <- function(rtdata, rt_col = "RTsec", rttol_id = 5, typeanalysis) {
#
#   typeanalysis <- match.arg(typeanalysis, choices = c("filtering", "prediction"))
#
#   required_cols <- c("LipidClass", "LipidSubclass", "C", "DB", rt_col)
#
#   ## ======================================================
#   ## FILTERING MODE
#   ## ======================================================
#   if (typeanalysis == "filtering") {
#
#     data <- rtdata$filtering$data
#
#     if (!all(c(required_cols, "adduct") %in% names(data))) {
#       stop("Missing required columns in filtering$data")
#     }
#
#     o <- order(
#       data$LipidClass,
#       data$LipidSubclass,
#       data$C,
#       data$DB,
#       data[[rt_col]],
#       na.last = TRUE
#     )
#
#     data <- data[o, ]
#     rownames(data) <- NULL
#
#     n <- nrow(data)
#
#     if (n == 0L) {
#       data$GroupID <- integer(0)
#       data$adducts_expanded <- character(0)
#       rtdata$filtering$data <- data
#       return(rtdata)
#     }
#
#     complete_row <- !is.na(data$LipidClass) &
#       !is.na(data$LipidSubclass) &
#       !is.na(data$C) &
#       !is.na(data$DB) &
#       !is.na(data[[rt_col]])
#
#     same_lipid <- c(
#       FALSE,
#       complete_row[-1] & complete_row[-n] &
#         data$LipidClass[-1]    == data$LipidClass[-n] &
#         data$LipidSubclass[-1] == data$LipidSubclass[-n] &
#         data$C[-1]             == data$C[-n] &
#         data$DB[-1]            == data$DB[-n]
#     )
#
#     rt_close <- c(
#       FALSE,
#       complete_row[-1] & complete_row[-n] &
#         abs(diff(data[[rt_col]])) <= rttol_id
#     )
#
#     new_group <- !complete_row | !(same_lipid & rt_close)
#
#     data$GroupID <- cumsum(new_group)
#
#     ## ======================================================
#     ## EXPAND ADDUCTS PER GROUPID
#     ## ======================================================
#     data$adducts_expanded <- vapply(
#       seq_len(nrow(data)),
#       function(i) {
#
#         gid <- data$GroupID[i]
#
#         all_adducts <- unique(
#           data$adduct[data$GroupID == gid & !is.na(data$adduct)]
#         )
#
#         current_adduct <- data$adduct[i]
#
#         ordered_adducts <- c(
#           current_adduct,
#           all_adducts[all_adducts != current_adduct]
#         )
#
#         paste(ordered_adducts, collapse = "|")
#       },
#       character(1)
#     )
#
#     rtdata$filtering$data <- data
#     return(rtdata)
#   }
#
#   ## ======================================================
#   ## PREDICTION MODE
#   ## ======================================================
#   if (typeanalysis == "prediction") {
#
#     filtering  <- rtdata$filtering$data
#     candidates <- rtdata$prediction$candidates
#
#     if (!all(c(required_cols, "GroupID", "adduct") %in% names(filtering))) {
#       stop("filtering$data must contain GroupID and adduct")
#     }
#
#     if (!all(c(required_cols, "adduct") %in% names(candidates))) {
#       stop("Missing required columns in prediction$candidates")
#     }
#
#     if (any(is.na(filtering$GroupID))) {
#       stop("filtering$data contains NA GroupID values")
#     }
#
#     ## ======================================================
#     ## HANDLE EMPTY candidates
#     ## ======================================================
#     if (nrow(candidates) == 0) {
#       candidates$GroupID <- integer(0)
#       candidates$adducts_expanded <- character(0)
#       rtdata$prediction$candidates <- candidates
#       return(rtdata)
#     }
#
#     ## ======================================================
#     ## INIT
#     ## ======================================================
#     current_max <- max(filtering$GroupID)
#     candidates$GroupID <- rep(NA_integer_, nrow(candidates))
#
#     chem_key <- paste(
#       filtering$LipidClass,
#       filtering$LipidSubclass,
#       filtering$C,
#       filtering$DB,
#       sep = "_"
#     )
#
#     chem_index <- split(filtering, chem_key)
#
#     ## ======================================================
#     ## MAIN LOOP
#     ## ======================================================
#     for (i in seq_len(nrow(candidates))) {
#
#       cand <- candidates[i, ]
#
#       if (any(is.na(cand[required_cols]))) {
#         current_max <- current_max + 1L
#         candidates$GroupID[i] <- current_max
#         next
#       }
#
#       key <- paste(
#         cand$LipidClass,
#         cand$LipidSubclass,
#         cand$C,
#         cand$DB,
#         sep = "_"
#       )
#
#       assigned <- FALSE
#
#       if (key %in% names(chem_index)) {
#
#         ref <- chem_index[[key]]
#         group_ids <- unique(ref$GroupID)
#
#         for (gid in group_ids) {
#
#           rts <- ref[[rt_col]][ref$GroupID == gid]
#
#           if (any(abs(rts - cand[[rt_col]]) <= rttol_id)) {
#             candidates$GroupID[i] <- gid
#             assigned <- TRUE
#             break
#           }
#         }
#       }
#
#       if (!assigned) {
#         current_max <- current_max + 1L
#         candidates$GroupID[i] <- current_max
#       }
#     }
#
#     ## ======================================================
#     ## EXPAND ADDUCTS USING SAME GROUPID
#     ## ======================================================
#     combined <- rbind(
#       filtering[, c("GroupID", "adduct")],
#       candidates[, c("GroupID", "adduct")]
#     )
#
#     candidates$adducts_expanded <- vapply(
#       seq_len(nrow(candidates)),
#       function(i) {
#
#         gid <- candidates$GroupID[i]
#
#         all_adducts <- unique(
#           combined$adduct[combined$GroupID == gid & !is.na(combined$adduct)]
#         )
#
#         current_adduct <- candidates$adduct[i]
#
#         ordered_adducts <- c(
#           current_adduct,
#           all_adducts[all_adducts != current_adduct]
#         )
#
#         paste(ordered_adducts, collapse = "|")
#       },
#       character(1)
#     )
#
#     rtdata$prediction$candidates <- candidates
#     return(rtdata)
#   }
# }
groupuniqueids <- function(rtdata, rt_col = "RTsec", rttol_id = 5, typeanalysis) {

  typeanalysis <- match.arg(typeanalysis, choices = c("filtering", "prediction"))
  required_cols <- c("LipidClass", "LipidSubclass", "C", "DB", rt_col)

  ## ======================================================
  ## FILTERING
  ## ======================================================
  if (typeanalysis == "filtering") {

    data <- rtdata$filtering$data

    if (!all(c(required_cols, "adduct") %in% names(data))) {
      stop("Missing required columns in filtering$data")
    }

    o <- order(
      data$LipidClass,
      data$LipidSubclass,
      data$C,
      data$DB,
      data[[rt_col]],
      na.last = TRUE
    )

    data <- data[o, ]
    rownames(data) <- NULL

    n <- nrow(data)

    if (n == 0L) {
      data$GroupID <- integer(0)
      data$adducts_expanded <- character(0)
      rtdata$filtering$data <- data
      return(rtdata)
    }

    complete_row <- !is.na(data$LipidClass) &
      !is.na(data$LipidSubclass) &
      !is.na(data$C) &
      !is.na(data$DB) &
      !is.na(data[[rt_col]])

    same_lipid <- c(
      FALSE,
      complete_row[-1] & complete_row[-n] &
        data$LipidClass[-1]    == data$LipidClass[-n] &
        data$LipidSubclass[-1] == data$LipidSubclass[-n] &
        data$C[-1]             == data$C[-n] &
        data$DB[-1]            == data$DB[-n]
    )

    rt_close <- c(
      FALSE,
      complete_row[-1] & complete_row[-n] &
        abs(diff(data[[rt_col]])) <= rttol_id
    )

    data$GroupID <- cumsum(!complete_row | !(same_lipid & rt_close))

    ## expand adducts
    data$adducts_expanded <- vapply(
      seq_len(nrow(data)),
      function(i) {
        gid <- data$GroupID[i]
        ads <- unique(data$adduct[data$GroupID == gid & !is.na(data$adduct)])
        paste(c(data$adduct[i], ads[ads != data$adduct[i]]), collapse = "|")
      },
      character(1)
    )

    rtdata$filtering$data <- data
    return(rtdata)
  }

  ## ======================================================
  ## PREDICTION
  ## ======================================================
  if (typeanalysis == "prediction") {

    filtering  <- rtdata$filtering$data
    candidates <- rtdata$prediction$candidates

    if (nrow(candidates) == 0) {
      candidates$GroupID <- integer(0)
      candidates$adducts_expanded <- character(0)
      rtdata$prediction$candidates <- candidates
      return(rtdata)
    }

    current_max <- max(filtering$GroupID)
    candidates$GroupID <- NA_integer_
    candidates$from_filtering <- FALSE  # 🔥 clave

    ## índice químico
    chem_key <- paste(
      filtering$LipidClass,
      filtering$LipidSubclass,
      filtering$C,
      filtering$DB,
      sep = "_"
    )
    chem_index <- split(filtering, chem_key)

    ## ======================================================
    ## 1️⃣ MATCH REAL CON FILTERING
    ## ======================================================
    for (i in seq_len(nrow(candidates))) {

      cand <- candidates[i, ]

      if (any(is.na(cand[required_cols]))) {
        next
      }

      key <- paste(
        cand$LipidClass,
        cand$LipidSubclass,
        cand$C,
        cand$DB,
        sep = "_"
      )

      if (key %in% names(chem_index)) {

        ref <- chem_index[[key]]
        ref <- ref[order(ref[[rt_col]]), ]

        gids <- unique(ref$GroupID)

        for (gid in gids) {

          rts <- ref[[rt_col]][ref$GroupID == gid]

          if (
            cand[[rt_col]] >= min(rts) - rttol_id &&
            cand[[rt_col]] <= max(rts) + rttol_id
          ) {
            candidates$GroupID[i] <- gid
            candidates$from_filtering[i] <- TRUE
            break
          }
        }
      }
    }

    ## ======================================================
    ## 2️⃣ AGRUPAR SOLO LOS NUEVOS ENTRE SÍ
    ## ======================================================
    new_idx <- which(!candidates$from_filtering)

    if (length(new_idx) > 0) {

      new_data <- candidates[new_idx, ]

      o <- order(
        new_data$LipidClass,
        new_data$LipidSubclass,
        new_data$C,
        new_data$DB,
        new_data[[rt_col]]
      )

      new_data <- new_data[o, ]
      n <- nrow(new_data)

      same_lipid <- c(FALSE,
                      new_data$LipidClass[-1] == new_data$LipidClass[-n] &
                        new_data$LipidSubclass[-1] == new_data$LipidSubclass[-n] &
                        new_data$C[-1] == new_data$C[-n] &
                        new_data$DB[-1] == new_data$DB[-n]
      )

      rt_close <- c(FALSE,
                    abs(diff(new_data[[rt_col]])) <= rttol_id
      )

      new_group <- c(TRUE,
                     !(same_lipid[-1] & rt_close[-1])
      )

      local_gid <- cumsum(new_group)

      new_ids <- current_max + local_gid

      ## devolver posiciones originales
      candidates$GroupID[new_idx[o]] <- new_ids

      current_max <- max(new_ids)
    }

    ## ======================================================
    ## 3️⃣ EXPAND ADDUCTS
    ## ======================================================
    combined <- rbind(
      filtering[, c("GroupID", "adduct")],
      candidates[, c("GroupID", "adduct")]
    )

    candidates$adducts_expanded <- vapply(
      seq_len(nrow(candidates)),
      function(i) {
        gid <- candidates$GroupID[i]
        ads <- unique(combined$adduct[combined$GroupID == gid & !is.na(combined$adduct)])
        paste(c(candidates$adduct[i], ads[ads != candidates$adduct[i]]), collapse = "|")
      },
      character(1)
    )

    rtdata$prediction$candidates <- candidates
    return(rtdata)
  }
}



# lipidnameexpanded
#' Expand lipid names based on chemically and temporally close references.
#'
#' This function assigns an expanded lipid name to each observation by
#' searching for chemically identical lipids (same LipidClass, LipidSubclass,
#' C and DB) within a specified retention time tolerance. Among all matching
#' candidates, the lipid name with the best (lowest) confidence level is selected.
#'
#' The function operates in two modes:
#' \itemize{
#'   \item \strong{filtering}: Expanded names are inferred using only the
#'   filtering dataset as reference.
#'   \item \strong{prediction}: Expanded names are inferred using both filtering
#'   data and predicted candidates as reference.
#' }
#'
#' If no suitable reference is found within the tolerance, the expanded lipid
#' name is left empty.
#'
#' @param rtdata A list containing filtering and/or prediction results.
#' @param typeanalysis Character. Either "filtering" or"prediction".
#' @param rt_col Character. Name of the retention time column.
#' @param rttol_name Numeric. Maximum allowed retention-time difference (in
#'   seconds) to consider two lipids as temporally close.
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
lipidnameexpanded <- function(rtdata,
                              typeanalysis = c("filtering", "prediction"),
                              rt_col = "RTsec",
                              rttol_name = 1) {

  typeanalysis <- match.arg(typeanalysis)

  required <- c("LipidClass", "LipidSubclass", "C", "DB",
                "RTsec", "lipidname", "confidenceLevel")

  ## ======================================================
  ## SELECT TARGET AND REFERENCE DATASETS
  ## ======================================================
  if (typeanalysis == "filtering") {
    target <- rtdata$filtering$data
    reference <- rtdata$filtering$data[, c(required, "raw_lipidname"), drop = FALSE]
  }

  if (typeanalysis == "prediction") {
    target <- rtdata$prediction$candidates

    ## HANDLE EMPTY TARGET EARLY
    if (nrow(target) == 0) {
      target$confidenceLevel <- factor(character(0),
                                       levels = c("FA position", "FA", "Subclass", "MS-only", "MSMS", "Predicted",
                                                  "Other", "Unsupported class", "Not pass QC", ""),
                                       ordered = TRUE
      )
      target$lipidname_expanded <- character(0)

      rtdata$prediction$candidates <- target
      return(rtdata)
    }

    ## SAFE ASSIGNMENT
    if (!"confidenceLevel" %in% colnames(target)) {
      target$confidenceLevel <- rep("Predicted", nrow(target))
    }

    target$confidenceLevel <- factor(
      target$confidenceLevel,
      levels = c("FA position", "FA", "Subclass", "MS-only", "MSMS", "Predicted",
                 "Other", "Unsupported class", "Not pass QC", ""),
      ordered = TRUE
    )

    reference <- rbind(
      rtdata$filtering$data[, required, drop = FALSE],
      target[, required, drop = FALSE]
    )
  }

  if (!all(required %in% names(target))) {
    stop("Missing required columns in target dataset")
  }
  if (!all(required %in% names(reference))) {
    stop("Missing required columns in reference dataset")
  }

  ## SAFE INIT
  target$lipidname_expanded <- rep("", nrow(target))

  ## ======================================================
  ## MAIN LOOP
  ## ======================================================
  for (i in seq_len(nrow(target))) {

    row_i <- target[i, ]

    if (row_i$confidenceLevel %in% c("Unsupported class", "Not pass QC")) {
      target$lipidname_expanded[i] <- row_i$raw_lipidname
      next
    }

    same_chem <-
      reference$LipidClass    == row_i$LipidClass &
      reference$LipidSubclass == row_i$LipidSubclass &
      reference$C             == row_i$C &
      reference$DB            == row_i$DB

    close_rt <-
      abs(reference[[rt_col]] - row_i[[rt_col]]) <= rttol_name

    idx <- which(same_chem & close_rt)
    if (length(idx) == 0) next

    ## Best confidence level (ordered factor)
    best_level <- min(reference$confidenceLevel[idx])

    ## All entries with that level
    best_idx <- idx[reference$confidenceLevel[idx] == best_level]

    ## Expanded names
    best_names <- unique(reference$lipidname[best_idx])

    target$lipidname_expanded[i] <- paste(best_names, collapse = "|")
  }

  ## ======================================================
  ## STORE RESULTS
  ## ======================================================
  if (typeanalysis == "filtering") {
    rtdata$filtering$data <- target
  }
  if (typeanalysis == "prediction") {
    rtdata$prediction$candidates <- target
  }

  return(rtdata)
}


# recoverlipidnamespredicted
#' This function recovers lipid predictions that were previously filtered out
#' when a higher-confidence annotation exists for the same lipid species.
#' It ensures that improved annotations are not lost during filtering steps,
#' by restoring lipid identifications when a matching prediction with higher
#' confidence is available.
#'
#' The function operates on the rtdata object and updates the annotation set
#' by reassigning lipid names where applicable, prioritising annotations with
#' higher confidence scores.
#'
#' @param rtdata rtdata object.
#'
#' @return rtdata object with recovered lipid predictions where higher-confidence
#' annotations are available.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' rtdata <- recoverlipidpredictions(rtdata)
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
recoverlipidnamespredicted <- function(rtdata) {

  # =========================
  # CHECK REQUIRED OBJECTS
  # =========================
  if (is.null(rtdata$prediction$candidates) ||
      is.null(rtdata$filtering$data)) {
    return(rtdata)
  }

  candidates <- rtdata$prediction$candidates
  filtering  <- rtdata$filtering$data

  # =========================
  # DETERMINE STATUS COLUMNS
  # =========================
  pred_status_col <- if ("Filter_final" %in% colnames(candidates)) {
    "Filter_final"
  } else {
    "checkcandidate"
  }

  filt_status_col <- if ("Filter_final" %in% colnames(filtering)) {
    "Filter_final"
  } else {
    "checkfilter"
  }

  # =========================
  # CONFIDENCE RANKING
  # =========================
  confidence_levels <- c(
    "FA position",
    "FA",
    "Subclass",
    "MS-only",
    "MSMS",
    "Predicted",
    "Other",
    "Unsupported class",
    "Not pass QC",
    ""
  )

  # =========================
  # KEPT PREDICTIONS
  # =========================
  kept_predictions <- candidates[[pred_status_col]] %in% c(
    "Kept",
    "Kept manually by the user"
  )

  # =========================
  # REMOVED FILTERING
  # =========================
  removed_filtering <- !filtering[[filt_status_col]] %in% c(
    "Kept",
    "Recovered with RT modelling",
    "Kept manually by the user"
  )

  # =========================
  # LOOP OVER PREDICTIONS
  # =========================
  idx_predictions <- which(kept_predictions)

  for (i in idx_predictions) {

    pred_row <- candidates[i, ]

    # =========================
    # MATCHING REMOVED LIPIDS
    # =========================
    matches <- which(
      removed_filtering &
        filtering$LipidRTid     == pred_row$LipidRTid &
        filtering$LipidSubclass == pred_row$LipidSubclass &
        filtering$C             == pred_row$C &
        filtering$DB            == pred_row$DB &
        filtering$adduct        == pred_row$adduct
    )

    # =========================
    # IF MATCHES EXIST
    # =========================
    if (length(matches) > 0) {

      matched_rows <- filtering[matches, , drop = FALSE]

      # =========================
      # ORDER BY CONFIDENCE
      # =========================
      matched_rows$confidence_rank <- match(
        matched_rows$confidenceLevel,
        confidence_levels
      )

      # NA -> worst confidence
      matched_rows$confidence_rank[
        is.na(matched_rows$confidence_rank)
      ] <- length(confidence_levels) + 1

      # keep highest confidence
      matched_rows <- matched_rows[
        order(matched_rows$confidence_rank),
        ,
        drop = FALSE
      ]

      best_match <- matched_rows[1, ]

      # =========================
      # REPLACE VALUES
      # =========================
      # candidates$lipidname_expanded[i] <- best_match$lipidname
      candidates$lipidname[i] <- best_match$lipidname

      if ("confidenceLevel" %in% colnames(candidates)) {
        candidates$confidenceLevel[i] <- best_match$confidenceLevel
      }
    }
  }

  # =========================
  # SAVE BACK
  # =========================
  rtdata$prediction$candidates <- candidates

  return(rtdata)
}


# Combine filtered and predicted lipid annotations
#' This function combines filtered lipid annotations and prediction results into
#' a single data frame. It is used to integrate the outputs generated during the
#' filtering and prediction steps, ensuring a unified representation of lipid
#' annotations for downstream analysis.
#'
#' Depending on the analysis workflow, the function can handle both filtering
#' and prediction outputs, merging them into a consistent format.
#'
#' @param rtdata rtdata object.
#' @param typeanalysis character. Type of analysis to combine. It can be
#' "filtering" or "prediction".
#'
#' @return rtdata object with combined lipid annotation and prediction results
#' in a single data frame.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' rtdata <- combinefilterpredictions(rtdata, typeanalysis = "filtering")
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
combinefilterpredictions <- function(rtdata, typeanalysis) {

  if (!("Filter_final" %in% colnames(rtdata$filtering$data))) {
    # rtdata$filtering$data$Filter_final <- rtdata$filtering$data$Filter_id
    rtdata$filtering$data$Filter_final <- rtdata$filtering$data$checkfilter
  }
  if (!("Filter_final" %in% colnames(rtdata$prediction$candidates))) {
    rtdata$prediction$candidates$Filter_final <- rtdata$prediction$candidates$checkcandidate
  }


  if (typeanalysis == "filtering") {

    data <- rtdata$filtering$data
    data_kept <- data[!is.na(data$lipidname) &
                        data$lipidname != "" &
                        data$Filter == "Kept", , drop = FALSE]

  } else if (typeanalysis == "prediction") {

    data_filter <- rtdata$filtering$data

    data_predictions <- rtdata$prediction$candidates
    data_predictions$Filter <- data_predictions$Filter_final
    data_predictions$checkfilter <- data_predictions$Filter_final

    cols <- intersect(colnames(data_predictions), colnames(data_filter))

    data_predictions <- data_predictions[, cols]
    data_filter <- data_filter[, cols]

    data <- rbind(data_predictions, data_filter)

    data_kept <- data[!is.na(data$lipidname) &
                        !is.na(data$checkfilter) &
                        data$lipidname != "" &
                        data$Filter != "Ignore duplicated" &
                        data$Filter_final %in% c("Kept", "Kept manually by the user",

                                                 "Recovered with RT modelling"), , drop = FALSE]
  }

  return(data_kept)

}
