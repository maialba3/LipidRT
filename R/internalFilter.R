# filterannotationsdata
#' Filtering of annotated data based on RT criteria that allows discerning
#' erroneous annotations.
#'
#' Filtering of annotated data that allows discerning erroneous annotations.
#' The filtering is done in two phases. The first filter is based on the order
#' in which the different classes of lipids should elute during chromatography,
#' allowing the differentiation between correct annotations and those associated
#' with in-source fragments, eliminating the latter. After this, a second filter will
#' be applied, using the number of carbons and double bonds of the lipids belonging
#' to each class to check if the annotated species follow a pattern consistent
#' with their composition, eliminating those that do not.
#'
#' @param rtdata rtdata object.
#' @param distC_filter integer. Minimum difference in number of carbons from which lipids
#' are considered distant from the rest of lipids from the same class.
#' @param minSpecies_filter numeric. minimum number of lipids in a cluster from which the cluster is considered correct from the rest.
#' @param distribution character. Data distribution for RT modelling. It can be "gamma" or "normal".
#' @param rtscorethr numeric. Minimum Score_rtpred value from which a prediction will be considered valid.
#' @param minC_predict numeric. Minimum number of carbon units from any valid model data or prediction beyond which predictions will be discarded. This value must be at least 3, with a default of 5.
#' @param rt_thr numeric. Threshold for the difference (in seconds) between the
#' observed retention time (RT) and the predicted retention time (RTsec_rtpred).
#' @param rttol_id Numeric. Maximum allowed retention-time difference (in
#' seconds) to consider two observations part of the same group.
#' @param filter_plot logical. If TRUE, filtering is plotted. If FALSE, filtering is not plotted
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' rtdata <- filterannotationsdata(rtdata)
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
filterannotationsdata <- function(rtdata, distC_filter, minSpecies_filter, distribution,
                                  rtscorethr, minC_predict, rt_thr,
                                  rttol_id, rttol_name, filter_plot) {

  # save rtscorethr parameter in rtdata
  rtdata$metadata$filtering$rtscorethr <- as.numeric(rtscorethr)

  # rtdata$filtering$data <- rtdata$data[!(is.na(rtdata$data$lipidname)) &
  #                                        rtdata$data$lipidname != "" &
  #                                        !(rtdata$data$confidenceLevel %in% c("Not pass QC", "Unsupported class")), , drop = FALSE]
  rtdata$filtering$data <- rtdata$data[
    (
      !(is.na(rtdata$data$lipidname)) &
        rtdata$data$lipidname != ""
    ) |
      rtdata$data$confidenceLevel %in% c("Unsupported class", "Not pass QC"),
    ,
    drop = FALSE
  ]

  ##############################################################################
  # Check for other common incorrect annotations
  rtdata$filtering$data <- checkisobaricspecies(data = rtdata$filtering$data)

  ##############################################################################
  # Check for in-source fragments to remove incorrect annotations
  rtdata <- checkinsourcefragments(rtdata = rtdata,
                                   typeanalysis = "filtering")

  ##############################################################################
  # Remove distant data points for RT modelling
  rtdata <- removedistantrtdatapoints(rtdata = rtdata,
                                      distC_filter = distC_filter,
                                      minSpecies_filter = minSpecies_filter,
                                      typeanalysis = "filtering")

  ##############################################################################
  # Check for RT coherence between different classes
  rtdata<- checkRTclasses(rtdata = rtdata,
                          typeanalysis = "filtering")

  ##############################################################################
  # Filtering of each lipid class according to its RT pattern based on its
  # composition of C and DB
  rtdata <- filterbyrtrules(rtdata = rtdata,
                            typeanalysis = "filtering")

  ##############################################################################
  # Check outliers in terms of number of carbons
  rtdata <- checkoutliers(rtdata = rtdata,
                          typeanalysis = "filtering")

  ##############################################################################
  # Modelling the retention time of each lipid class
  rtdata <- checkremovedids(rtdata = rtdata,
                            distribution = distribution,
                            rtscorethr = rtscorethr)
  # Group unique IDs
  rtdata <- groupuniqueids(rtdata,
                           typeanalysis = "filtering",
                           rttol_id = rttol_id)

  # Expand candidates lipid name
  rtdata <- lipidnameexpanded(rtdata,
                              typeanalysis = "filtering",
                              rttol_name = rttol_name)

  #============================================================================#
  # Plot filtering results
  #============================================================================#
  if (filter_plot) {
    rtdata <- filteringgraph(rtdata = rtdata)
  }

  if (is.null(rtdata$filtering$data$comments)) {
    rtdata$filtering$data$comments <- ""
  }

  return(rtdata)
}


# proposemzcandidates
#' Propose candidates for the features without annotations
#'
#' Propose candidates for the features without annotations from the classes of
#' interest, based on the mz of the features, considering the mz tolerance value
#' provided by the user. This value can be given in ppm (parts per million) or
#' absolute mass error
#'
#' @param rtdata rtdata
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
proposemzcandidates <- function(rtdata,
                                filedbs = NULL,
                                sep = ",",
                                classes = NULL,
                                adducts = NULL,
                                dmz = 5,
                                ppm = TRUE) {

  polarity <- rtdata$metadata$polarity

  if (is.null(polarity) || polarity == "both") {
    return(rtdata)
  }

  # -------------------------------------------------------------------------
  # Identify LipidRTid containing ether phospholipid annotations
  # PX-O, PX-P, LPX-O, LPX-P
  # -------------------------------------------------------------------------
  if(!"comments" %in% colnames(rtdata$data)) {
    rtdata$data$comments <- ""
  }
  data <- rtdata$data


  # lipid_pattern <- "^(P[A-Z]|P[A-Z]-(O|P)|LP[A-Z]-(O|P))$"
  lipid_pattern <- "^(P[A-Z]-(O|P)|LP[A-Z]-(O|P))$"

  lipid_rows <- data[grepl(lipid_pattern, data$LipidSubclass, perl = TRUE), , drop = FALSE]

  if (nrow(lipid_rows) == 0) {
    rtdata$filtering$isomeric_plasmalogens <- data.frame()
    return(rtdata)
  }



  # -------------------------------------------------------------------------
  # Keep LipidRTid with unique mz and RT
  # -------------------------------------------------------------------------
  lipidrtids <- unique(lipid_rows$LipidRTid)

  target_rows <- list()

  for (id in lipidrtids) {

    tmp <- data[data$LipidRTid == id, , drop = FALSE]

    unique_mz <- unique(round(tmp$mz, 3))
    unique_rt <- unique(round(tmp$RTsec, 3))

    unique_subclasses <- unique(tmp$LipidSubclass)
    unique_confidence <- unique(tmp$confidenceLevel)

    # Require: unique mz, unique RT, more than one subclass,                                  # more than one subclass,           # all annotations must have confidenceLevel == "Subclass"
    if (length(unique_mz) == 1 && length(unique_rt) == 1) { # && length(unique_subclasses) != 1 # && length(unique_confidence) == 1 && unique_confidence == "Subclass"
      target_rows[[as.character(id)]] <- tmp[1, , drop = FALSE]
    }
  }

  if (length(target_rows) == 0) {
    rtdata$filtering$isomeric_plasmalogens <- data.frame()
    return(rtdata)
  }

  targets <- do.call(rbind, target_rows)

  # -------------------------------------------------------------------------
  # Generate mz databases
  # -------------------------------------------------------------------------
  generateadductsmzdbs <- createadductsmzdbs(filedbs = filedbs,
                                             sep = sep,
                                             classes = classes,
                                             adducts = adducts)

  filedbs <- generateadductsmzdbs$filedbs
  adductsmzdbs <- generateadductsmzdbs$adductsmzdbs

  if (is.character(adductsmzdbs)) {
    return(adductsmzdbs)
  }

  # -------------------------------------------------------------------------
  # Valid adducts according to polarity
  # -------------------------------------------------------------------------
  if (polarity == "positive") {
    valid_adducts <- adductsmzdbs$adductsTable$adduct[adductsmzdbs$adductsTable$charge > 0]
  } else if (polarity == "negative") {
    valid_adducts <- adductsmzdbs$adductsTable$adduct[adductsmzdbs$adductsTable$charge < 0]
  }

  invalid_adducts <- adductsmzdbs$adductsTable$adduct[!adductsmzdbs$adductsTable$adduct %in% valid_adducts]

  # -------------------------------------------------------------------------
  # Check adducts by LipidSubclass
  # -------------------------------------------------------------------------
  # subclasses <- unique(targets$LipidSubclass)
  target_subclasses <- unique(targets$LipidSubclass)

  adduct_check <- list()
  messages <- character()

  for (subclass in target_subclasses) {

    cand.subclass <- unique(
      ifelse(grepl("-O$", subclass), sub("-O$", "-P", subclass),
        ifelse(grepl("-P$", subclass), sub("-P$", "-O", subclass),
          NA_character_
        )
      )
    )
    cand.subclass <- cand.subclass[!is.na(cand.subclass)]

    target_adducts <- unique(targets$adduct[targets$LipidSubclass == subclass])

    # Aductos inválidos por polaridad
    invalid_polarity <- target_adducts[!target_adducts %in% valid_adducts]

    # Aductos no presentes en la base de datos de esa clase
    missing_db <- character()

    if (!is.null(adductsmzdbs[[subclass]])) {
      missing_db <- target_adducts[!target_adducts %in% colnames(adductsmzdbs[[cand.subclass]])]
    } else {
      # Si no existe la clase en la base, todos faltan
      missing_db <- target_adducts
    }

    # Guardar resultados
    adduct_check[[subclass]] <- list(
      invalid_polarity = invalid_polarity,
      missing_db = missing_db
    )

    # Construir mensaje
    msg <- character()

    if (length(invalid_polarity) > 0) {
      msg <- c(
        msg,
        paste0(
          subclass,
          ": invalid adduct(s) for ",
          polarity,
          " mode: ",
          paste(invalid_polarity, collapse = ", ")
        )
      )
    }

    if (length(missing_db) > 0) {
      msg <- c(
        msg,
        paste0(
          subclass,
          ": adduct(s) not available in database: ",
          paste(missing_db, collapse = ", "),
          ". If you want to use these adduct(s), add them to the corresponding database through the 'filedbs' parameter by completing the appropriate Excel file."
        )
      )
    }

    messages <- c(messages, msg)
  }

  # Mostrar mensaje
  if (length(messages) > 0) {
    message(paste(messages, collapse = "\n"))
  }


  # -------------------------------------------------------------------------
  # Propose ether candidates
  # -------------------------------------------------------------------------

  candidates <- data.frame()

  for (i in seq_len(nrow(targets))) {

    LipidRTid <- targets$LipidRTid[i]
    mz <- as.numeric(targets$mz[i])
    rt <- as.numeric(targets$RTsec[i])
    adduct_lipid <- targets$adduct[i]

    # existing_subclasses <- unique(data$LipidSubclass[data$LipidRTid == LipidRTid])
    #
    # # Only search ether subclasses
    # lipid_subclasses <- names(adductsmzdbs)[grepl(lipid_pattern, names(adductsmzdbs), perl = TRUE)]
    #
    # # Remove already assigned subclasses
    # lipid_subclasses <- setdiff(lipid_subclasses, existing_subclasses)

    existing_subclasses <- unique(
      data$LipidSubclass[data$LipidRTid == LipidRTid]
    )

    lipid_subclasses <- unique(unlist(lapply(existing_subclasses, function(x) {

      if (grepl("-O$", x)) {
        sub("-O$", "-P", x)
      } else if (grepl("-P$", x)) {
        sub("-P$", "-O", x)
      } else {
        NULL
      }

    })))

    lipid_subclasses <- intersect(lipid_subclasses, names(adductsmzdbs))

    new_rows <- data.frame()
    roworder_prev <- as.numeric(max(rtdata$data$roworder[rtdata$data$LipidRTid %in% LipidRTid]))

    for (LipidSubclass.cand in lipid_subclasses) {

      mz_db <- adductsmzdbs[[LipidSubclass.cand]]

      # adduct_columns <- names(mz_db[!(names(mz_db) %in% c("formula", "total",
      #                                                     "Mass", "ID",
      #                                                     invalid_adducts)) &
      #                                 names(mz_db) %in% adduct_lipid])

      adduct_columns <- adduct_lipid # que tengan mismo aducto

      if (length(adduct_columns) == 0 | !(adduct_columns) %in% colnames(mz_db)) {
        next
      }

      # ---------------------------------------------------------------------
      # mz matching
      # ---------------------------------------------------------------------

      if (ppm == TRUE) {
        matches <- apply(mz_db[adduct_columns], 2, function(x) {
          (((abs(as.numeric(x) - mz)) / x) * 1e6 <= as.numeric(dmz))
        })
      } else {
        matches <- apply( mz_db[adduct_columns], 2, function(x) {
          (abs(as.numeric(x) - mz) <= as.numeric(dmz))
        })
      }

      # ---------------------------------------------------------------------
      # Store candidates
      # ---------------------------------------------------------------------

      if (length(which(matches)) > 0) {

        for (cand in which(matches)) {

          CDB <- rep(mz_db$total, length(adduct_columns))[cand]

          lipidname.cand <- paste0(LipidSubclass.cand, "(", CDB, ")")

          lipidname.cand <- gsub("\\s+", "", lipidname.cand)

          LipidClass.cand <- sub("-(O|P)$", "", LipidSubclass.cand)

          C.cand <- as.numeric(unlist(strsplit(CDB, ":"))[1])
          DB.cand <- as.numeric(unlist(strsplit(CDB, ":"))[2])

          adducts.cand <- adduct_columns[ceiling(cand / dim(mz_db)[1])]

          # Check whether this lipid is already annotated
          # for the same LipidRTid
          already_exists <- any(data$LipidRTid == LipidRTid &
                                  data$LipidSubclass == LipidSubclass.cand &
                                  data$C == C.cand &
                                  data$DB == DB.cand &
                                  data$adduct == adducts.cand)

          # Skip if already annotated
          if (already_exists) {
            next
          }

          # ------------------------------------------------------------
          # Cross ether/plasmalogen chain borrowing (O <-> P)
          # ------------------------------------------------------------

          adjusted_lipidname <- lipidname.cand
          adjusted_DB <- DB.cand
          adjusted_confidenceLevel <- "Subclass"

          # helper dataset (features of same LipidRTid)
          feat_subset <- data[data$LipidRTid == LipidRTid, , drop = FALSE]

          # -----------------------------
          # CASE 1: candidate is O -> look for P
          # Example:
          # PE-O(16:0_22:6) <- from PE-P(16:0_22:5)
          # First chain DB changes (+1)
          # -----------------------------
          if (grepl("-O$", LipidSubclass.cand) & DB.cand > 0) {

            match_feat <- feat_subset[
              feat_subset$LipidSubclass == sub("-O$", "-P", LipidSubclass.cand) &
                feat_subset$C == C.cand &
                feat_subset$DB == (DB.cand - 1) &
                feat_subset$adduct == adducts.cand &
                feat_subset$confidenceLevel != "Subclass",
              ,
              drop = FALSE
            ]

            if (nrow(match_feat) > 0) {

              # Extract chains from matched feature
              p_chain <- gsub(".*\\((.*)\\)", "\\1", match_feat$lipidname[1])
              p_chain <- gsub("P-", "", p_chain)

              # Detect chain separator ("_" or "/")
              separator <- if (grepl("_", p_chain, fixed = TRUE)) {
                "_"
              } else if (grepl("/", p_chain, fixed = TRUE)) {
                "/"
              } else {
                NULL
              }

              if (!is.null(separator)) {

                parts <- strsplit(p_chain, separator, fixed = TRUE)[[1]]

                if (length(parts) == 2) {

                  # Modify FIRST chain DB
                  first <- strsplit(parts[1], ":")[[1]]

                  first_db <- as.numeric(first[2]) + 1

                  if (!(first_db < 0)) {
                    parts[1] <- paste0(first[1], ":", first_db)

                    adjusted_lipidname <- paste0(LipidClass.cand,"(O-",paste(parts, collapse = separator),")")
                  } else {
                    next
                  }

                  adjusted_DB <- DB.cand

                  adjusted_confidenceLevel <- match_feat$confidenceLevel
                }
              }
            }
          }

          # -----------------------------
          # CASE 2: candidate is P -> look for O
          # Example:
          # PE-P(16:0_22:5) <- from PE-O(16:1_22:5)
          # First chain DB changes (-1)
          # -----------------------------
          if (grepl("-P$", LipidSubclass.cand)) {

            match_feat <- feat_subset[
              feat_subset$LipidSubclass == sub("-P$", "-O", LipidSubclass.cand) &
                feat_subset$C == C.cand &
                feat_subset$DB == (DB.cand + 1) &
                feat_subset$adduct == adducts.cand &
                feat_subset$confidenceLevel != "Subclass",
              ,
              drop = FALSE
            ]

            if (nrow(match_feat) > 0) {

              o_chain <- gsub(".*\\((.*)\\)", "\\1", match_feat$lipidname[1])
              o_chain <- gsub("O-", "", o_chain)

              # Detect chain separator ("_" or "/")
              separator <- if (grepl("_", o_chain, fixed = TRUE)) {
                "_"
              } else if (grepl("/", o_chain, fixed = TRUE)) {
                "/"
              } else {
                NULL
              }

              if (!is.null(separator)) {

                parts <- strsplit(o_chain, separator, fixed = TRUE)[[1]]

                if (length(parts) == 2) {

                  # Modify FIRST chain DB
                  first <- strsplit(parts[1], ":")[[1]]

                  first_db <- as.numeric(first[2]) - 1

                  if (!(first_db < 0)) {
                    parts[1] <- paste0(first[1], ":", first_db)

                    adjusted_lipidname <- paste0(LipidClass.cand,"(P-",paste(parts, collapse = separator),")")
                  } else {
                    next
                  }

                  adjusted_DB <- DB.cand

                  adjusted_confidenceLevel <- match_feat$confidenceLevel
                }
              }
            }
          }

          new_rows <- rbind(new_rows,
                            data.frame(
                              mz = mz,
                              RTsec = rt,
                              LipidRTid = LipidRTid,
                              raw_lipidname = lipidname.cand, # de momento esto por si al agrupar se usa ese o algo -- comprobar
                              lipidname = adjusted_lipidname,
                              adduct = adducts.cand,
                              IS = "",
                              confidenceLevel = adjusted_confidenceLevel,
                              # confidenceLevel = "Subclass",
                              LipidClass = LipidClass.cand,
                              LipidSubclass = LipidSubclass.cand,
                              C = C.cand,
                              DB = DB.cand,
                              Filter = "Kept",
                              stringsAsFactors = FALSE
                            ))

        }
      }
    }

    if (nrow(new_rows) == 0) next
    new_rows$roworder <- seq(roworder_prev + 1, roworder_prev + nrow(new_rows))

    candidates <- rbind(candidates, new_rows)
  }

  # -------------------------------------------------------------------------
  # Finalize
  # -------------------------------------------------------------------------

  if (nrow(candidates) > 0) {
    candidates <- candidates[order(candidates$LipidRTid, candidates$LipidSubclass,
                                   candidates$C, candidates$DB), , drop = FALSE]

    # añadir indice
    candidates$index <- seq(max(rtdata$data$index) + 1,
                            max(rtdata$data$index) + nrow(candidates))

    # grouplipidclass
    # candidates <- grouplipidclasses(candidates,  map_type = "class")

    candidates$comments <- "Proposed isomeric species"

  } else {

    candidates <- data.frame(
      mz = numeric(0),
      RTsec = numeric(0),
      LipidRTid = integer(0),
      raw_lipidname = character(0),
      lipidname = character(0),
      adduct = character(0),
      IS = character(0),
      confidenceLevel = character(0),
      LipidClass = character(0),
      LipidSubclass = character(0),
      C = integer(0),
      DB = integer(0),
      Filter = character(0),
      index = numeric(0)

    )
  }

  rtdata$filtering$isomeric_plasmalogens <- candidates

  rtdata$data <- rbind(rtdata$data, candidates)

  return(rtdata)
}


# checkcommonerrors
#' Check for other common incorrect annotations in specific classes
#'
#' Check for other common incorrect annotations in specific classes
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkcommonerrors <- function(data){

  data_kept <- data[!is.na(data$lipidname) &
                      data$lipidname != "" &
                      data$Filter == "Kept", , drop = FALSE]

  # If a feature is identified as MG or Carnitines, we keep only Carnitines
  to_remove <- c()
  for (g in unique(data_kept$LipidRTid[data_kept$LipidClass == "MG"])) {
    groupids <- data_kept[data_kept$LipidRTid == g,]
    if(nrow(groupids) > 1) {
      if ("Carnitine" %in% groupids$LipidClass) {
        to_remove <- c(to_remove, groupids$index[groupids$LipidClass=="MG"])
      }
    }
  }

  # if a feature is identified as PXx and LPXx, keep only LPXx
  fragment_class <- c("PC-O", "PC-P", "PE-O", "PE-P")
  for (class in fragment_class) {
    compare_class <- paste0(paste0("L", substr(class, 1, 2)), c("", "-O", "-P"))
    for (g in unique(data_kept$LipidRTid[data_kept$LipidClass == class])) {
      groupids <- data_kept[data_kept$LipidRTid == g,]
      if(nrow(groupids) > 1) {
        if (any(compare_class %in% groupids$LipidClass)) {
          to_remove <- c(to_remove, groupids$index[groupids$LipidClass==class])
        }
      }
    }
  }

  # Check for other common incorrect annotations in plasmalogens
  # for PLp/PLp keep only those with the highest confidence level
  # fragment_class <- c("PC-O", "PC-P", "PE-O", "PE-P")
  # confLevels <- LipidRT::confLevels
  # for (g in unique(data_kept$LipidRTid[data_kept$LipidClass %in% fragment_class])){
  #   groupids <- data_kept[data_kept$LipidRTid == g,]
  #   if(nrow(groupids) > 1){
  #     conflev <- confLevels$order[match(groupids$confidenceLevel, confLevels$level)]
  #     to_remove <- c(to_remove, groupids$index[conflev < max(conflev)])
  #   }
  # }

  data$Filter <- ifelse(data$index %in% to_remove, "Confidence level among isobars", data$Filter)

  return(data)
}


# checkisobaricspecies
#' Check for other common incorrect annotations in specific classes
#'
#' Check for other common incorrect annotations in specific classes
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkisobaricspecies <- function(data) {

  data_kept <- data[!is.na(data$lipidname) &
                      data$lipidname != "" &
                      data$Filter %in% c("Kept", "Kept manually by the user"), , drop = FALSE]

  confLevels <- LipidRT::confLevels
  confLevels$isobaric_order <- confLevels$order
  confLevels$isobaric_order[confLevels$isobaric_order == 4] <- 3

  to_remove <- c()

  for (g in unique(data_kept$LipidRTid)) {
    groupids <- data_kept[data_kept$LipidRTid == g, ]

    if (nrow(groupids) > 1) {

      combs <- combn(seq_len(nrow(groupids)), 2)

      has_diff <- any(apply(combs, 2, function(idx) {
        vals <- groupids[idx, c("LipidSubclass", "adduct")]

        if (anyNA(vals)) return(FALSE)

        !identical(vals[1, ], vals[2, ])
      }))

      if (isTRUE(has_diff)) {
        conflev <- confLevels$isobaric_order[
          match(groupids$confidenceLevel, confLevels$level)
        ]

        if (any(conflev == 3, na.rm = TRUE) &&
            any(conflev < 3, na.rm = TRUE)) {

          to_remove <- c(
            to_remove,
            groupids$index[conflev < 3 & !is.na(conflev)]
          )
        }
      }
    }
  }

  data$Filter <- ifelse(
    data$index %in% to_remove,
    "Confidence level among isobars",
    data$Filter
  )

  return(data)
}



# checkinsourcefragments
#' Check for in-source fragments to remove incorrect annotations
#'
#' Check for in-source fragments to remove incorrect annotations
#'
#' @param rtdata rtdata
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkinsourcefragments <- function(rtdata, typeanalysis = "filtering") {

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)


  # Remove in-source fragments (Glycerolipids, Lysophospholipids).
  to_remove <- c()
  fragment_subclass <- intersect(unique(data_kept$LipidSubclass),
                                 c("FA", "MG", "DG",
                                   "LPE", "LPC", "LPG", "LPI", "LPS", "LPA",
                                   "Cer_d", "Cer_t", "HexCer_d", "HexCer_t"))

  for (subclass in fragment_subclass) {
    if (subclass == "MG") {
      compare_subclass <- c("DG", "TG")

      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = FALSE)
      to_remove <- append(to_remove, rem)
    }

    if (subclass == "DG") {
      compare_subclass <- "TG"
      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = FALSE)
      to_remove <- append(to_remove, rem)

      compare_subclass <- "PL"
      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = TRUE)
      to_remove <- append(to_remove, rem)

      compare_subclass <- c("MGDG", "DGDG")
      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = TRUE)
      to_remove <- append(to_remove, rem)

    }

    if (subclass %in% fragment_subclass[grep("^LP", fragment_subclass)]) {
      compare_subclass <- substr(subclass, 2, nchar(subclass))
      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = FALSE)
      to_remove <- append(to_remove, rem)
    }

    if (subclass == "FA") {
      compare_subclass <- "LPL"
      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = TRUE)
      to_remove <- append(to_remove, rem)

      compare_subclass <- "PL"
      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = FALSE)
      to_remove <- append(to_remove, rem)

    }

    if (subclass %in% c("HexCer_d", "HexCer_t")) {
      if (subclass %in% c("HexCer_d")) {
        compare_subclass <- c("SHexCer_d", "AHexCer_d")
      } else if (subclass %in% c("HexCer_t")){
        compare_subclass <- c("SHexCer_t", "AHexCer_t")
      }

      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = TRUE) # same composition
      to_remove <- append(to_remove, rem)
    }

    if (subclass %in% c("Cer_d", "Cer_t")) {
      if (subclass %in% c("Cer_d")) {
        compare_subclass <- c("SHexCer_d", "HexCer_d", "SM_d", "CerP_d")
      } else if (subclass %in% c("Cer_t")){
        compare_subclass <- c("SHexCer_t", "HexCer_t", "SM_t", "CerP_t")
      }

      rem <- getrtclusters(data_kept,
                           lipid_subclass = subclass,
                           compare_subclass = compare_subclass,
                           sameC = TRUE) # same composition
      to_remove <- append(to_remove, rem)
    }

  }

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "In-source fragments",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - In-source fragments",
                                                        # "Discarded Filter",
                                                  rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "In-source fragments",
                                                 rtdata$filtering$data$Filter_final)
  }

  return(rtdata)
}


# getrtclusters
#' Detect erroneous annotations due to in-source fragments by RT clustering
#'
#' Detect erroneous annotations due to in-source fragments by RT clustering
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#' @param lipid_subclass character. Specific lipid class on which I want to detect erroneous annotations
#' @param compare_subclass Specific lipid class that will be used to compare retention times
#' @param maxdist maximum RT distance between points within a cluster
#' @param sameC logical. In-source fragments are supposed to have the same number of carbons
#' as the original lipid class.
#' @param rtdiff difference RT.
#' @param dist_percentile distance percentile from all species to be used as maxdist.
#'
#' @details If in-source fragments are expected to have the same number of carbons
#' than the original class, species with the same C:DB composition are checked
#' with a two second windows. If any overlap is found, potential fragments are removed.
#' If same number of carbons is not expected, unsupervised clustering is performed
#' based on RT and number of C. Subsequently, the retention time of each cluster
#' will be compared with the retention time of the reference lipid class. Clusters
#' will be eliminated if at least 70% of the species fall within the retention time
#' interval of the 95% of the data of  reference class.
#'
#' @return to_remove
#' vector containing the indices of the data points that are removed
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
getrtclusters <- function(data, lipid_subclass, compare_subclass, maxdist = 90,
                          sameC, rtdiff = 1, dist_percentile = 0.25) {

  # rename classes
  if (any(compare_subclass == "LPL")) {
    compare_subclass <- c("LPC", "LPE", "LPI", "LPS", "LPG", "LPC-O", "LPC-P",
                          "LPE-O", "LPE-P", "LPA")
  } else if (any(compare_subclass == "PL")) {
    compare_subclass <- c("PC", "PE", "PI", "PS", "PG", "PA")
  }

  if (lipid_subclass %in% unique(data$LipidSubclass) &
      any(compare_subclass %in% unique(data$LipidSubclass))) {

    data_subclass <- data[data$LipidSubclass == lipid_subclass, c("index", "RTsec", "C", "DB")]

    # if comparing classes are expected to have the same C:DB composition
    if (sameC) {
      perc_inside <- c()
      for (c in data_subclass$index) {
        x <- data_subclass[data_subclass$index == c,]
        dataindex <- (data$LipidSubclass %in% compare_subclass &
                        data$C == x$C &
                        data$DB == x$DB)
        # check coelution based on rtdiff
        if (any(dataindex)) {
          perc_inside <- append(perc_inside,
                                as.numeric(any(abs(data$RTsec[dataindex] - x$RTsec) <= rtdiff)))
        } else {
          perc_inside <- append(perc_inside, 0)
        }
      }
      # if different C:DB composition
    } else {
      # clusterize the class that could contain in-source fragments
      clusters <- clust_insourcefrags(rt = data_subclass$RTsec,
                                      ncarbons = data_subclass$C,
                                      scale.data = TRUE,
                                      rt_weight = 2,
                                      nc_weight = 0.5,
                                      maxdist = NULL,
                                      linkage = "center",
                                      dist_percentile = dist_percentile)

      nclusters <- length(unique(clusters))
      data_subclass$Cluster <- as.numeric(as.factor(clusters))

      # check overlapping between clusters and reference class
      dataindex <- data$LipidSubclass %in% compare_subclass
      percentil_95 <- c(quantile(data$RTsec[dataindex], 0.025),
                        quantile(data$RTsec[dataindex], 0.975))
      perc_inside <- tapply(data_subclass$RTsec, data_subclass$Cluster,
                            function(x) mean(x > percentil_95[1] &
                                               x < percentil_95[2]))
    }
    # apply results
    if (sameC) {
      iremove <- data_subclass$index[which(perc_inside == 1)]
    } else {
      iremove <- data_subclass$index[data_subclass$Cluster %in% which(perc_inside >= 0.7)]
    }

    if (length(iremove) > 0) {
      to_remove <- data_subclass$index[data_subclass$index %in% iremove]
    } else {
      to_remove <- c()
    }
  } else {
    to_remove <- c()
  }

  return(to_remove)
}



# clust_insourcefrags
#' clustering for in-source fragments
#'
#' clustering for in-source fragments
#'
#' @param rt numeric vector for retention times
#' @param ncarbons integer vector for number of carbons
#' @param maxdist max distance intra cluster
#' @param scale.data logical. scale data to make RT sn ncarbons comparable
#' @param dist_percentile distance percentile from all species to be used as maxdist
#' @param rt_weight weight for scaled retention time
#' @param nc_weight weight for scaled number of carbons
#' @param linkage method for distance calculation
#'
#' @return vector
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
clust_insourcefrags <- function(rt, ncarbons, maxdist = NULL, scale.data = TRUE,
                                dist_percentile = 0.5, rt_weight = 2, nc_weight = 0.5,
                                linkage = c("center", "single", "complete")) {

  linkage <- match.arg(linkage)

  # scaling and weighting
  features <- cbind(rt, ncarbons)
  if (scale.data) {
    features[, 1] <- minmaxscale(features[, 1])
    features[, 2] <- minmaxscale(features[, 2])
  }
  features[, 1] <- features[, 1] * sqrt(rt_weight)
  features[, 2] <- features[, 2] * sqrt(nc_weight)

  original_features <- features  # Keep for pairwise distance calculations

  # estimate maxdist if missing
  dists <- stats::dist(features)
  if (is.null(maxdist)) {
    maxdist <- quantile(dists, probs = dist_percentile)
  }

  mins <- maxs <- features
  n_points <- nrow(features)

  if (n_points == 1) return(1)

  clust <- rep(0, n_points)
  n <- rep(1, n_points)
  atclust <- lapply(1:n_points, function(i) i)

  distmatrix <- compute_distmatrix(linkage, mins, maxs, original_features, atclust)

  while (nrow(features) > 1) {
    if(all(is.na(distmatrix))) {
      break
    }

    pos <- which(distmatrix == min(distmatrix, na.rm = TRUE), arr.ind = TRUE)[1, ]
    n1 <- pos[2]
    n2 <- pos[1]

    if (distmatrix[n2, n1] > maxdist) {
      distmatrix[n2, n1] <- NA
      next
    }

    idx_keep <- min(n1, n2)
    idx_remove <- max(n1, n2)

    mins[idx_keep, ] <- pmin(mins[n1, , drop = FALSE], mins[n2, , drop = FALSE])
    maxs[idx_keep, ] <- pmax(maxs[n1, , drop = FALSE], maxs[n2, , drop = FALSE])
    features[idx_keep, ] <- (features[n1, , drop = FALSE] * n[n1] +
                               features[n2, , drop = FALSE] * n[n2]) / (n[n1] + n[n2])
    n[idx_keep] <- n[n1] + n[n2]
    atclust[[idx_keep]] <- c(atclust[[n1]], atclust[[n2]])

    mins <- mins[-idx_remove, , drop = FALSE]
    maxs <- maxs[-idx_remove, , drop = FALSE]
    features <- features[-idx_remove, , drop = FALSE]
    n <- n[-idx_remove]
    atclust <- atclust[-idx_remove]

    distmatrix <- compute_distmatrix(linkage, mins, maxs, original_features, atclust)
  }

  for (c in seq_along(atclust)) {
    clust[atclust[[c]]] <- c
  }

  return(clust)
}

# minmaxscale
#' scale numeric vector between 0 and 1 using min and max values
#'
#' scale numeric vector between 0 and 1 using min and max values
#'
#' @param x numeric vector
#'
#' @return vector
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
minmaxscale <- function(x) {
  if (length(unique(x)) > 1) {
    (x - min(x)) / (max(x) - min(x))
  } else {
    rep(1, length(x))
  }
}

# compute_distmatrix
#' compute distance matriix
#'
#' compute distance matriix
#'
#' @param linkage method for distance calculation
#' @param mins numeric vector
#' @param maxs numeric vector
#' @param original_features matrix with two columns: rt and ncarbons
#' @param atclust cluster of each feature at a given point
#'
#' @return matrix
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
compute_distmatrix <- function(linkage, mins, maxs, original_features, atclust) {
  if (linkage == "center") {
    return(clustdist_center(mins, maxs))
  } else if (linkage == "single") {
    return(clustdist_single(original_features, atclust))
  } else if (linkage == "complete") {
    return(clustdist_complete(original_features, atclust))
  }
}

# clustdist_center
#' distance between centers of two clusters
#'
#' distance between centers of two clusters
#'
#' @param mins numeric vector
#' @param maxs numeric vector
#'
#' @return matrix
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
clustdist_center <- function(mins, maxs) {
  n <- nrow(mins)
  if (n < 2) {
    return(matrix(NA, n, n))
  }
  distmat <- matrix(NA, n, n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      center_i <- (mins[i, ] + maxs[i, ]) / 2
      center_j <- (mins[j, ] + maxs[j, ]) / 2
      distmat[i, j] <- sqrt(sum((center_i - center_j)^2))
    }
  }
  return(distmat)
}

# clustdist_single
#' min distance between two clusters
#'
#' min distance between two clusters
#'
#' @param features matrix with two columns: rt and ncarbons
#' @param atclust cluster of each feature at a given point
#'
#' @return matrix
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
clustdist_single <- function(features, atclust) {
  n <- length(atclust)
  if (n < 2) {
    return(matrix(NA, n, n))
  }
  distmat <- matrix(NA, n, n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      points_i <- features[atclust[[i]], , drop = FALSE]
      points_j <- features[atclust[[j]], , drop = FALSE]
      min_dist <- Inf
      for (p1 in 1:nrow(points_i)) {
        for (p2 in 1:nrow(points_j)) {
          d <- sqrt(sum((points_i[p1, ] - points_j[p2, ])^2))
          if (d < min_dist){
            min_dist <- d
          }
        }
      }
      distmat[i, j] <- min_dist
    }
  }
  return(distmat)
}

# clustdist_complete
#' max distance between two clusters
#'
#' max distance between two clusters
#'
#' @param features matrix with two columns: rt and ncarbons
#' @param atclust cluster of each feature at a given point
#'
#' @return matrix
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
clustdist_complete <- function(features, atclust) {
  n <- length(atclust)
  if (n < 2) return(matrix(NA, n, n))
  distmat <- matrix(NA, n, n)
  for (i in 2:n){
    for (j in 1:(i - 1)) {
      points_i <- features[atclust[[i]], , drop = FALSE]
      points_j <- features[atclust[[j]], , drop = FALSE]
      max_dist <- -Inf
      for (p1 in 1:nrow(points_i)) {
        for (p2 in 1:nrow(points_j)) {
          d <- sqrt(sum((points_i[p1, ] - points_j[p2, ])^2))
          if (d > max_dist){
            max_dist <- d
          }
        }
      }
      distmat[i, j] <- max_dist
    }
  }
  return(distmat)
}


# removedistantrtdatapoints
#' Remove distant data points for RT modelling
#'
#' Remove distant data points for RT modeling that are more than distC away from
#' the lipids located at the median of C
#'
#' @param rtdata rtdata
#' @param distC_filter numeric. minimum number of carbons from which lipids are considered distant from the rest
#' @param minSpecies_filter numeric. minimum number of lipids in a cluster from which the cluster is considered correct from the rest.
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
removedistantrtdatapoints <- function(rtdata, distC_filter = 5, minSpecies_filter = 5,
                                      typeanalysis = "filtering") {

  data <- combinefilterpredictions(rtdata, typeanalysis)

  for (lipidsubclass in unique(data$LipidSubclass)) {
    ss_index <- data$LipidSubclass == lipidsubclass
    data_class <- data[ss_index, , drop = FALSE]
    data_class <- data_class[data_class$Filter == "Kept", , drop = FALSE]

    if (nrow(data_class) > 0) {
      total_db <- unique(data_class$DB)[order(unique(data_class$DB))]

      # Remove distant data points
      for (db in total_db) { # for each DB serie
        c_values <- unique(data_class$C[data_class$DB == db])
        c_values <- c_values[order(c_values, decreasing = FALSE)]

        # Detect a break greater than distC
        index <- which(abs(diff(c_values)) > distC_filter)

        # Index refers to the positions of the C that have a distance greater than distC to the next C
        if (length(c_values) > 1 & length(index) > 0) { # If there is any C within that DB that is more than distC carbons away...

          # clusterize based on the breaks detected and assign reference cluster:
          c_clust <- rep(0, length(c_values) + 1)
          c_clust[index + 1] <- 1
          c_clust <- c_clust[1:length(c_values)]
          c_clust <- cumsum(c_clust) + 1
          dist_clust <- split(c_values, c_clust)
          median_clust <- median(data_class$C[data_class$DB == db])
          median_clust_general <- median(data_class$C)
          clust_median <- c_clust[which.min(abs(median_clust - c_values))]

          # Check which cluster has the fewest C atoms above the minimum required to avoid elimination
          # If all have fewer C atoms than minC, all would be eliminated except the one at the median
          # If any have a higher number, that one and the closest ones to the median cluster are kept
          left_clusters <- sort(unique(c_clust[c_clust < clust_median]))
          right_clusters <- sort(unique(c_clust[c_clust > clust_median]), decreasing = TRUE)

          c_values_outliers <- c()
          if (length(left_clusters) != 0) { # If there is a cluster to the left of the median cluster
            for (c in as.character(left_clusters)) { # For each cluster
              if (length(dist_clust[[c]]) < minSpecies_filter) { # If they have fewer than minC, the C to be eliminated are stored
                if (!any(abs(dist_clust[[c]] - median_clust_general) < distC_filter)) {
                  c_values_outliers <- c(c_values_outliers, dist_clust[[c]])
                }
              } else {
                break
              }
            }
          }

          if (length(right_clusters) != 0) { # If there is a cluster to the right of the median cluster
            for (c in as.character(right_clusters)) { # The clusters to the right of the median
              if (length(dist_clust[[c]]) < minSpecies_filter) {
                if (!any(abs(dist_clust[[c]] - median_clust_general) < distC_filter)) {
                  c_values_outliers <- c(c_values_outliers, dist_clust[[c]])
                }
              } else {
                break
              }
            }
          }

          to_remove <- data_class[data_class$DB == db & data_class$C %in%
                                        c_values_outliers, "index"]

          # data$Filter <- ifelse(data$index %in% index_distant,
          #                       "Removed distant data",
          #                       data$Filter)

          if(typeanalysis == "filtering") {
            rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                                   "Removed distant data",
                                                   rtdata$filtering$data$Filter)
          } else if (typeanalysis == "prediction") {
            rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                                  !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                                "Discarded - Distant data",
                                                                # "Discarded Filter",
                                                                rtdata$prediction$candidates$Filter_final)

            rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                           !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                         "Removed distant data",
                                                         rtdata$filtering$data$Filter_final)
          }


        }
      }
    }
  }

  return(rtdata)
}

# get_lipid_name_column
#' Get lipid annotation column
#'
#' Select the lipid annotation column used for chain-level
#' comparisons. Preference is given to 'lipidname_expanded'
#' when available.
#'
#' @param x data frame containing lipid annotations
#'
#' @return Character string with the selected column name
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
get_lipid_name_column <- function(x) {
  if ("lipidname_expanded" %in% colnames(x)) {
    return("lipidname_expanded")
  }

  if ("lipidname" %in% colnames(x)) {
    return("lipidname")
  }

  NA_character_
}

# has_defined_chain_position
#' Check whether FA positions are defined
#'
#' Determine whether a lipid annotation contains defined FA
#' positions. Returns TRUE only for annotations containing "/"
#' and not "_".
#'
#' @param lipidname lipid annotation
#'
#' @return Logical value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
has_defined_chain_position <- function(lipidname) {
  lipidname <- as.character(lipidname)

  grepl("/", lipidname, fixed = TRUE) &&
    !grepl("_", lipidname, fixed = TRUE)
}

# get_chain_key
#' Generate chain composition key
#'
#' Generate a standardized chain-level key from a lipid
#' annotation. For ether lipids, a double bond shift can be
#' applied to the ether-linked chain to compare plasmalogens
#' and ethers against their corresponding diacyl species.
#'
#' @param lipidname lipid annotation
#' @param db_shift number of double bonds removed from the
#'   ether-linked chain
#'
#' @return Character vector containing chain composition keys
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
get_chain_key <- function(lipidname, db_shift = 0) {

  lipidname <- as.character(lipidname)
  key <- rep(NA_character_, length(lipidname))

  has_chains <- grepl("\\([^()]+/[^()]+\\)", lipidname) &
    !grepl("_", lipidname, fixed = TRUE)

  if (!any(has_chains)) {
    return(key)
  }

  chains <- sub("^.*\\(([^()]*)\\).*$", "\\1", lipidname[has_chains])

  key[has_chains] <- vapply(chains, function(x) {

    z <- strsplit(x, "/", fixed = TRUE)[[1]]

    if (length(z) < 2) {
      return(NA_character_)
    }

    first_chain <- z[1]
    is_ether_chain <- grepl("^P-|^O-", first_chain)

    first_chain <- sub("^P-", "", first_chain)
    first_chain <- sub("^O-", "", first_chain)

    z[-1] <- sub("^P-", "", z[-1])
    z[-1] <- sub("^O-", "", z[-1])

    if (db_shift != 0 && is_ether_chain) {

      db1 <- suppressWarnings(as.integer(sub("^.*:", "", first_chain)))

      if (is.na(db1)) {
        return(NA_character_)
      }

      db1 <- db1 - db_shift

      if (db1 < 0) {
        return(NA_character_)
      }

      first_chain <- sub(":\\d+$", paste0(":", db1), first_chain)
    }

    z[1] <- first_chain

    paste(sort(z), collapse = "/")

  }, character(1))

  key
}

# get_reference_subset
#' Get reference lipids for RT comparison
#'
#' Select the reference lipids against which a query lipid is
#' compared. When chain-level information is available, only
#' lipids with the corresponding chain composition are retained.
#'
#' @param ss data frame containing lipids to compare
#' @param ss_query query lipid
#' @param ref_class reference lipid subclass
#' @param db_shift number of double bonds removed from the
#'   ether-linked chain
#'
#' @return Data frame containing reference lipids
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
get_reference_subset <- function(ss, ss_query, ref_class, db_shift = 0) {

  ss_ref <- ss[ss$LipidSubclass == ref_class, , drop = FALSE]

  if (!nrow(ss_ref)) {
    return(ss_ref)
  }

  lipid_col <- get_lipid_name_column(ss)

  if (is.na(lipid_col)) {
    return(ss_ref)
  }

  query_name <- ss_query[[lipid_col]]

  # If the query has no defined chain position, compare against all
  # reference lipids in the corresponding sum-composition group,
  # including reference lipids with chain-level annotations.
  if (!has_defined_chain_position(query_name)) {
    return(ss_ref)
  }

  # If the query has defined chains, compare only against matching
  # chain-level references. If the shifted chain key is invalid,
  # do not compare against anything.
  query_key <- get_chain_key(query_name, db_shift = db_shift)

  if (is.na(query_key)) {
    return(ss_ref[FALSE, , drop = FALSE])
  }

  ref_key <- get_chain_key(ss_ref[[lipid_col]], db_shift = 0)

  ss_ref[!is.na(ref_key) & ref_key == query_key, , drop = FALSE]
}

# get_removal_indices
#' Identify lipids violating RT coherence rules
#'
#' Evaluate RT relationships between query and reference lipid
#' subclasses and return lipids violating the expected RT order.
#'
#' @param ss data frame containing lipids to compare
#' @param query_class query lipid subclass
#' @param ref_class reference lipid subclass
#' @param direction comparison direction. Allowed values are
#'   "greater_equal" and "less_equal"
#' @param db_shift number of double bonds removed from the
#'   ether-linked chain
#'
#' @return Integer vector containing lipid indices to remove
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
get_removal_indices <- function(ss, query_class, ref_class, direction,
                                db_shift = 0) {

  ss_query_all <- ss[ss$LipidSubclass == query_class, , drop = FALSE]

  if (!nrow(ss_query_all)) {
    return(integer(0))
  }

  rem <- c()

  for (j in seq_len(nrow(ss_query_all))) {

    ss_query <- ss_query_all[j, , drop = FALSE]

    ss_ref <- get_reference_subset(
      ss = ss,
      ss_query = ss_query,
      ref_class = ref_class,
      db_shift = db_shift
    )

    ref_rt <- ss_ref$RTsec
    ref_rt <- ref_rt[!is.na(ref_rt)]

    if (!length(ref_rt) || is.na(ss_query$RTsec)) {
      next
    }

    if (direction == "greater_equal" && all(ss_query$RTsec >= ref_rt)) {
      rem <- c(rem, ss_query$index)
    }

    if (direction == "less_equal" && all(ss_query$RTsec <= ref_rt)) {
      rem <- c(rem, ss_query$index)
    }
  }

  rem
}

# has_chain_level_comparison
#' Check whether chain-level comparison is possible
#'
#' Determine whether chain-level comparison can be performed
#' between two lipid subclasses. Returns TRUE only when at
#' least one valid chain-level match exists.
#'
#' @param ss data frame containing lipids to compare
#' @param query_class query lipid subclass
#' @param ref_class reference lipid subclass
#' @param db_shift number of double bonds removed from the
#'   ether-linked chain
#'
#' @return Logical value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
has_chain_level_comparison <- function(ss, query_class, ref_class,
                                       db_shift = 0) {

  lipid_col <- get_lipid_name_column(ss)

  if (is.na(lipid_col)) {
    return(FALSE)
  }

  ss_query <- ss[ss$LipidSubclass == query_class, , drop = FALSE]
  ss_ref <- ss[ss$LipidSubclass == ref_class, , drop = FALSE]

  if (!nrow(ss_query) || !nrow(ss_ref)) {
    return(FALSE)
  }

  query_has_chains <- vapply(
    ss_query[[lipid_col]],
    has_defined_chain_position,
    logical(1)
  )

  if (!any(query_has_chains)) {
    return(FALSE)
  }

  query_keys <- get_chain_key(ss_query[[lipid_col]][query_has_chains],
                              db_shift = db_shift)
  ref_keys <- get_chain_key(ss_ref[[lipid_col]], db_shift = 0)

  any(!is.na(query_keys) & query_keys %in% ref_keys[!is.na(ref_keys)])
}

# checkplasmalogens
#' Check RT coherence between diacyl, ether and plasmalogen lipids
#'
#' Evaluate the expected retention-time order between diacyl,
#' alkyl-ether and plasmalogen lipid subclasses.
#'
#' The following rules are evaluated:
#' 1) PL-P(X:Y+1) < PL-O(X:Y+1) < PL(X:Y)
#' 2) PL-P(X:Y+2) < PL-O(X:Y+2) < PL(X:Y)
#' 3) PL(X:Y) < PL-P(X:Y) < PL-O(X:Y)
#'
#' When chain-level annotations are available, comparisons are
#' performed using the corresponding FA composition. Otherwise,
#' comparisons are performed using sum composition (C:DB).
#'
#' The PL-O versus PL-P comparison is only evaluated when
#' chain-level information is available.
#'
#' @param rtdata LipidRT object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return Updated LipidRT object
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
checkplasmalogens <- function(rtdata, typeanalysis = "filtering") {

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  lipidsubclasses <- c("PC", "PE", "PG", "PI", "PS", "PA",
                       "LPC", "LPE", "LPG", "LPI", "LPS", "LPA")

  # Rule 1: PL-P(X:Y+1) < PL-O(X:Y+1) < PL(X:Y)
  to_remove <- c()

  for (subclass in lipidsubclasses) {

    all_classes <- c(subclass, paste0(subclass, "-P"), paste0(subclass, "-O"))
    data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, ,
                            drop = FALSE]

    data_class$comparisonDB <- data_class$DB

    is_ether <- data_class$LipidSubclass %in% all_classes[2:3]
    data_class$comparisonDB[is_ether] <- data_class$DB[is_ether] - 1

    data_class <- data_class[data_class$comparisonDB >= 0, , drop = FALSE]

    index_list <- split(
      data_class$index,
      paste(data_class$C, data_class$comparisonDB, sep = "_")
    )

    index_list <- index_list[lengths(index_list) >= 2]

    index_list <- Filter(function(idx) {
      length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
    }, index_list)

    for (i in index_list) {

      ss <- data_class[data_class$index %in% i, , drop = FALSE]

      if (all(c(all_classes[1], all_classes[3]) %in% ss$LipidSubclass)) {
        to_remove <- c(
          to_remove,
          get_removal_indices(ss, all_classes[3], all_classes[1],
                              "greater_equal", db_shift = 1)
        )
      }

      if (all(c(all_classes[1], all_classes[2]) %in% ss$LipidSubclass)) {
        to_remove <- c(
          to_remove,
          get_removal_indices(ss, all_classes[2], all_classes[1],
                              "greater_equal", db_shift = 1)
        )
      }
    }
  }

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  # Rule 2: PL-P(X:Y+2) < PL-O(X:Y+2) < PL(X:Y)
  to_remove <- c()

  for (subclass in lipidsubclasses) {

    all_classes <- c(subclass, paste0(subclass, "-P"), paste0(subclass, "-O"))
    data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, ,
                            drop = FALSE]

    data_class$comparisonDB <- data_class$DB

    is_ether <- data_class$LipidSubclass %in% all_classes[2:3]
    data_class$comparisonDB[is_ether] <- data_class$DB[is_ether] - 2

    data_class <- data_class[data_class$comparisonDB >= 0, , drop = FALSE]

    index_list <- split(
      data_class$index,
      paste(data_class$C, data_class$comparisonDB, sep = "_")
    )

    index_list <- index_list[lengths(index_list) >= 2]

    index_list <- Filter(function(idx) {
      length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
    }, index_list)

    for (i in index_list) {

      ss <- data_class[data_class$index %in% i, , drop = FALSE]

      if (all(c(all_classes[1], all_classes[3]) %in% ss$LipidSubclass)) {
        to_remove <- c(
          to_remove,
          get_removal_indices(ss, all_classes[3], all_classes[1],
                              "greater_equal", db_shift = 2)
        )
      }

      if (all(c(all_classes[1], all_classes[2]) %in% ss$LipidSubclass)) {
        to_remove <- c(
          to_remove,
          get_removal_indices(ss, all_classes[2], all_classes[1],
                              "greater_equal", db_shift = 2)
        )
      }
    }
  }

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  # Rule 3: PL < PL-P < PL-O for the same C:DB composition.
  # The PL-O vs PL-P comparison is only applied when chain-level
  # information is available.
  to_remove <- c()

  for (subclass in lipidsubclasses) {

    all_classes <- c(subclass, paste0(subclass, "-P"), paste0(subclass, "-O"))
    data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, ,
                            drop = FALSE]

    index_list <- split(
      data_class$index,
      paste(data_class$C, data_class$DB, sep = "_")
    )

    index_list <- index_list[lengths(index_list) >= 2]

    index_list <- Filter(function(idx) {
      length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
    }, index_list)

    for (i in index_list) {

      ss <- data_class[data_class$index %in% i, , drop = FALSE]

      if (all(c(all_classes[1], all_classes[2]) %in% ss$LipidSubclass)) {
        to_remove <- c(
          to_remove,
          get_removal_indices(ss, all_classes[2], all_classes[1],
                              "less_equal", db_shift = 0)
        )
      }

      if (all(c(all_classes[1], all_classes[3]) %in% ss$LipidSubclass)) {
        to_remove <- c(
          to_remove,
          get_removal_indices(ss, all_classes[3], all_classes[1],
                              "less_equal", db_shift = 0)
        )
      }

      if (all(all_classes %in% ss$LipidSubclass) &&
          has_chain_level_comparison(ss, all_classes[3], all_classes[2],
                                     db_shift = 0)) {
        to_remove <- c(
          to_remove,
          get_removal_indices(ss, all_classes[3], all_classes[2],
                              "less_equal", db_shift = 0)
        )
      }
    }
  }

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }

  return(rtdata)
}

# checksphingolipids
#' Check elution order between sphingolipids
#'
#' Check elution order between sphingolipids
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checksphingolipids <- function(rtdata, typeanalysis = "filtering") {

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  to_remove <- c()

  for (suffix in c("_d", "_t")) {

    all_classes <- paste0(c("SHexCer", "HexCer", "Cer", "AHexCer", "SM"), suffix)

    data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, ]
    if (nrow(data_class) == 0) next

    data_class <- data_class[order(data_class$C, data_class$DB), ]

    index_list <- split(data_class$index, paste(data_class$C, data_class$DB, sep = "_"))

    index_list <- index_list[lengths(index_list) >= 2]
    index_list <- Filter(function(idx) {
      length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
    }, index_list)

    for (i in index_list) {

      ss <- data_class[data_class$index %in% i, ]

      rt_classes <- sapply(split(ss$RTsec, ss$LipidSubclass), median, na.rm = TRUE)

      # SHexCer_<d/t> < HexCer_<d/t>, remove SHexCer_<d/t>
      if (all(c(paste0("SHexCer", suffix), paste0("HexCer", suffix)) %in% ss$LipidSubclass)) {
        ss2 <- ss[ss$LipidSubclass == paste0("SHexCer", suffix), ]
        to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes[paste0("HexCer", suffix)]])
      }

      # HexCer_<d/t> < Cer_<d/t>, remove HexCer_<d/t>
      if (all(c(paste0("HexCer", suffix), paste0("Cer", suffix)) %in% ss$LipidSubclass)) {
        ss2 <- ss[ss$LipidSubclass == paste0("HexCer", suffix), ]
        to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes[paste0("Cer", suffix)]])
      }

      # Cer_<d/t> < AHexCer_<d/t>, remove AHexCer_<d/t> (Cer reference)
      if (all(c(paste0("Cer", suffix), paste0("AHexCer", suffix)) %in% ss$LipidSubclass)) {
        ss2 <- ss[ss$LipidSubclass == paste0("AHexCer", suffix), ]
        to_remove <- c(to_remove, ss2$index[ss2$RTsec < rt_classes[paste0("Cer", suffix)]])
      }

      if (all(c(paste0("Cer", suffix), paste0("SM", suffix)) %in% ss$LipidSubclass)) {
        ss2 <- ss[ss$LipidSubclass == paste0("Cer", suffix), ]
        to_remove <- c(to_remove, ss2$index[ss2$RTsec < rt_classes[paste0("SM", suffix)]])
      }
    }
  }

  # data$Filter <- ifelse(
  #   data$index %in% unique(to_remove),
  #   "RT coherence between subclasses",
  #   data$Filter
  # )

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }


  # data_kept <- data[!is.na(data$lipidname) &
  #                     data$lipidname != "" &
  #                     data$Filter == "Kept", ,
  #                   drop = FALSE]

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  to_remove <- c()

  # Rule 2: Expected order between -t/-d: t < d
  lipidsubclasses <- c("Cer", "SM", "HexCer", "SHexCer", "AHexCer", "Sph")

  for (subclass in lipidsubclasses) {

    all_classes <- paste0(subclass, c("_d", "_t"))

    data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, ]
    if (nrow(data_class) == 0) next

    data_class <- data_class[order(data_class$C, data_class$DB), ]

    index_list <- split(data_class$index, paste(data_class$C, data_class$DB, sep = "_"))

    index_list <- index_list[lengths(index_list) >= 2]
    index_list <- Filter(function(idx) {
      length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
    }, index_list)

    for (i in index_list) {

      ss <- data_class[data_class$index %in% i, ]

      rt_classes <- sapply(split(ss$RTsec, ss$LipidSubclass), median, na.rm = TRUE)

      if (all(all_classes %in% ss$LipidSubclass)) {
        ss2 <- ss[ss$LipidSubclass == all_classes[2], ]
        to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes[all_classes[1]]])
      }

    }
  }

  # data$Filter <- ifelse(
  #   data$index %in% unique(to_remove),
  #   "RT coherence between subclasses",
  #   data$Filter
  # )

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }

  # data_kept <- data[!is.na(data$lipidname) &
  #                     data$lipidname != "" &
  #                     data$Filter == "Kept", ,
  #                   drop = FALSE]

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  to_remove <- c()

  # Rule 1: Expected order between -t and -d with shifted DB:
  # subclass_t(X:Y) < subclass_d(X:Y+1)
  lipidsubclasses <- c("Cer", "SM", "HexCer", "SHexCer", "AHexCer", "Sph")

  for (subclass in lipidsubclasses) {

    all_classes <- paste0(subclass, c("_d", "_t"))

    data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, , drop = FALSE]
    if (nrow(data_class) == 0) next

    data_class <- data_class[order(data_class$C, data_class$DB), , drop = FALSE]

    data_class$comparisonDB <- data_class$DB

    is_d <- data_class$LipidSubclass == all_classes[1]

    data_class$comparisonDB[is_d] <-
      data_class$DB[is_d] - 1

    data_class <- data_class[data_class$comparisonDB >= 0, , drop = FALSE]

    index_list <- split(
      data_class$index,
      paste(data_class$C, data_class$comparisonDB, sep = "_")
    )

    index_list <- index_list[lengths(index_list) >= 2]

    index_list <- Filter(function(idx) {
      length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
    }, index_list)

    for (i in index_list) {

      ss <- data_class[data_class$index %in% i, , drop = FALSE]

      rt_classes <- sapply(
        split(ss$RTsec, ss$LipidSubclass),
        median,
        na.rm = TRUE
      )

      #  subclass_t(X:Y) should not elute after subclass_d(X:Y+1)
      if (all(all_classes %in% ss$LipidSubclass)) {
        ss2 <- ss[ss$LipidSubclass == all_classes[2], , drop = FALSE]
        rem <- ss2$index[which(ss2$RTsec >= rt_classes[all_classes[1]])]
        to_remove <- c(to_remove, rem)
      }
    }
  }

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }

  return(rtdata)
}



# checkglycolipids_DG
#' Check elution order bertween DG, MGDG and DGDG
#'
#'  Check elution order bertween DG, MGDG and DGDG
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkglycolipids_DG <- function(rtdata, typeanalysis = "filtering") {

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  to_remove <- c()

  all_classes <- c("DGDG", "MGDG", "DG")
  data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, ]

  data_class <- data_class[order(data_class$C, data_class$DB), ]

  index_list <- split(data_class$index, paste(data_class$C, data_class$DB, sep = "_"))

  index_list <- index_list[lengths(index_list) >= 2]
  index_list <- Filter(function(idx) {
    length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
  }, index_list)

  for (i in index_list) {

    ss <- data_class[data_class$index %in% i, ]

    rt_classes <- sapply(split(ss$RTsec, ss$LipidSubclass), median,  na.rm = TRUE)


    # DGDG < MGDG, remove DGDG
    if (all(c("DGDG", "MGDG") %in% ss$LipidSubclass)) {
      ss2 <- ss[ss$LipidSubclass == "DGDG", ]
      to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes["MGDG"]])
    }

    # MGDG < DGDG, remove DGDG
    if (all(c("DGDG", "DG") %in% ss$LipidSubclass)) {
      ss2 <- ss[ss$LipidSubclass == "DGDG", ]
      to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes["DG"]])
    }

    # DG < MGDG, remove MGDG (DG reference)
    if (all(c("MGDG", "DG") %in% ss$LipidSubclass)) {
      ss2 <- ss[ss$LipidSubclass == "MGDG", ]
      to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes["DG"]])
    }

  }

  # data$Filter <- ifelse(data$index %in% to_remove, "RT coherence between subclasses", data$Filter)

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }


  return(rtdata)
}



# checkglycolipids_DG
#' Check elution order bertween DG, MGDG and DGDG
#'
#'  Check elution order bertween DG, MGDG and DGDG
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkglycolipids_DG <- function(rtdata, typeanalysis = "filtering") {

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  to_remove <- c()

  all_classes <- c("DGDG", "MGDG", "DG")
  data_class <- data_kept[data_kept$LipidSubclass %in% all_classes, ]

  data_class <- data_class[order(data_class$C, data_class$DB), ]

  index_list <- split(data_class$index, paste(data_class$C, data_class$DB, sep = "_"))

  index_list <- index_list[lengths(index_list) >= 2]
  index_list <- Filter(function(idx) {
    length(unique(data_class$LipidSubclass[data_class$index %in% idx])) > 1
  }, index_list)

  for (i in index_list) {

    ss <- data_class[data_class$index %in% i, ]

    rt_classes <- sapply(split(ss$RTsec, ss$LipidSubclass), median,  na.rm = TRUE)


    # DGDG < MGDG, remove DGDG
    if (all(c("DGDG", "MGDG") %in% ss$LipidSubclass)) {
      ss2 <- ss[ss$LipidSubclass == "DGDG", ]
      to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes["MGDG"]])
    }

    # MGDG < DGDG, remove DGDG
    if (all(c("DGDG", "DG") %in% ss$LipidSubclass)) {
      ss2 <- ss[ss$LipidSubclass == "DGDG", ]
      to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes["DG"]])
    }

    # DG < MGDG, remove MGDG (DG reference)
    if (all(c("MGDG", "DG") %in% ss$LipidSubclass)) {
      ss2 <- ss[ss$LipidSubclass == "MGDG", ]
      to_remove <- c(to_remove, ss2$index[ss2$RTsec > rt_classes["DG"]])
    }

  }

  # data$Filter <- ifelse(data$index %in% to_remove, "RT coherence between subclasses", data$Filter)

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "RT coherence between subclasses",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - RT coherence between subclasses",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "RT coherence between subclasses",
                                                 rtdata$filtering$data$Filter_final)
  }


  return(rtdata)
}



# checkRTclasses
#' Check for RT coherence between classes
#'
#' Check for RT coherence between classes
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
checkRTclasses <- function(rtdata, typeanalysis = "filtering") {

  # Check elution order bertween diacyl an plasmalogens of PL
  rtdata <- checkplasmalogens(rtdata, typeanalysis)

  # Check elution order bertween Cer and HexCer
  rtdata <- checksphingolipids(rtdata, typeanalysis)

  # Check elution order bertween DG, MGDG and DGDG
  rtdata <- checkglycolipids_DG(rtdata, typeanalysis)


  return(rtdata)
}


# firstfilterbyrtrules
#' Remove those annotation that do not match the specified RT rules
#'
#' Remove those annotation that do not match the specified RT rules based on the
#' expected elution order.
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
firstfilterbyrtrules <- function(rtdata, typeanalysis = "filtering") {

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  for (lipidsubclass in unique(data_kept$LipidSubclass)) {
    data_class <- data_kept[data_kept$LipidSubclass == lipidsubclass, , drop = FALSE]

    # The data is filtered the same number of times as length(total_db)
    # randomizing the order of DB series
    set.seed(2024)
    to_remove <- c()
    total_db <- sort(unique(data_class$DB))
    for (t in 1:length(total_db)){
      db_values <- sample(total_db)
      # apply RT rules following the randomized DB series order
      filter_db <- filterdata(data_class, db_values)
      to_remove <- c(to_remove, as.numeric(filter_db))
    }
    # finally, remove only ids that have been removed the 80% of the times
    threshold <- ceiling(length(total_db) * 0.8)
    index_to_remove <- as.numeric(names(which(table(to_remove) >= threshold)))

    # data$Filter <- ifelse(data$index %in% index_to_remove, "Removed by RT rules",
    #                       data$Filter)

    if(typeanalysis == "filtering") {
      rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% index_to_remove,
                                             "Removed by RT rules",
                                             rtdata$filtering$data$Filter)
    } else if (typeanalysis == "prediction") {
      rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% index_to_remove &
                                                            !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                          "Discarded - Removed by RT rules",
                                                          rtdata$prediction$candidates$Filter_final)

      rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% index_to_remove &
                                                     !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                   "Removed by RT rules",
                                                   rtdata$filtering$data$Filter_final)
    }
  }

  return(rtdata)
}

# secondfilterbyrtrules
#' Repeat rt rules-based filtering with more strict parameters
#'
#' Repeat rt rules-based filtering with more strict parameters
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
secondfilterbyrtrules <- function(rtdata, typeanalysis = "filtering") {

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  for (lipidsubclass in unique(data_kept$LipidSubclass)) {
    data_class <- data_kept[data_kept$LipidSubclass == lipidsubclass, , drop = FALSE]

    # Set the DB series order based on the number of records for each serie
    db_values <- as.numeric(names(sort(table(data_class$DB))))
    # apply RT rules following the DB series order
    index_to_remove <- filterdata(data_class, db_values)

    # data$Filter <- ifelse(data$index %in% index_to_remove, "Removed by RT rules",
    #                       data$Filter)

    if(typeanalysis == "filtering") {
      rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% index_to_remove,
                                             "Removed by RT rules",
                                             rtdata$filtering$data$Filter)
    } else if (typeanalysis == "prediction") {
      rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% index_to_remove &
                                                            !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                          "Discarded - Removed by RT rules",
                                                          rtdata$prediction$candidates$Filter_final)

      rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% index_to_remove &
                                                     !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                   "Removed by RT rules",
                                                   rtdata$filtering$data$Filter_final)
    }
  }

  return(rtdata)
}

# filterbyrtrules
#' Filtering of each lipid class according to its RT pattern based on its
#' composition of C and DB
#'
#' Filtering of each lipid class according to its RT pattern based on its
#' composition of C and DB
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return rtdata.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
filterbyrtrules <- function(rtdata, typeanalysis = "filtering") {

  # PHASE 1: Strict filtering (elimination of those excluded 80% of the time)
  rtdata <- firstfilterbyrtrules(rtdata)

  # PHASE 2: Less strict filtering (ascending order of data within each DB)
  rtdata <- secondfilterbyrtrules(rtdata)

  return(rtdata)
}

# safe_min
#' min value of a numeric vector returning NA if vector has no elements
#'
#' min value of a numeric vector returning NA if vector has no elements
#'
#' @param x numeric vector
#'
#' @return numeric
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
safe_min <- function(x) {
  if (length(x[!is.na(x)]) > 0) min(x, na.rm = TRUE) else NA
}

# safe_max
#' max value of a numeric vector returning NA if vector has no elements
#'
#' max value of a numeric vector returning NA if vector has no elements
#'
#' @param x numeric vector
#'
#' @return numeric
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
safe_max <- function(x) {
  if (length(x[!is.na(x)]) > 0) max(x, na.rm = TRUE) else NA
}

# filtercondition
#' The conditions that a lipid must meet according to its RT
#'
#' The conditions that a lipid must meet according to its RT in order not to be
#' eliminated during the filtering
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#' @param db_value numeric. value of DB on which the filtering will be performed
#' @param c_values vector. possible values of C within the DB on which filtering will be performed
#'
#' @return to_remove
#' vector containing the indices of the data points that are removed during filtering
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
filtercondition <- function(data, db_value, c_values) {

  to_remove <- c()
  for (c_value in c_values) {
    # Filter the data based on the DB and C conditions
    subset_data <- data[data$DB == db_value & data$C == c_value, ]

    # For each C value, determine <C and >C, and <DB and >DB, to enable comparisons
    # Retrieve the C value located to the left of c_value with the same DB (<C, DB)
    c_right <- safe_min(data$C[data$DB == db_value & data$C > c_value])
    # Retrieve the C value located to the right of c_value with the same DB (>C, DB)
    c_left <- safe_max(data$C[data$DB == db_value & data$C < c_value])
    # Retrieve the DB value that is below db_value with the same c_value (C, <DB)
    db_right <- safe_min(data$DB[data$C == c_value & data$DB > db_value])
    # Retrieve the DB value that is above db_value with the same c_value (C, >DB)
    db_left <- safe_max(data$DB[data$C == c_value & data$DB < db_value])

    # Calculate the RT medians for each of the C and DB values used to compare the lipid of interest
    if (!is.na(c_right)) {
      medianC_right <- median(data$RTsec[data$C == c_right & data$DB == db_value])
    }
    if (!is.na(c_left)) {
      medianC_left <- median(data$RTsec[data$C == c_left & data$DB == db_value])
    }
    if (!is.na(db_right)) {
      medianDB_right <- median(data$RTsec[data$C == c_value & data$DB == db_right])
    }
    if (!is.na(db_left)) {
      medianDB_left <- median(data$RTsec[data$C == c_value & data$DB == db_left])
    }

    # For a lipid to be retained, its RT must meet the following conditions:
    # The RT of the lipid (C, DB) must be > than the RT of the median of lipids (< C, = DB)
    # The RT of the lipid (C, DB) must be > than the RT of the median of lipids (= C, > DB)
    # The RT of the lipid (C, DB) must be < than the RT of the median of lipids (> C, = DB)
    # The RT of the lipid (C, DB) must be < than the RT of the median of lipids (= C, < DB)

    # For each group of lipids with a given c_value and db_value, the conditions
    # above are checked by comparing them with the surrounding lipids.

    # The code considers all possible cases where lipids may be missing for comparison.
    # Different possibilities for each data point and the criteria to follow.
    if (is.na(db_left) & is.na(db_right) & is.na(c_left) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianC_right]
    } else if (is.na(db_left) & is.na(db_right) & is.na(c_right) & !is.na(c_left)) {
      index_rem <- subset_data$index[subset_data$RTsec < medianC_left]
    } else if (is.na(db_right) & is.na(c_left) & is.na(c_right) & !is.na(db_left)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianDB_left]
    } else if (is.na(db_left) & is.na(c_left) & is.na(c_right) & !is.na(db_right)) {
      index_rem <- subset_data$index[subset_data$RTsec < medianDB_right]
    } else if (is.na(db_left) & is.na(db_right) & !is.na(c_left) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianC_right |
                                       subset_data$RTsec < medianC_left]
    } else if (is.na(db_right) & is.na(c_left) & !is.na(db_left) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianC_right |
                                       subset_data$RTsec > medianDB_left]
    } else if (is.na(db_left) & is.na(c_left) & !is.na(db_right) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianC_right |
                                       subset_data$RTsec < medianDB_right]
    } else if (is.na(db_right) & is.na(c_right) & !is.na(db_left) & !is.na(c_left)) {
      index_rem <- subset_data$index[subset_data$RTsec < medianC_left |
                                       subset_data$RTsec > medianDB_left]
    } else if (is.na(db_left) & is.na(c_right) & !is.na(db_right) & !is.na(c_left)) {
      index_rem <- subset_data$index[subset_data$RTsec < medianC_left |
                                       subset_data$RTsec < medianDB_right]
    } else if (is.na(c_left) & is.na(c_right) & !is.na(db_left) & !is.na(db_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianDB_left |
                                       subset_data$RTsec < medianDB_right]
    } else if (is.na(db_left) & !is.na(db_right) & !is.na(c_left) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec < medianDB_right |
                                       subset_data$RTsec < medianC_left |
                                       subset_data$RTsec > medianC_right]
    } else if (is.na(db_right) & !is.na(db_left) & !is.na(c_left) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianDB_left |
                                       subset_data$RTsec < medianC_left |
                                       subset_data$RTsec > medianC_right]
    } else if (is.na(c_left) & !is.na(db_left) & !is.na(db_right) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianDB_left |
                                       subset_data$RTsec < medianDB_right |
                                       subset_data$RTsec > medianC_right]
    } else if (is.na(c_right) & !is.na(db_left) & !is.na(db_right) & !is.na(c_left)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianDB_left |
                                       subset_data$RTsec < medianDB_right |
                                       subset_data$RTsec < medianC_left]
    } else if (!is.na(db_left) & !is.na(db_right) & !is.na(c_left) & !is.na(c_right)) {
      index_rem <- subset_data$index[subset_data$RTsec > medianDB_left |
                                       subset_data$RTsec < medianDB_right |
                                       subset_data$RTsec < medianC_left |
                                       subset_data$RTsec > medianC_right]
    } else {
      index_rem <- c()
    }

    to_remove <- unique(c(to_remove, index_rem))

    # Remove those data points for the next iteration
    data <- data[!data$index %in% index_rem, ]

    # Other conditions being considered (< DB, > C)
    # Retrieve the DB value that is above db_value regardless of c_value (<DB)
    db_c_left <- safe_max(data$DB[which(data$DB < db_value)]) # < DB
    if (!is.na(db_c_left)) {
      # If there is a <DB, retrieve the highest C value greater than c_value
      # within that <DB (the top-right point relative to c_value, db_value)
      c_db_left <- safe_min(data$C[data$DB == db_c_left & data$C > c_value]) # > C
      if (!is.na(c_db_left)) {
        # If there is a lipid with < DB and > C compared to C:DB, retrieve its RT median and compare
        # The RT of the lipid (C, DB) must be < than the RT of the median of lipids (> C, < DB)
        medianRT <- median(data$RTsec[data$C == c_db_left & data$DB == db_c_left])
        index_rem <- subset_data$index[subset_data$RTsec > medianRT]
        if (length(index_rem) != 0) {
          to_remove <- unique(c(to_remove, index_rem))
          data <- subset(data, !(data$index %in% index_rem))
        }
      }
    }

    # Other conditions being considered (> DB, < C)
    # Retrieve the DB value that is below db_value regardless of c_value (>DB)
    db_c_right <- safe_min(data$DB[which(data$DB > db_value)]) # > DB
    if (!is.na(db_c_left)) {
      # If there is a > DB, retrieve the lowest C value less than c_value within
      # that > DB (the bottom-left point relative to c_value, db_value)
      c_db_right <- safe_max(data$C[data$DB == db_c_right & data$C < c_value]) # < C
      if (!is.na(c_db_right)) {
        # If there is a lipid with > DB and < C compared to C:DB, retrieve its RT median and compare
        # The RT of the lipid (C, DB) must be > than the RT of the median of lipids (< C, > DB)
        medianRT <- median(data$RTsec[data$C == c_db_right & data$DB == db_c_right])
        index_rem <- subset_data$index[subset_data$RTsec < medianRT]
        if (length(index_rem) != 0) {
          to_remove <- unique(c(to_remove, index_rem))
          data <- subset(data, !(data$index %in% index_rem))
        }
      }
    }
  }

  return(to_remove)
}

# filterdata
#' Filtering annotations based on the RT of lipids.
#'
#' Filtering annotations using the number of carbons and double bonds of the
#' lipids belonging to each class to check if the annotated species follow a
#' pattern consistent with their composition, eliminating those that do not.
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#' @param db_values vector. All possible values of DB on this data
#'
#' @return vector containing the indices of the data points that are removed
#' during filtering
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
filterdata <- function(data, db_values) {

  total_removed <- c()

  # Data belonging to each DB group are filtered
  for (db_value in db_values) {

    # set data filtering order
    c_values <- as.numeric(names(sort(table(data$C[data$DB == db_value]))))

    # Filter data
    to_remove <- filtercondition(data, db_value, c_values)

    total_removed <- c(total_removed, to_remove)

    # Remove those data points for the next iteration
    data <- subset(data, !(data$index %in% to_remove))
  }

  return(total_removed)

}


# alldifferences
#' calculate difference between all the elements of a vector
#'
#' calculate difference between all the elements of a vector
#'
#' @param vec numeric vector
#'
#' @return vector
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
alldifferences <- function(vec) {
  if (length(vec) < 2) return(numeric(0))
  combs <- combn(vec, 2)
  abs(combs[1, ] - combs[2, ])
}

# checkoutliers
#' Check outliers to remove incorrect annotations
#'
#' Check outliers to remove incorrect annotations
#'
#' @param rtdata rtdata object
#' @param typeanalysis type of analysis. Either "filtering" or "prediction"
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
checkoutliers <- function(rtdata, typeanalysis = "filtering") {

  # data_kept <- data[!is.na(data$lipidname) &
  #                     data$lipidname != "" &
  #                     data$Filter == "Kept", ,
  #                   drop = FALSE]

  data_kept <- combinefilterpredictions(rtdata, typeanalysis)

  # Remove outliers
  to_remove <- c()
  lipidsubclass <- unique(data_kept$LipidSublass)

  for (class in lipidsubclass) {
    data_class <- data_kept[data_kept$LipidSubclass == class, c("index", "RTsec", "C", "DB")]

    clusters <- clust_insourcefrags(data_class$RTsec, data_class$C,
                                    scale.data = TRUE,
                                    rt_weight = 2,
                                    nc_weight = 0.5,
                                    dist_percentile = 0.75,
                                    linkage = "complete")

    nclusters <- length(unique(clusters))

    data_class$Cluster <- as.numeric(as.factor(clusters))

    perc <- table(data_class$Cluster)/length(clusters)

    ############################################################################
    # consider outliers if:
    # 1) cluster has n < 3 and covers less than 10% of the total class
    # iremove <- data_class$index[data_class$Cluster %in% which(table(data_class$Cluster) < 3) &
    #                               data_class$Cluster %in% which(perc < 0.1)]
    iremove <- data_class$index[data_class$Cluster %in% which(perc < 0.1)]

    # 2) rt is far enough
    if (length(iremove) > 0) {
      for (i in iremove){
        dists <- abs(data_class$RTsec[data_class$index == i] -
                       data_class$RTsec[data_class$index != i])
        all_dists <- alldifferences(data_class$RTsec[data_class$index != i])
        if (min(dists) < quantile(all_dists, probs = c(0.9))){
          iremove <- iremove[iremove != i]
        }
      }
    }

    rem <- data_class$index[data_class$index %in% iremove]
    to_remove <- append(to_remove, rem)
  }

  # data$Filter <- ifelse(data$index %in% to_remove, "Outlier", data$Filter)

  if(typeanalysis == "filtering") {
    rtdata$filtering$data$Filter <- ifelse(rtdata$filtering$data$index %in% to_remove,
                                           "Outlier",
                                           rtdata$filtering$data$Filter)
  } else if (typeanalysis == "prediction") {
    rtdata$prediction$candidates$Filter_final <- ifelse(rtdata$prediction$candidates$index %in% to_remove &
                                                          !(rtdata$prediction$candidates$Filter_final %in% c("Kept manually by the user")),
                                                        "Discarded - Outlier",
                                                        rtdata$prediction$candidates$Filter_final)

    rtdata$filtering$data$Filter_final <- ifelse(rtdata$filtering$data$index %in% to_remove &
                                                   !(rtdata$filtering$data$Filter_final %in% c("Kept manually by the user")),
                                                 "Outlier",
                                                 rtdata$filtering$data$Filter_final)
  }

  return(rtdata)
}

# checkremovedids
#' Predicts the retention time of lipids from the data and retrieves those that
#' were removed during filtering if their retention time is consistent with the
#' observed retention time for their lipid class.
#'
#' Predicts the retention time of lipids from the data and retrieves those that
#' were removed during filtering if their retention time is consistent with the
#' observed retention time for their lipid class.
#'
#' @param rtdata rtdata object.
#' @param distribution character. Data distribution for RT modelling. It can be "gamma" or "normal".
#' @param rtscorethr numeric. Minimum Score_rtpred value from which a prediction will be considered valid.
#'
#' @return A list containing: data, predictions, models
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkremovedids <- function(rtdata, distribution, rtscorethr) {

  # Annotated peaklist
  data <- rtdata$filtering$data

  # Modelling the retention time
  rtdata <- rtmodelling(rtdata = rtdata,
                        distribution = distribution,
                        classes = NULL,
                        typeanalysis = "filtering")

  # Predict retention time of lipids in order to recover incorrectly removed data
  rtdata <- predictannotations(rtdata = rtdata,
                               rtscorethr = rtscorethr,
                               typeanalysis = "filtering")

  return(rtdata)
}

# filteringgraph
#' Plot the filtered data for each class of lipids
#'
#' Plot the filtered data for each class of lipids
#'
#' @param rtdata rtdata object.
#' @param classes classes
#'
#' @return plots
#' list where each element is a plot of the filtering for each lipid class
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
filteringgraph <- function(rtdata, classes = NULL) {

  data <- rtdata$filtering$data
  # data <- data[!(data$Filter %in% c("Ignore duplicated", "Not pass QC",
  #                                     "Unsupported class")), , drop = FALSE]
  data <- data[!(data$Filter %in% c("Ignore duplicated", "Not pass QC",
                                    "Unsupported class")) &
                 !is.na(data$RTsec), , drop = FALSE]

  # All valid lipid classes in the dataset
  classes_validas <- unique(na.omit(data$LipidSubclass[data$LipidSubclass != ""]))
  if (is.null(classes)) {
    classes <- classes_validas
  } else {
    classes <- classes[classes %in% classes_validas]
  }


  # predictions <- rtdata$filtering$predictions
  models <- rtdata$filtering$models

  plots <- c()

  # Define levels manually
  shape_levels <- c("Kept" = 16,
                    "Kept manually by the user" = 17,
                    "Recovered with RT modelling" = 15,
                    "Removed with RT modelling" = 0,
                    "In-source fragments" = 8,
                    "Removed distant data" = 5,
                    "Confidence level among isobars" = 4,
                    "Removed by RT rules" = 1,
                    "Removed manually by the user" = 2,
                    "Outlier" = 6,
                    # "RT coherence between plasmalogens" = 3
                    "RT coherence between subclasses" = 3)

  # Plot the filtered data for each class of lipids
  for (lipidsubclass in classes) {
    data_class <- data[data$LipidSubclass == lipidsubclass, , drop = FALSE]

    if("Filter_final" %in% colnames(data_class)) {
      data_class$checkfilter <-data_class$Filter_final
    }

    # information about the classes used when building the model to include it in the subtitle of the graph
    lipid_subclasses <- NULL

    if (!is.null(models[[lipidsubclass]]$model$data) && !is.null(models[[lipidsubclass]]$model$data$LipidSubclass)) {
      lipid_subclasses <- unique(models[[lipidsubclass]]$model$data$LipidSubclass)
    } else if (!is.null(models[[lipidsubclass]]$model$model) && !is.null(models[[lipidsubclass]]$model$model$LipidSubclass)) {
      lipid_subclasses <- unique(models[[lipidsubclass]]$model$model$LipidSubclass)
    } else {
    }


    subtitle_text <- character(0)
    if (!is.null(lipid_subclasses)) {
      if (length(lipid_subclasses) == 1) {
        subtitle_text <- paste("This model was built using the lipid class:", lipid_subclasses)
      } else {
        subtitle_text <- paste("This model was built using the following lipid classes:\n   -",
                               paste0(lipid_subclasses, collapse = ", "))
      }
    } else {
      subtitle_text <- "No model for this lipid class"
    }


    data_class$Filter.1 <- ifelse(data_class$checkfilter %in%
                                    c("Kept manually by the user",
                                      "Removed manually by the user"),
                                  data_class$checkfilter, data_class$Filter)



    # Check levels
    data_class$Filter.1 <- factor(data_class$Filter.1, levels = names(shape_levels),
                                  exclude = NULL)
    data_class$checkfilter <- factor(data_class$checkfilter, levels = names(shape_levels),
                                     exclude = NULL)


    min_c <- min(unique(data_class$C)) - 1
    max_c <- max(unique(data_class$C)) + 1
    min_c <- ifelse(min_c %% 2 == 0, min_c - 2, min_c - 1)  # If it is even, subtract 2. If it is odd, subtract 1.
    max_c <- ifelse(max_c %% 2 == 0, max_c + 2, max_c + 1)  # If it is even, add 2. If it is odd, add 1.

    if (all(is.na(data_class$RTsec))) {
      warning("RTsec contains only NA values.")
      min_rt <- NA
      max_rt <- NA

    } else {
      min_rt <- floor(min(data_class$RTsec, na.rm = TRUE))
      if (min_rt < 0) min_rt <- 0
      max_rt <- ceiling(max(data_class$RTsec, na.rm = TRUE))
    }

    # check if prediction layer is required
    prediction_layer <- predictionlayergraph(rtdata = rtdata,
                                             lipidsubclass = lipidsubclass,
                                             min_c = min_c,
                                             max_c = max_c,
                                             min_rt = min_rt,
                                             max_rt = max_rt,
                                             alpha = 0.6,
                                             typeanalysis = "filtering")

    # Graph1: Show the data that is kept (filled circles) and those that
    # will be eliminated after filtering
    C <- RTsec <- DB <- Filter.1 <- RTsec_rtpred <- checkfilter <- comments <- c()
    # Build plot
    graph1 <- ggplot() +
      prediction_layer +
      geom_point(data = data_class,
                 aes(x = C,
                     y = RTsec,
                     color = factor(DB),
                     shape = factor(Filter.1)),
                 size = 2.5, stroke = 1.5) +
      labs(title = paste0("Filtering for '", lipidsubclass, "' class"),
           subtitle = subtitle_text,
           x = "Carbon number", y = "Retention time (s)", shape = "Filtered data") +
      theme_minimal() +
      theme(plot.subtitle = element_text(color = "#5A8D78")) +
      scale_color_manual(name = "Double \nbonds",
                         values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                    "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                    "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                    "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                    "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                    "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                    "#D9FF68",
                                    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                    "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                    "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                    "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                    "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                    "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                    "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                    "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                    "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                    "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                    "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                    "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                    "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                    "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                    "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                    "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                    "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                         breaks = seq(0:24)-1) +
      scale_shape_manual(values = shape_levels) +
      scale_x_continuous(breaks = seq(min_c, max_c, by = 2)) +
      scale_y_continuous()


    if ("comments" %in% colnames(data_class) &&
      any(!is.na(data_class$comments) & nzchar(trimws(data_class$comments)))) {
      graph1 <- graph1 +

        geom_point(
          data = subset(data_class, !is.na(comments) & comments != ""),
          aes(C, RTsec),
          shape = 22,
          size = 5,
          fill = "black",
          alpha = 0.25
        )
    }

    graph1 <- graph1 +
      coord_cartesian(
        xlim = c(min_c, max_c),
        ylim = c(min_rt, max_rt)
      )

    # Graph2: Show the data that is kept (filled circles)
    filtered_data <- data_class[data_class$checkfilter %in% c("Kept",
                                                              "Kept manually by the user",
                                                              "Recovered with RT modelling"), , drop = FALSE]

    kept_layer <- list()
    if (nrow(filtered_data) != 0) {
      kept_layer <- list(
        geom_point(data = filtered_data,
                   aes(x = C,
                       y = RTsec,
                       color = factor(DB),
                       shape = factor(checkfilter)),
                   size = 2.5, stroke = 1.5),
        scale_color_manual(name = "Double \nbonds",
                           values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                      "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                      "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                      "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                      "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                      "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                      "#D9FF68",
                                      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                      "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                      "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                      "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                      "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                      "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                      "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                      "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                      "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                      "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                      "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                      "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                      "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                      "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                      "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                      "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                      "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                      "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                           breaks = seq(0:24)-1),
        scale_shape_manual(values = shape_levels)
      )
    }

    # Build plot
    graph2 <- ggplot() +
      kept_layer +
      # labs(title = "Data kept according to the RT", x = "Carbon number", y = "", shape = "Data kept") +
      labs(title = "Data kept according to the RT", x = "Carbon number", y = "") +
      theme_minimal() +
      guides(color = "none") +
      scale_x_continuous(limits = c(min_c, max_c), breaks = seq(min_c, max_c, by = 2)) +
      scale_y_continuous(limits = c(min_rt, max_rt))

    if ("comments" %in% colnames(filtered_data) &&
        any(!is.na(filtered_data$comments) & nzchar(trimws(filtered_data$comments)))) {
      graph2 <- graph2 +

        geom_point(
          data = subset(filtered_data, !is.na(comments) & comments != ""),
          aes(C, RTsec),
          shape = 22,
          size = 5,
          fill = "black",
          alpha = 0.25
        )
    }

    # save plot
    grDevices::pdf(NULL) # use a pdf NULL device to save plots to an object
    grDevices::dev.control(displaylist = "enable")

    graph3 <- graph1 + graph2 +
      patchwork::plot_layout(ncol = 2, nrow = 1, guides = "collect")
    print(graph3)

    plots[[lipidsubclass]] <- grDevices::recordPlot() # save plot
    invisible(grDevices::dev.off()) # close pdf NULL device
  }

  rtdata$filtering$plots <- plots

  return(rtdata)
}



# filteringcombinedgraph
#' Plot the filtered data for each class of lipids
#'
#' Plot the filtered data for each class of lipids
#'
#' @param rtdata rtdata object.
#' @param classes classes
#'
#' @return plots
#' list where each element is a plot of the filtering for each lipid class
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
filteringcombinedgraph <- function(rtdata, classes = NULL) {

  data <- rtdata$filtering$data

  data <- data[!(data$Filter %in% c("Ignore duplicated",
                                    "Not pass QC",
                                    "Unsupported class")), , drop = FALSE]

  classes_validas <- unique(na.omit(data$LipidSubclass[data$LipidSubclass != ""]))

  if (is.null(classes)) {
    lipidsubclasses <- classes_validas
  } else {
    lipidsubclasses <- classes[classes %in% classes_validas]
  }

  if (length(lipidsubclasses) == 0) {
    stop("No valid classes.")
  }

  if (length(lipidsubclasses) > 7) {
    stop("Maximum 7 classes.")
  }

  lipidsubclasses_combined <- paste0(lipidsubclasses, collapse = "|")

  class_shapes <- c(0,1,2,3,5,6,7)

  shape_map <- setNames(class_shapes[seq_along(lipidsubclasses)], lipidsubclasses)

  kept_levels <- c("Kept",
                   "Kept manually by the user",
                   "Recovered with RT modelling")

  removed_levels <- c("Removed with RT modelling",
                      "In-source fragments",
                      "Removed distant data",
                      "Confidence level among isobars",
                      "Removed by RT rules",
                      "Removed manually by the user",
                      "Outlier",
                      "RT coherence between subclasses")

  data_classes <- data[data$LipidSubclass %in% lipidsubclasses, , drop = FALSE]

  data_classes$LipidSubclass <- factor(data_classes$LipidSubclass, levels = lipidsubclasses)

  data_classes$Filter.1 <- ifelse(data_classes$checkfilter %in% c("Kept manually by the user",
                                                                  "Removed manually by the user"),
                                  data_classes$checkfilter, data_classes$Filter)

  min_c <- min(data_classes$C) - 1
  max_c <- max(data_classes$C) + 1

  min_rt <- max(0, floor(min(data_classes$RTsec)))
  max_rt <- ceiling(max(data_classes$RTsec))


  removed_data <- data_classes[data_classes$Filter.1 %in% removed_levels, , drop = FALSE]


  C <- RTsec <- DB <- LipidSubclass <- c()

  graph1 <- ggplot() +
    # puntos base
    geom_point(data = data_classes, aes(x = C, y = RTsec, color = factor(DB), shape = LipidSubclass),
               size = 2.5, stroke = 1.5) +

    # overlay eliminados
    geom_point(data = removed_data, aes(x = C, y = RTsec, color = factor(DB)),
               size = 5, shape = 22, alpha = 0.25, stroke = 1.2) +

    scale_color_manual(name = "Double \nbonds",
                       values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                  "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                  "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                  "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                  "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                  "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                  "#D9FF68",
                                  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                  "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                  "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                  "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                  "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                  "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                  "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                  "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                  "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                  "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                  "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                  "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                  "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                  "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                  "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                       breaks = seq(0:24)-1) +
    scale_shape_manual(values = shape_map) +
    labs(title = paste0("Combined filtering view: ", lipidsubclasses_combined),
         x = "Carbon number", y = "Retention time (s)") +
    theme_minimal() +
    coord_cartesian(xlim = c(min_c, max_c), ylim = c(min_rt, max_rt))


  graph2 <- ggplot() +
    # puntos base
    geom_point(data = data_classes[data_classes$checkfilter %in% kept_levels,],
               aes(x = C, y = RTsec, color = factor(DB), shape = LipidSubclass),
               size = 2.5, stroke = 1.5) +
    scale_color_manual(name = "Double \nbonds",
                       values = c("#212121", "#EE2C2C", "#FF9500", "#FFE100",
                                  "#90EE90", "#228B22", "#B0E2FF", "#00B2EE",
                                  "#3A5FCD", "#BF3EFF", "#68228B", "#795548",
                                  "#9E9E9E", "#FFB6C1", "#FFD39B", "#FFF68F",
                                  "#B8860B", "#8B864E", "#B87F9E", "#7FFF00",
                                  "#FF69B4", "#DEBF4B", "#51FFFC", "#EAC4FF",
                                  "#D9FF68",
                                  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                                  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                                  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
                                  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
                                  "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
                                  "#E7969C", "#CEDB9C", "#9C9EDE", "#DBDB8D", "#C49C94",
                                  "#F7B6D2", "#C7C7C7", "#AEC7E8", "#FFBB78", "#98DF8A",
                                  "#FF9896", "#C5B0D5", "#C49C94", "#FDD0A2", "#BDBDBD",
                                  "#5254A3", "#6B6ECF", "#9C9EDE", "#637939", "#B5CF6B",
                                  "#CEDB9C", "#8CA252", "#BD9E39", "#E7BA52", "#E7CB94",
                                  "#843C39", "#AD494A", "#D6616B", "#E7969C", "#7B4173",
                                  "#A55194", "#CE6DBD", "#DE9ED6", "#3182BD", "#6BAED6",
                                  "#9ECAE1", "#C6DBEF", "#E6550D", "#FD8D3C", "#FDAE6B",
                                  "#FDD0A2", "#31A354", "#74C476", "#A1D99B", "#C7E9C0",
                                  "#756BB1", "#9E9AC8", "#BCBDDC", "#DADAEB", "#636363",
                                  "#969696", "#BDBDBD", "#D9D9D9", "#393B79", "#5254A3",
                                  "#6B6ECF", "#9C9EDE", "#637939", "#8CA252", "#B5CF6B",
                                  "#CEDB9C", "#8C6D31", "#BD9E39", "#E7BA52", "#E7CB94"),
                       breaks = seq(0:24)-1) +
    scale_shape_manual(values = shape_map) +
    labs(x = "Carbon number", y = "Retention time (s)") +
    theme_minimal() +
    coord_cartesian(xlim = c(min_c, max_c), ylim = c(min_rt, max_rt))


  grDevices::pdf(NULL)
  grDevices::dev.control(displaylist = "enable")

  graph3 <- graph1 + graph2 +
    patchwork::plot_layout(ncol = 2, nrow = 1, guides = "collect")
  print(graph3)

  plots <- grDevices::recordPlot()

  invisible(grDevices::dev.off())

  return(plots)
}



# get_filtering_summary
#' Get a summary of RT-based filtering results to the console.
#'
#' This function create a formatted summary of the filtering step,
#' reporting the number of lipid identifications that were kept,
#' removed, recovered, or ignored (did not pass QC) for each lipid class.
#' The output is printed directly to the console, with lipid class names
#' highlighted (if the console supports ANSI formatting).
#'
#' @param rtdata rtdata object.
#'
#' @return information about filtering
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
get_filtering_summary <- function(rtdata) {

  lines <- c("Filtering Summary:\n")

  for (lipid in unique(rtdata$filtering$information$LipidSubclass)) {

    kept <- rtdata$filtering$information$keep[
      rtdata$filtering$information$LipidSubclass == lipid
    ] +
      rtdata$filtering$information$recover[
        rtdata$filtering$information$LipidSubclass == lipid
      ]

    removed <- rtdata$filtering$information$remove[
      rtdata$filtering$information$LipidSubclass == lipid
    ]

    ignored <- rtdata$filtering$information$ignored[
      rtdata$filtering$information$LipidSubclass == lipid
    ]

    lines <- c(
      lines,
      sprintf(
        "  - %-10s | %4d kept IDs | %4d removed IDs | %4d IDs did not pass QC",
        lipid, kept, removed, ignored
      )
    )
  }

  return(lines)
}


# printfilteringinfo
#' Print a summary of RT-based filtering results to the console.
#'
#' This function displays a formatted summary of the filtering step,
#' reporting the number of lipid identifications that were kept,
#' removed, recovered, or ignored (did not pass QC) for each lipid class.
#' The output is printed directly to the console, with lipid class names
#' highlighted (if the console supports ANSI formatting).
#'
#' @param rtdata rtdata object.
#'
#' @return information about filtering
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
printfilteringinfo <- function(rtdata) {
  cat(paste(get_filtering_summary(rtdata), collapse = "\n"))
  cat("\n")
}




# checkfiltervalue
#' Modify the filter tag generated during the filtering process.
#'
#' Modify the filter tag generated during the filtering process in a new variable
#' "Filter_id", ensuring that all lipids receive their corresponding tag at the
#' end. Initially, within the same group of lipids with the same class, number
#' of carbons, double bonds, and adducts, a tag is assigned to only one of them,
#' while the others are marked as ignored duplicates. However, at the end of the
#' process, it is necessary to consider all lipids and assign them their
#' corresponding tag.
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
checkfiltervalue <- function(data) {
  data$Filter_id <- NA

  LipidRTids <- unique(data$LipidRTid)

  for (LipidRTid in LipidRTids) {
    if ("Filter_final" %in% colnames(data)) {
      filtervalues <- unique(data[data$LipidRTid == LipidRTid &
                                    data$Filter_final != "Ignore duplicated", ])
    } else if ("checkfilter" %in% colnames(data)){
      filtervalues <- unique(data[data$LipidRTid == LipidRTid &
                                    data$checkfilter != "Ignore duplicated", ])
    }


    for (row in 1:nrow(filtervalues)) {
      if ("Filter_final" %in% colnames(data)) {
        value <- filtervalues$Filter_final[row]
      } else if ("checkfilter" %in% colnames(data)){
        value <- filtervalues$checkfilter[row]
      }
      lipidsubclass <- filtervalues$LipidSubclass[row]
      adduct <- filtervalues$adduct[row]

      condition <- data$LipidRTid == LipidRTid &
        data$LipidSubclass == lipidsubclass &
        data$adduct == adduct

      if (value %in% c("Kept", "Kept manually by the user", "Recovered with RT modelling")) {
        data$Filter_id[condition] <- "Kept"
      } else if (value %in% c("Not pass QC")) {
        data$Filter_id[condition] <- "Not pass QC"
      } else if (value == "Unsupported class") {
        data$Filter_id[condition] <- "Unsupported class"
      } else {
        # data$Filter_id[condition] <- "remove"
        data$Filter_id[condition] <- "Removed"
      }
    }
  }

  return(data)
}

# groupidfilter
#' Update annotated peaklist with the filtered data from external data
#'
#' Update the annotated peaklist with the filtered data, using the "Filter_id" tag.
#' Add the two new filtering columns: "filter_label" and "lipidname_filter"
#'
#' @param rtdata rtdata
#'
#' @return featuresFilter
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
groupidfilter <- function(rtdata) {

  features <- rtdata$features
  source   <- rtdata$metadata$source
  filteringdata <- rtdata$filtering$data

  filteringdata$confidenceLevel <- factor(
    filteringdata$confidenceLevel,
    levels = c("FA position", "FA", "Subclass", "MS-only", "MSMS", "Predicted",
               "Other", "Unsupported class", "Not pass QC", ""),
    ordered = TRUE)

  ## -----------------------
  ## Source-specific config
  ## -----------------------
  sep.lipid <- "|"

  if (source == "lipidms") {
    position_cols <- c("mz", "minmz", "maxmz", "RT", "minRT", "maxRT",
                       "iniRT", "endRT", "RTminutes", "npeaks", "group",
                       "isotope", "isoGroup", "LipidMSid", "Adduct",
                       "confidenceLevel", "Score", "ScoreInt",
                       "nsamples", "LipidRTid")

  } else if (source == "msdial") {
    position_cols <- c("Alignment ID", "Average Rt(min)", "Average Mz",
                       "Metabolite name", "Adduct type", "LipidRTid")

  } else {
    if ("RTsec" %in% colnames(features)) {
      position_cols <- c("mz", "RTsec", "lipidname",
                         "adduct", "confidenceLevel", "LipidRTid")
    } else {
      position_cols <- c("mz", "RTmin", "lipidname",
                         "adduct", "confidenceLevel", "LipidRTid")
    }
  }

  ## -----------------------
  ## Main loop by LipidRTid
  ## -----------------------
  filter_data <- list()
  LipidRTids <- unique(filteringdata$LipidRTid[order(filteringdata$LipidRTid)])

  # if ("Filter_final" %in% colnames(filteringdata)) { # ya se hace en algún lado, Filter_id cambia con Filter_final
  #   # filteringdata$checkfilter <- filteringdata$Filter_final
  #   if ("Filter_id" %in% colnames(filteringdata)) {
  #     filteringdata$Filter_id[
  #       filteringdata$Filter_id != "Kept" &
  #         filteringdata$Filter_final %in% c(
  #           "Kept",
  #           "Kept manually by the user",
  #           "Recovered with RT modelling"
  #         )
  #     ] <- "Kept"
  #   }
  # }

  for (LipidRTid in LipidRTids) {

    subset_data <- filteringdata[filteringdata$LipidRTid == LipidRTid, ]
    subset_data <- subset_data[order(subset_data$roworder), ]

    ## Preallocate
    n <- nrow(subset_data)

    id <- lipidclass <- lipidsubclass <- c_val <- db <- adducts <-
      confidencelevel <- score <- label <- comments <- groupid <- lipidnameexp <-
      character(n)

    if (!("Filter_id" %in% colnames(subset_data))) {
      kept        <- subset_data$checkfilter %in% c("Kept",
                                                    "Kept manually by the user",
                                                    "Recovered with RT modelling")
      notqc       <- subset_data$checkfilter == "Not pass QC"
      unsupported <- subset_data$checkfilter == "Unsupported class"
      other       <- !(kept | notqc | unsupported)
    } else {
      kept        <- subset_data$Filter_id == "Kept"
      notqc       <- subset_data$Filter_id == "Not pass QC"
      unsupported <- subset_data$Filter_id == "Unsupported class"
      other       <- !(kept | notqc | unsupported)
    }


    ## Kept
    id[kept]              <- subset_data$lipidname[kept]
    lipidclass[kept]      <- subset_data$LipidClass[kept]
    lipidsubclass[kept]   <- subset_data$LipidSubclass[kept]
    c_val[kept]           <- subset_data$C[kept]
    db[kept]              <- subset_data$DB[kept]
    adducts[kept]         <- subset_data$adduct[kept]
    confidencelevel[kept] <- as.character(subset_data$confidenceLevel[kept])
    score[kept]           <- subset_data$Score_rtpred[kept]
    comments[kept]        <- subset_data$comments[kept]
    groupid[kept]         <- subset_data$GroupID[kept]
    lipidnameexp[kept]     <- subset_data$lipidname_expanded[kept]

    ## Not pass QC
    confidencelevel[notqc] <- as.character(subset_data$confidenceLevel[notqc])
    comments[notqc]        <- subset_data$comments[notqc]
    groupid[notqc]         <- subset_data$GroupID[notqc]
    lipidnameexp[notqc]     <- subset_data$lipidname_expanded[notqc]
    label[notqc]           <- "Not pass QC"
    comments[notqc]        <- subset_data$comments[notqc]

    ## Unsupported class
    id[unsupported]              <- subset_data$lipidname[unsupported]
    adducts[unsupported]         <- subset_data$adduct[unsupported]
    confidencelevel[unsupported] <- as.character(subset_data$confidenceLevel[unsupported])
    comments[unsupported]        <- subset_data$comments[unsupported]
    groupid[unsupported]         <- subset_data$GroupID[unsupported]
    lipidnameexp[unsupported]     <- subset_data$lipidname_expanded[unsupported]
    label[unsupported]           <- "Unsupported class"
    comments[unsupported]        <- subset_data$comments[unsupported]

    ## Other
    label[other]    <- subset_data$checkfilter[other]
    comments[other] <- subset_data$comments[other]

    ## -----------------------
    ## Collapse logic
    ## -----------------------
    if (all(id == "")) {
      id_str <- lipidclass_str <- lipidsubclass_str <- c_str <- db_str <-
        adducts_str <- confidencelevel_str <- score_str <-
        groupid_str <- lipidnameexp_str <- comments_str <- ""
    } else {

      id_str <- paste(id[id != ""], collapse = sep.lipid)

      same_core <-
        length(unique(lipidclass[lipidclass != ""])) == 1 &&
        length(unique(lipidsubclass[lipidsubclass != ""])) == 1 &&
        length(unique(c_val[c_val != ""])) == 1 &&
        length(unique(db[db != ""])) == 1  &&
        length(unique(adducts[adducts != ""])) == 1

      lipidclass_str     <- collapse_conditional(lipidclass, sep.lipid, !same_core)
      lipidsubclass_str  <- collapse_conditional(lipidsubclass, sep.lipid, !same_core)
      c_str              <- collapse_conditional(c_val, sep.lipid, !same_core)
      db_str             <- collapse_conditional(db, sep.lipid, !same_core)
      adducts_str             <- collapse_conditional(adducts, sep.lipid, !same_core)
      confidencelevel_str             <- collapse_conditional(confidencelevel, sep.lipid, !same_core)
      score_str             <- collapse_conditional(score, sep.lipid, !same_core)
      groupid_str             <- collapse_conditional(groupid, sep.lipid, !same_core)
      lipidnameexp_str   <- collapse_conditional(lipidnameexp, sep.lipid, !same_core)
      comments_str <- collapse_conditional(comments, sep.lipid, !same_core)
    }

    label_str    <- if (all(label == "")) "" else paste(label, collapse = sep.lipid)
    # comments_str <- if (all(comments == "")) "" else paste(comments, collapse = sep.lipid)

    filter_data[[length(filter_data) + 1]] <-
      data.frame(
        GroupID = groupid_str,
        LipidRTid = LipidRTid,
        filter_label = label_str,
        comments_filter = comments_str,
        lipidname_filter = id_str,
        expandedlipidname_filter = lipidnameexp_str,
        LipidClass_filter = lipidclass_str,
        LipidSubclass_filter = lipidsubclass_str,
        C_filter = c_str,
        DB_filter = db_str,
        adduct_filter = adducts_str,
        confidenceLevel_filter = confidencelevel_str,
        Score_filter = score_str,
        check.names = FALSE,
        stringsAsFactors = FALSE

      )
  }

  filter_data <- do.call(rbind, filter_data)

  ## -----------------------
  ## Merge back to features
  ## -----------------------
  filter_cols <- c(
    "GroupID","filter_label","comments_filter","lipidname_filter", "expandedlipidname_filter",
    "LipidClass_filter","LipidSubclass_filter","C_filter","DB_filter",
    "adduct_filter","confidenceLevel_filter","Score_filter"
  )

  features <- features[, !names(features) %in% filter_cols, drop = FALSE]
  features <- merge(features, filter_data, by = "LipidRTid", all.x = TRUE)

  col_order <- c(position_cols, filter_cols,
                 setdiff(colnames(features), c(position_cols, filter_cols)))

  features <- features[, col_order[col_order %in% colnames(features)], drop = FALSE]

  for (col in filter_cols) {
    features[[col]][is.na(features[[col]]) | features[[col]] == "NA"] <- ""
  }

  rtdata$features <- features


  # summarizedFeatures
  cols <- c("LipidRTid", "GroupID", "mz", "RTsec", "RTmin", "lipidname_filter", "expandedlipidname_filter",
            "LipidClass_filter", "LipidSubclass_filter", "C_filter", "DB_filter", "adduct_final",
            "confidenceLevel_filter", "Score_filter", "comments_filter")

  # Encontrar la última columna cuyo nombre termina en "_final"
  final_cols <- which(grepl("_filter$", names(features)))

  # Columnas que aparecen después de la última "_final"
  if (length(final_cols) > 0 && final_cols[length(final_cols)] < ncol(features)) {

    new_cols <- names(features)[
      (final_cols[length(final_cols)] + 1):ncol(features)
    ]

    # Añadir únicamente las que todavía no están en cols
    cols <- c(
      cols,
      new_cols[!new_cols %in% cols]
    )
  }

  summ <- rtdata$features

  if ("RTsec" %in% colnames(summ)) {
    summ$RTsec <- summ$RTsec
  } else if ("RT" %in% colnames(summ)) {
    summ$RTsec <- summ$RT
  } else if ("RTmin" %in% colnames(summ)) {
    summ$RTsec <- as.numeric(summ$RTmin) * 60
    summ$RTmin <- summ$RTmin
  } else if ("Average Rt(min)" %in% colnames(summ)) {
    summ$RTsec <- as.numeric(summ$`Average Rt(min)`) * 60
    summ$RTmin <- summ$`Average Rt(min)`
  } else {
    stop("Neither RTsec, RTmin nor RT column found in features")
  }

  # RTmin
  if (!"RTmin" %in% colnames(summ)) {
    summ$RTmin <- as.numeric(summ$RTsec) / 60
  }

  # RTsec
  if (!"RTsec" %in% colnames(summ) &
      "RTmin" %in% colnames(summ)) {
    summ$RTsec <- as.numeric(summ$RTmin) * 60
  }

  # Seleccionar columnas disponibles respetando el orden
  # final_cols <- cols[cols %in% colnames(summ)]
  final_cols <- cols[cols %in% colnames(summ)]

  rtdata$summarizedFeatures <- summ[, final_cols, drop = FALSE]

  filter_idx <- grepl("_filter$", colnames(rtdata$summarizedFeatures))
  colnames(rtdata$summarizedFeatures)[filter_idx] <-
    sub("_filter$", "_final", colnames(rtdata$summarizedFeatures)[filter_idx])

  return(rtdata)
}

