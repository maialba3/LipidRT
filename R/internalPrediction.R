# rtpredictannotationsdata
#' Retention time modeling to propose new identifications for the features.
#'
#' Retention time modeling to propose new identifications for the features.
#' New identifications are suggested for features that have not been previously
#' annotated, using potential candidates selected based on the m/z value.
#' Once potential candidates for the features are obtained, retention time is
#' modeled for each lipid class using previous annotations. Subsequently,
#' retention time is predicted for each candidate.
#'
#' @param rtdata rtdata
#' @param distribution character. Data distribution for RT modelling. It can be "gamma" or "normal".
#' @param classes classes
#' @param distC_predict integer. Minimum difference in number of carbons from which lipids
#' are considered distant from the rest of lipids from the same class.
#' @param minSpecies_predict numeric. Minimum number of carbon units from any valid model data or prediction beyond which predictions will be discarded. This value must be at least 3, with a default of 5.
#' @param rtscorethr numeric. Minimum Score_rtpred value from which a prediction will be considered valid.
#' @param rt_thr numeric. Threshold for the difference (in seconds) between the
#' observed retention time (RT) and the predicted retention time (RTsec_rtpred).
#' @param rttol_id Numeric. Maximum allowed retention-time difference (in
#' seconds) to consider two observations part of the same group.
#' @param rttol_name Numeric. Maximum allowed retention-time difference (in
#' seconds) to consider two lipids as temporally close.
#' @param predict_plot logical. If TRUE, predictions are plotted. If FALSE, predictions are not plotted.
#'
#' @return rtdata
#'
#' @examples
#' \dontrun{
#' rtdata <- rtpredictannotationsdata(rtdata, distribution = "gamma", dmz = 5,
#'                                    ppm = TRUE,  rtscorethr = 0.5, minC_predict = 5)
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
rtpredictannotationsdata <- function(rtdata, distribution, classes = NULL,
                                     distC_predict = 5, minSpecies_predict = 5,
                                     rtscorethr = 0.5,
                                     rt_thr = 30, rttol_id = 5, rttol_name = 1,
                                     predict_plot = TRUE) {

    # save rtscorethr parameter in rtdata
    rtdata$metadata$prediction$rtscorethr <- as.numeric(rtscorethr)

    # Annotated peaklist
    data <- rtdata$filtering$data

    # All valid classes in the dataset
    classes_validas <- unique(na.omit(data$LipidSubclass[data$LipidSubclass != ""]))
    if (is.null(classes)) {
      classes <- classes_validas
    } else {
      classes <- classes[classes %in% classes_validas]
    }

    # Add "comments" column
    if (is.null(rtdata$prediction$candidates$comments)) {
      if (nrow(rtdata$prediction$candidates) > 0) {
        rtdata$prediction$candidates$comments <- ""
      } else {
        rtdata$prediction$candidates$comments <- character(0)
      }
    }


    ############################################################################
    # RT modelling and predictions
    rtdata <- rtmodelling(rtdata = rtdata,
                          distribution = distribution,
                          classes = classes,
                          typeanalysis = "prediction")

    rtdata <- predictannotations(rtdata = rtdata,
                                 minC_predict = minSpecies_predict,
                                 rtscorethr = rtscorethr,
                                 rt_thr = rt_thr,
                                 typeanalysis = "prediction")


    # Group unique IDs
    rtdata <- groupuniqueids(rtdata,
                             typeanalysis = "prediction",
                             rttol_id = rttol_id)

    # Expand candidates lipid name
    rtdata <- lipidnameexpanded(rtdata,
                                typeanalysis = "prediction",
                                rttol_name = rttol_name)

    # filter lipid predictions based on the level of adduct support
    rtdata <- filterpredictionsbyadducts(rtdata)

    # Recover lipid names of filtered lipids
    rtdata <- recoverlipidnamespredicted(rtdata)

    #==========================================================================#
    # Filter prediction results
    #==========================================================================#
    rtdata <- filterpredictionsdata(rtdata,
                                    distC_predict = distC_predict,
                                    minSpecies_predict = minSpecies_predict,
                                    plot = TRUE)

    #==========================================================================#
    # Plot prediction results
    #==========================================================================#
    if (predict_plot == TRUE) {
      rtdata <- predictgraph(rtdata = rtdata,
                             classes = classes,
                             eval_predict_plot = TRUE)
    }

  return(rtdata)
}


# targetpredictions
#' Retention time modeling to propose new identifications for the features
#'
#' Retention time modeling to propose new identifications for the features.
#' New identifications are suggested for features that have not been previously
#' annotated, using potential candidates selected based on the m/z value.
#' Once potential candidates for the features are obtained, retention time is
#' modeled for each lipid class using previous annotations. Subsequently,
#' retention time is predicted for each candidate.
#'
#' @param predictions rtdata$prediction$predictions
#'
#' @details Once
#'
#' @return predictions
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
targetpredictions <- function(predictions) {

  if (is.null(predictions) || !is.data.frame(predictions) || nrow(predictions) == 0) {
    return(data.frame(
      lipidname = character(0),
      RTsec_rtpred = numeric(0),
      CI_lower = numeric(0),
      CI_upper = numeric(0)
    ))
  }

  predictions$lipidname <- ifelse(
    grepl("_d$", predictions$LipidSubclass),
    paste0(gsub("_d$", "", predictions$LipidSubclass), "(d", predictions$C, ":", predictions$DB, ")"),
    ifelse(
      grepl("_t$", predictions$LipidSubclass),
      paste0(gsub("_t$", "", predictions$LipidSubclass), "(t", predictions$C, ":", predictions$DB, ")"),
      paste0(predictions$LipidSubclass, "(", predictions$C, ":", predictions$DB, ")")
    )
  )

  predictions$lipidname <- ifelse(
    grepl("-O$", predictions$LipidSubclass),
    paste0(gsub("-O$", "", predictions$LipidSubclass), "(O-", predictions$C, ":", predictions$DB, ")"),
    ifelse(
      grepl("-P$", predictions$LipidSubclass),
      paste0(gsub("-P$", "", predictions$LipidSubclass), "(P-", predictions$C, ":", predictions$DB, ")"),
      paste0(predictions$LipidSubclass, "(", predictions$C, ":", predictions$DB, ")")
    )
  )

  predictions$RTsec_rtpred <- round(as.numeric(predictions$RTsec_rtpred), 3)
  predictions$CI_lower     <- round(as.numeric(predictions$CI_lower), 3)
  predictions$CI_upper     <- round(as.numeric(predictions$CI_upper), 3)

  predictions <- predictions[,
                             c("lipidname", "RTsec_rtpred", "CI_lower", "CI_upper"),
                             drop = FALSE
  ]

  return(predictions)
}


# getmetrics
#' Model metrics
#'
#' Model metrics
#'
#' @param model model
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
# getmetrics <- function(model, data) {
#   pred <- predict(model, type = "response")
#   actual <- data$RTsec
#   rmse <- sqrt(mean((actual - pred)^2))
#   pseudoR2 <- 1 - (model$deviance / model$null.deviance)
#   list(AIC = AIC(model),
#        logLik = as.numeric(logLik(model)),
#        RMSE = rmse,
#        pseudoR2 = pseudoR2)
# }
getmetrics <- function(model, data) {
  # Use model data to ensure alignment
  mf <- model.frame(model)
  pred <- predict(model, type = "response")
  actual <- mf$RTsec

  rmse <- sqrt(mean((actual - pred)^2))
  pseudoR2 <- 1 - (model$deviance / model$null.deviance)
  list(AIC = AIC(model),
       logLik = as.numeric(logLik(model)),
       RMSE = rmse,
       pseudoR2 = pseudoR2)
}


# rtmodelling
#' Model the retention time for each lipid class
#'
#' Model the retention time for each lipid class
#'
#' @param rtdata annotated peaklist, with each annotation separated into a row (filtering$data)
#' @param distribution character. Data distribution for modeling. It can be 'normal' or 'gamma'.
#' @param classes A vector containing the lipid classes for which retention time
#' is to be predicted. If classes = NULL, predictions will be made for all classes
#' where there are candidates and data to create a model.
#' @param typeanalysis Character. Specifies whether the function is being used for
#' filtering or RT predictions. Possible values are "filtering" or "prediction".
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
rtmodelling <- function(rtdata, distribution = "gamma", classes = NULL,
                        typeanalysis) {

  data <- rtdata$filtering$data

  if ("checkfilter" %in% colnames(data)) {
    keep_values <- c("Kept","Kept manually by the user","Recovered with RT modelling")

    data <- data[!is.na(data$lipidname) & data$lipidname != "" &
                   !is.na(data$LipidSubclass) &
                   data$checkfilter %in% keep_values, ]

    if (!is.null(classes)) {
      data <- data[data$LipidSubclass %in% classes, ]
    }

  } else {
    data <- data[!is.na(data$lipidname) & data$lipidname != "" &
                   !is.na(data$LipidSubclass) &
                   data$Filter == "Kept", ]
  }

  data <- grouplipidclasses(data, map_type = "class")
  lipidsubclasses <- unique(na.omit(data$LipidSubclass))

  models <- list()

  for (lipidsubclass in lipidsubclasses) {

    lipidclass <- unique(data$LipidClass[data$LipidSubclass == lipidsubclass])

    data_class <- data[data$LipidClass == lipidclass, ]
    data_subclass <- data[data$LipidSubclass == lipidsubclass, ]

    subclasses  <- unique(data_class$LipidSubclass)
    category    <- unique(data_subclass$groupClass)

    model_out <- NULL

    # ===========================================================
    # 1. Class/Subclass level
    # ===========================================================
    if (has_enough_data(data = data_class, mode = "class")) {

      if (length(subclasses) > 1 & category != "PX") {
        model_out <- create_model(data = data_class,
                                  distribution = distribution,
                                  include_subclass = TRUE)
      } else {
        model_out <- create_model(data = data_subclass,
                                  distribution = distribution,
                                  include_subclass = FALSE)
      }

      if (!is.null(model_out) && length(subclasses) > 1) {

        # subclases _d y _t
        subclass_d <- subclasses[grepl("_d$", subclasses)]
        subclass_t <- subclasses[grepl("_t$", subclasses)]

        if (length(subclass_d) == 1 && length(subclass_t) == 1) {

          coef_names <- names(coef(model_out))

          coef_d_name <- paste0("LipidSubclass", subclass_d)
          coef_t_name <- paste0("LipidSubclass", subclass_t)

          coef_d <- if (coef_d_name %in% coef_names) coef(model_out)[coef_d_name] else NA
          coef_t <- if (coef_t_name %in% coef_names) coef(model_out)[coef_t_name] else NA

          wrong_order <-
            (!is.na(coef_d) && coef_d <= 0) ||
            (!is.na(coef_t) && coef_t >= 0)

          if (wrong_order) {
            warning_msg <- "Incorrect RT elution order between subclasses"

            if (typeanalysis == "prediction") {
              rtdata$prediction$candidates$comments <-
                unique(c(rtdata$prediction$candidates$comments, warning_msg))
            }

            if (typeanalysis == "filtering") {
              rtdata$filtering$data$comments <-
                unique(c(rtdata$filtering$data$comments, warning_msg))
            }
          }
        }
      }

    } else {

      # ===========================================================
      # 2. Category level (work at subclass level directly)
      # ===========================================================
      model_map <- grouplipidclasses(data = NULL, return_map = TRUE, map_type = "model")

      model_classes <- model_map[[category]]

      if (!is.null(model_classes)) {

        # Get all data belonging to the category classes
        data_category <- data[data$LipidSubclass %in% model_classes, ]

        # Identify subclasses within this category
        category_subs <- unique(data_category$LipidSubclass)

        # Keep only subclasses with enough data
        valid_subclasses <- Filter(function(s) {
          df_sub <- data_category[data_category$LipidSubclass == s, ]
          has_enough_data(data = df_sub, mode = "category")
        }, category_subs)

        if (length(valid_subclasses) > 0) {

          valid_subclasses <- unique(c(lipidsubclass, valid_subclasses))

          # Build dataset only with valid subclasses
          data_group <- data_category[data_category$LipidSubclass %in% valid_subclasses, ]

          # Decide if subclass structure is meaningful in category model
          include_sub <- length(unique(data_group$LipidSubclass)) > 1

          model_out <- create_model(data = data_group,
                                    distribution = distribution,
                                    include_subclass = include_sub)
        }
      }
    }

    # ===========================================================
    # Save model
    # ===========================================================
    if (!is.null(model_out)) {
      models[[lipidsubclass]] <- model_out
    }
  }

  # ===============================================================
  # 5) GUARDAR RESULTADOS
  # ===============================================================
  if (typeanalysis == "filtering") rtdata$filtering$models <- models
  if (typeanalysis == "prediction") rtdata$prediction$models <- models

  return(rtdata)
}


# has_enough_data
#' Model the retention time for each lipid class
#'
#' Model the retention time for each lipid class
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#' @param mode character. Data distribution for modeling. It can be 'normal' or 'gamma'.
#'
#' @return condition. logical
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
has_enough_data <- function(data, mode = c("class", "category")) {
  mode <- match.arg(mode)

  # --------------------------------------------------
  # CLASS MODE
  # --------------------------------------------------
  if (mode == "class") {

    subclasses <- unique(data$LipidSubclass)

    for (sc in subclasses) {

      subdata <- data[data$LipidSubclass == sc, , drop = FALSE]

      n_db <- sapply(unique(subdata$DB), function(db) {
        length(unique(subdata$C[subdata$DB == db]))
      })

      cond0 <- (sum(n_db >= 3) >= 1 && sum(n_db) >= 5)
      cond1 <- (sum(n_db >= 2) >= 3 && sum(n_db) > 5)

      if (cond0 || cond1) {
        return(TRUE)
      }
    }

    return(FALSE)
  }

  # --------------------------------------------------
  # CATEGORY MODE (unchanged)
  # --------------------------------------------------
  if (mode == "category") {

    n_db <- sapply(unique(data$DB), function(db) {
      length(unique(data$C[data$DB == db]))
    })

    cond0 <- (sum(n_db >= 3) >= 1 && sum(n_db) >= 5)
    cond1 <- (length(which(n_db >= 2)) >= 2 && sum(n_db) >= 5)
    cond2 <- (sum(n_db >= 2) >= 3)

    return(cond0 || cond1 || cond2)
  }
}

# create_model
#' Model the retention time for each lipid class
#'
#' Model the retention time for each lipid class
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#' @param distribution distribut
#' @param include_subclass character. Data distribution for modeling. It can be 'normal' or 'gamma'.
#'
#' @return model
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
create_model <- function(data, distribution, include_subclass = FALSE) {
  form <- RTsec ~ I(C^2) + C + I(1/(1+exp(0.1*C))) + DB
  if (include_subclass) form <- stats::update(form, . ~ . + LipidSubclass)

  fam <- if (tolower(distribution) == "normal")
    gaussian(link = "identity") else Gamma(link = "log")

  model <- glm(form, family = fam, data = data)

  return(list(
    model = model,
    # metrics = getmetrics(model, data)
    metrics = getmetrics(model)
  ))
}


# predictannotations_glm
#' predict retention time for each of the potential candidates
#'
#' predict retention time for each of the potential candidates
#'
#' @param model model
#' @param newdata dataframe. data to predict
#' @param level level = 1 - alpha. ICs lower, upper
#'
#' @return predictions
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
predictannotations_glm <- function(model, newdata, level = 0.95) {

  pred_link <- predict(model, newdata, type = "link", se.fit = TRUE)

  alpha <- 1 - level
  z <- qnorm(1 - alpha/2)

  fit_link   <- pred_link$fit
  se_link    <- pred_link$se.fit
  lower_link <- fit_link - z * se_link
  upper_link <- fit_link + z * se_link

  inv_link <- stats::family(model)$linkinv

  fit_resp   <- inv_link(fit_link)
  lower_resp <- inv_link(lower_link)
  upper_resp <- inv_link(upper_link)

  data.frame(
    fit = as.numeric(fit_resp),
    lwr = as.numeric(lower_resp),
    upr = as.numeric(upper_resp)
  )
}


# predictannotations
#' predict retention time for each of the potential candidates
#'
#' predict retention time for each of the potential candidates
#'
#' @param rtdata rtdata
#' @param rtscorethr numeric. Minimum Score_rtpred value from which a prediction will be considered valid.
#' @param rt_thr numeric. Threshold for the difference (in seconds) between the
#' observed retention time (RT) and the predicted retention time (RTsec_rtpred).
#' @param minC_predict numeric. Minimum number of carbon units from any valid model
#' data or prediction beyond which predictions will be discarded. This value must
#' be at least 3, with a default of 5.
#' @param typeanalysis Character. Specifies whether the function is being used for
#' filtering or RT predictions. Possible values are "filtering" or "prediction".
#'
#' @return predictions
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
predictannotations <- function(rtdata, rtscorethr, rt_thr, minC_predict, typeanalysis) {

  typeanalysis <- match.arg(typeanalysis, choices = c("filtering", "prediction"))

  data <- rtdata$filtering$data

  if (typeanalysis == "filtering") {
    models = rtdata$filtering$models
  } else if (typeanalysis == "prediction") {
    candidates = rtdata$prediction$candidates
    models = rtdata$prediction$models
  }

  predictions <- data.frame()
  if (length(models) > 0) {
    for (lipidsubclass in names(models)) {

      # Data on which to predict whether filtering or predicting on candidates
      if (typeanalysis == "filtering") { # Filtering ...........................

        data_subclass <- data[data$LipidSubclass == lipidsubclass &
                             data$Filter == "Kept", , drop = FALSE ]

        #  Data for prediction using a multi-class model
        if (length(unique(models[[lipidsubclass]][["model"]]$data$LipidSubclass)) > 1 |
            length(unique(models[[lipidsubclass]][["model"]]$model$LipidSubclass)) > 1) {
          newdat <-  expand.grid(LipidSubclass = unique(data_subclass$LipidSubclass),
                                 C = seq(min(data_subclass$C) - 1,
                                         max(data_subclass$C) + 1, by = 1),
                                 DB = seq(min(data_subclass$DB),
                                          max(data_subclass$DB), by = 1))
        } else {
          # Data for prediction using a single-class model
          newdat <- expand.grid(LipidSubclass = unique(data_subclass$LipidSubclass),
                                C = seq(min(data_subclass$C) - 1,
                                        max(data_subclass$C) + 1, by = 1),
                                DB = seq(min(data_subclass$DB),
                                         max(data_subclass$DB), by = 1))
        }

      } else if (typeanalysis == "prediction") { # Prediction ..................

        data_subclass <- data[data$LipidSubclass == lipidsubclass &
                             data$checkfilter %in% c("Kept",
                                                     "Kept manually by the user",
                                                     "Recovered with RT modelling"), , drop = FALSE]


        candidates_class <- candidates[candidates$LipidSubclass == lipidsubclass, , drop = FALSE]


        if (nrow(candidates_class) > 0) {
          #  Data for prediction using a multi-class model
          if (length(unique(models[[lipidsubclass]][["model"]]$data$LipidSubclass)) > 1 |
              length(unique(models[[lipidsubclass]][["model"]]$model$LipidSubclass)) > 1) {

            newdat <- expand.grid(LipidSubclass = unique(c(data_subclass$LipidSubclass, candidates_class$LipidSubclass)),
                                  C = seq(min(union(data_subclass$C, candidates_class$C)) - 1,
                                          max(union(data_subclass$C, candidates_class$C)) + 1, by = 1),
                                  DB = seq(min(union(data_subclass$DB, candidates_class$DB)),
                                           max(union(data_subclass$DB, candidates_class$DB)), by = 1))
          } else {
            # Data for prediction using a single-class model
            newdat <- expand.grid(LipidSubclass = unique(c(data_subclass$LipidSubclass, candidates_class$LipidSubclass)),
                                  C = seq(min(union(data_subclass$C, candidates_class$C)) - 1,
                                          max(union(data_subclass$C, candidates_class$C)) + 1, by = 1),
                                  DB = seq(min(data_subclass$DB),
                                           max(data_subclass$DB), by = 1))
          }


        } else {
          #  Data for prediction using a multi-class model
          if (length(unique(models[[lipidsubclass]][["model"]]$data$LipidSubclass)) > 1 |
              length(unique(models[[lipidsubclass]][["model"]]$model$LipidSubclass)) > 1) {
            newdat <- expand.grid(LipidSubclass = unique(data_subclass$LipidSubclass),
                                  C = seq(min(data_subclass$C) - 1,
                                          max(data_subclass$C) + 1, by = 1),
                                  DB = seq(min(data_subclass$DB),
                                           max(data_subclass$DB), by = 1))
          } else {
            # Data for prediction using a single-class model
            newdat <- expand.grid(LipidSubclass = unique(data_subclass$LipidSubclass),
                                  C = seq(min(data_subclass$C) - 1,
                                          max(data_subclass$C) + 1, by = 1),
                                  DB = seq(min(data_subclass$DB),
                                           max(data_subclass$DB), by = 1))
          }

        }

      } # ......................................................................

      # Predictions
      pred_ci <- predictannotations_glm(model = models[[lipidsubclass]][["model"]],
                                        newdata = newdat,
                                        level = 0.95)

      prediction <- data.frame(newdat,
                               "RTsec_rtpred" = round(pred_ci$fit, 3),
                               "CI_lower" = round(pred_ci$lwr, 3),
                               "CI_upper" = round(pred_ci$upr, 3),
                               stringsAsFactors = FALSE)


      prediction <- prediction[, c("LipidSubclass", "C", "DB",
                                   "RTsec_rtpred", "CI_lower", "CI_upper")]


      prediction <- prediction[order(prediction$C, prediction$DB), , drop = FALSE]

      # Calculate the difference in RT between the different DB with the same C for the score
      prediction$diffrt <- NA
      for (c in unique(prediction$C)) { # For each possible C value of the predictions...
        subset_data <- prediction[prediction$C == c, ]
        distances <- numeric(length(unique(subset_data$DB)) - 1)
        for (db in unique(subset_data$DB)) { # The possible DB values for that C value
          for (i in 1:length(distances)) { # Calculate the distance between the predictions of that C value and the different DB values
            distances[i] <- (abs(as.numeric(subset_data$RTsec_rtpred[i+1]) -
                                   as.numeric(subset_data$RTsec_rtpred[i])))/
              (as.numeric(subset_data$DB[i+1]) -
                 as.numeric(subset_data$DB[i]))
          }
        }
        # The average of the distances between the different DB values for a given C value will be used to calculate the Score_p
        prediction$diffrt[prediction$C == c] <- mean(distances)
      }

      predictions <- rbind(predictions, prediction)
    }

    # If it is a joint model and there are no data to calculate the diffRT for that class,
    # take the average from the different classes used in the model with the same C.
    if (any(is.na(predictions$diffrt))) {
      lipidsubclasses <- unique(predictions$LipidSubclass[is.na(predictions$diffrt)])

      for (lipidsubclass in lipidsubclasses) {
        subclasses_model <- unique(models[[lipidsubclass]][["model"]]$model$LipidSubclass)

        if (length(subclasses_model) > 1) {
          # For each C value with NA in that class
          C_na <- unique(predictions$C[predictions$LipidSubclass == lipidsubclass & is.na(predictions$diffrt)])

          for (c_val in C_na) {
            # Filter data with those model classes and that C
            mask <- predictions$LipidSubclass %in% subclasses_model & predictions$C == c_val & !is.na(predictions$diffrt)
            diffrt_mean <- mean(predictions$diffrt[mask])

            # Assign the mean of diffrt to the corresponding NA positions
            predictions$diffrt[predictions$LipidSubclass == lipidsubclass & predictions$C == c_val & is.na(predictions$diffrt)] <- diffrt_mean
          }
        }
      }
    }

    if (typeanalysis == "filtering") {
      rtdata$filtering$predictions <- predictions

    } else if (typeanalysis == "prediction") {
      rtdata$prediction$predictions <- predictions
    }

  } else { # If no models are available
    if (typeanalysis == "filtering") {

      # Predictions results null
      predictions <- data.frame()
      for (lipidsubclass in unique(data$LipidSubclass[data$LipidSubclass != ""])) {
        data_subclass <- data[data$LipidSubclass == lipidsubclass, ]
        newdat <- unique(data_subclass[, c("LipidSubclass", "C", "DB")])

        # No models are available for prediction, so all values are set to NA
        newdat$RTsec_rtpred <- NA
        newdat$CI_lower <- NA
        newdat$CI_upper <- NA
        newdat$RTsec_rtpred_prev <- NA
        newdat$diffrt <- NA
        newdat$slope <- NA

        predictions <- rbind(predictions, newdat)
      }

      rtdata$filtering$predictions <- predictions

    } else if (typeanalysis == "prediction") {

      # predictions results null
      predictions <- data.frame(LipidSubclass = character(),
                                C = integer(),
                                DB = integer(),
                                RTsec_rtpred = double(),
                                CI_lower = double(),
                                CI_upper = double(),
                                diffrt = double(),
                                RTsec_rtpred_prev = double(),
                                slope = double(),
                                stringsAsFactors = FALSE)

      rtdata$prediction$predictions <- predictions
    }
  }


  # Calculate predictions score
  if (typeanalysis == "prediction") { # Prediction .............................
    # Add the RT predictions to the proposed candidates
    if ("RTsec_rtpred" %in% colnames(candidates)) {
      candidates <- candidates[ , !(names(candidates) %in% c("RTsec_rtpred", "CI_lower", "CI_upper",
                                                             "RTsec_rtpred_prev", "diffrt", "slope", "range",
                                                             "Score_p", "Score_rtpred", "checkcandidate")) ]
    }

    # Apply RT bounds if local maxima and minima are available
    rtdata <- restrictrtpred(rtdata = rtdata,
                             typeanalysis = "prediction")

    candidates <- merge(candidates, rtdata$prediction$predictions,
                        by = c("LipidSubclass", "C", "DB"), all.x = TRUE)

    # Calculate the prediction score (Score_p)
    rtdata$prediction$candidates <- predictionscore(data = candidates)

    # Look for distant predictions (range = "out of range")  -- lo quitamos
    # rtdata <- distantprediction(rtdata = rtdata,
    #                             rtscorethr = rtscorethr,
    #                             minC_predict = minC_predict)

    # Generate Score_rtpred
    rtdata <- addscorertpred(rtdata = rtdata,
                             rtscorethr = rtscorethr,
                             rt_thr = rt_thr,
                             typeanalysis = "prediction")

  } else if (typeanalysis == "filtering") { # Filtering ........................
    if (!is.null(data$RTsec_rtpred)) {
      data <- data[ , !(names(data) %in% c("RTsec_rtpred", "CI_lower", "CI_upper",
                                           "RTsec_rtpred_prev", "diffrt", "slope",
                                           "Score_p", "Score_rtpred")) ]
    }

    # Apply RT bounds if local maxima and minima are available
    rtdata <- restrictrtpred(rtdata = rtdata,
                             typeanalysis = "filtering")

    data <- merge(data, rtdata$filtering$predictions,
                  by = c("LipidSubclass", "C", "DB"), all.x = TRUE)


    # Calculate the prediction score
    rtdata$filtering$data <- predictionscore(data = data)

    # Generate Score_rtpred
    rtdata <- addscorertpred(rtdata = rtdata,
                             rtscorethr = rtscorethr,
                             rt_thr = NULL,
                             typeanalysis = "filtering")

  } # ..........................................................................

  return(rtdata)
}

# predictionscore
#' Calculate the prediction score
#'
#' Calculate the prediction score. To do this, the distance between predictions
#' with the same C and different DB, as well as the predicted retention times, is used.
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#'
#' @return data
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
predictionscore <- function(data) {
  # (1 - (|RTtheoretical - RTobserved|)) / (diffRT_db * 2)
  data$Score_p <- round(ifelse(
    is.na(data$RTsec_rtpred), NA,
    1 - ((abs((as.numeric(data$RTsec) - as.numeric(data$RTsec_rtpred)))) /
           (as.numeric(data$diffrt*2)))), 3)
  return(data)
}

# distantprediction
#' Find the candidates that are distant from any correctly annotated lipid
#'
#' Find the candidates that are distant from any correctly annotated lipid or a
#' candidate with a Score_rtpred greater than the established threshold (parameter
#' score), from which we consider them as possible new annotations.
#'
#' @param rtdata rtdata
#' @param rtscorethr numeric. Minimum Score_rtpred value from which a prediction will be considered valid.
#' @param minC_predict numeric. Minimum number of carbon units from any valid model data or prediction beyond which predictions will be discarded. This value must be at least 3, with a default of 5.
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
distantprediction <- function(rtdata, rtscorethr, minC_predict) {
  candidates <- rtdata$prediction$candidates

  if (nrow(candidates) > 0) {
    data <- rtdata$filtering$data
    models <- rtdata$prediction$models

    if (as.numeric(minC_predict) < 3) {
      warning("The value of 'minC_predict' is less than 3. It will be automatically set to 3 to meet the minimum required by the function.\n")
      minC_predict <- 3
    }

    data <- data[!is.na(data$lipidname) & data$lipidname != "" &
                   data$checkfilter %in% c("Kept",
                                           "Kept manually by the user",
                                           "Recovered with RT modelling"), , drop = FALSE]

    # Discard predictions that are more than 5 carbon units away from any
    # valid model data or prediction (> score).
    candidates$range <- 'within the range'

    for (lipidsubclass in names(models)) {
      # Retrieve the unique values of DB for each class for which a model has been generated
      total_db <- unique(candidates$DB[candidates$LipidSubclass == lipidsubclass])

      for (db in total_db) {
        # Retrieve the candidates from the lipid class with those DB values that
        # have a Score_p greater than the specified score
        c_prediction <- unique(candidates$C[candidates$LipidSubclass == lipidsubclass &
                                              candidates$DB == db &
                                              candidates$Score_p > rtscorethr])

        # Values of C for that lipid class and DB
        c_values <- unique(data$C[data$LipidSubclass == lipidsubclass &
                                    data$DB == db])
        c_values <- c_values[order(c_values)]

        # Check which predicted lipids accepted as new annotations have C values
        # lower or higher than those already present in the initial data
        if (any(!is.na(c_values))) {
          c_left <- min(c_values, na.rm = TRUE) # Minimum C value from the model data
          # C values smaller than c_left from the accepted candidates
          c_left_values <- c(c_left, sort(c_prediction[c_prediction < c_left], decreasing = TRUE))
          c_right <- max(c_values, na.rm = TRUE) # Maximum C value from the model data
          # C values greater than c_left from the accepted candidates
          c_right_values <- c(c_right, sort(c_prediction[c_prediction > c_right], decreasing = FALSE))
        } else { # If there are no values, return NA
          c_left_values <- NA
          c_right_values <- NA
        }

        # Check from the model data towards the data with the smallest C of the
        # accepted candidates to see if they are more than 5 C units away
        if (length(c_left_values) == 1) {
          min_c <- c_left_values[1] # If there is only 1 C, it is the minimum C from the model data
        } else {
          min_c <- c_left_values[1] # If not, retain the minimum C value from the accepted candidates
          for (i in 2:length(c_left_values)) {
            # Compare each C value of the accepted candidates with the next smallest
            # C value. If they are separated by less than 6 C, consider the candidate valid
            if (abs(c_left_values[i] - min(c_left_values[i - 1], min_c, na.rm = TRUE)) <= as.numeric(minC_predict)) {
              # If the distance does not meet the allowed threshold,
              # update the value of min_c for the accepted candidates
              min_c <- c_left_values[i]
            } else {
              break
            }
          }
        }

        # Check from the data with the largest C of the accepted candidates towards
        # the model data to see if they are more than 5 C units away
        if (length(c_right_values) == 1) {
          max_c <- c_right_values[1] # If there is only 1 C, it is the maximum C from the model data
        } else {
          max_c <- c_right_values[1] # If not, retain the maximum C value from the accepted candidates
          for (i in 2:length(c_right_values)) {
            # Compare each C value of the accepted candidates with the next largest
            # C value of the candidates. If they are separated by less than 6 C units,
            # consider the candidate valid
            if (abs(c_right_values[i] - max(c_right_values[i - 1], max_c, na.rm = TRUE)) <= as.numeric(minC_predict)) {
              # If the distance does not meet the allowed threshold, update the value of min_c for the accepted candidates
              max_c <- c_right_values[i]
            } else {
              break
            }
          }
        }

        # Change "range" for the candidates that are outside of [min_c, max_c]
        out_of_range <- (
          (is.na(min_c) & is.na(max_c)) |
            (!is.na(min_c) & candidates$C < min_c) |
            (!is.na(max_c) & candidates$C > max_c)
        )

        candidates$range <- ifelse(candidates$DB == db &
                                     candidates$LipidSubclass == lipidsubclass &
                                     out_of_range, "out of range",
                                   candidates$range)
      }
    }

    # If no prediction has been made, change "range" as well
    candidates$range[is.na(candidates$RTsec_rtpred)] <- "out of range"

  } else {
    if (!"range" %in% names(candidates)) {
      candidates[["range"]] <- rep(NA, nrow(candidates))
    }
  }

  rtdata$prediction$candidates <- candidates

  return(rtdata)
}


# restrictrtpred
#' Restrict predicted retention times (RT) by local extrema and model-aware caps
#'
#' Restrict predicted retention times (RT) by local extrema and model-aware caps
#'
#' @param rtdata rtdata
#' @param typeanalysis typeanalysis
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
restrictrtpred <- function(rtdata, typeanalysis) {

  # -----------------------------
  # Argument validation
  # -----------------------------
  typeanalysis <- match.arg(typeanalysis, choices = c("filtering", "prediction"))

  if (typeanalysis == "filtering") {
    data        <- rtdata$filtering$data
    predictions <- rtdata$filtering$predictions
    models      <- rtdata$filtering$models
  } else if (typeanalysis == "prediction") {
    data        <- rtdata$prediction$candidates
    predictions <- rtdata$prediction$predictions
    models      <- rtdata$prediction$models
  }

  # -----------------------------
  # Column preparation
  # -----------------------------
  # Ensure "current" RT column exists; keep a "previous" copy to detect local extrema
  if (!"RTsec_rtpred" %in% names(data)) {
    data[["RTsec_rtpred"]] <- rep(NA, nrow(data))
  }
  if (!"RTsec_rtpred" %in% names(predictions)) {
    predictions[["RTsec_rtpred"]] <- rep(NA, nrow(predictions))
  }

  data$RTsec_rtpred_prev        <- data$RTsec_rtpred
  predictions$RTsec_rtpred_prev <- predictions$RTsec_rtpred


  # -----------------------------
  # Result data frames
  # -----------------------------
  df_max <- data.frame(
    class            = character(),
    DB               = integer(),
    C_max            = integer(),
    RTsec_rtpred_max = double(),
    distance         = double(),
    max              = character(),
    stringsAsFactors = FALSE
  )

  df_min <- data.frame(
    class             = character(),
    DB                = integer(),
    C_min             = integer(),
    RTsec_rtpred_min  = double(),
    distance          = double(),
    min               = character(),
    slope             = double(),
    stringsAsFactors  = FALSE
  )

  # -----------------------------
  # Local maxima/minima detection by class and DB
  # -----------------------------
  lipid_classes <- unique(predictions$LipidSubclass)
  lipid_classes <- lipid_classes[!is.na(lipid_classes)]

  for (cls in lipid_classes) {

    class_pred <- predictions[predictions$LipidSubclass == cls, , drop = FALSE]
    if (nrow(class_pred) == 0) next

    # Retrieve model data for the class (if available)
    class_model_data <- tryCatch(models[[cls]]$model$data, error = function(e) NULL)
    if (is.null(class_model_data) || nrow(class_model_data) == 0) next

    # Compute model median C (fallback to class predictions if needed)
    class_model_data$C <- as.numeric(class_model_data$C)
    model_C_median <- median(class_model_data$C, na.rm = TRUE)
    if (is.na(model_C_median)) {
      model_C_median <- median(class_pred$C, na.rm = TRUE)
    }
    if (is.na(model_C_median)) next # skip if no usable median

    dbs <- sort(unique(class_pred$DB))
    dbs <- dbs[!is.na(dbs)]

    for (db in dbs) {

      db_pred <- class_pred[class_pred$DB == db, , drop = FALSE]
      if (nrow(db_pred) == 0) next

      # Order by C (critical for neighbor-based extrema detection)
      ord <- order(db_pred$C)
      db_pred <- db_pred[ord, , drop = FALSE]

      # Neighbor vectors for local extrema checks
      rt_curr <- db_pred$RTsec_rtpred_prev
      rt_prev <- c(NA, rt_curr[-length(rt_curr)])
      rt_next <- c(rt_curr[-1], NA)

      # -------------------------
      # LOCAL MAXIMA
      # -------------------------
      idx_local_max <- which(
        !is.na(rt_prev) & !is.na(rt_next) &
          !is.na(rt_curr) &
          rt_curr > rt_prev & rt_curr > rt_next
      )

      if (length(idx_local_max) > 0) {
        df_max_class <- data.frame(
          class            = cls,
          DB               = db,
          C_max            = db_pred$C[idx_local_max],
          RTsec_rtpred_max = db_pred$RTsec_rtpred_prev[idx_local_max],
          stringsAsFactors = FALSE
        )
        # Distance to model median C
        df_max_class$distance <- abs(df_max_class$C_max - model_C_median)
        df_max_class$max <- "Ignore"

        # Choose the nearest local max on the RIGHT side of the model median (C > median)
        right_mask <- df_max_class$C_max > model_C_median
        if (any(right_mask, na.rm = TRUE)) {
          idx_rel <- which.min(df_max_class$distance[right_mask])
          # Map relative index back to absolute position in df_max_class
          idx_abs <- which(right_mask)[idx_rel]
          df_max_class$max[idx_abs] <- "max"
        }

        df_max <- rbind(df_max, df_max_class)
      }

      # -------------------------
      # LOCAL MINIMA
      # -------------------------
      # Re-derive neighbors (safety; db_pred is ordered by C)
      rt_curr <- db_pred$RTsec_rtpred_prev
      rt_prev <- c(NA, rt_curr[-length(rt_curr)])
      rt_next <- c(rt_curr[-1], NA)

      idx_local_min <- which(
        !is.na(rt_prev) & !is.na(rt_next) &
          !is.na(rt_curr) &
          rt_curr < rt_prev & rt_curr < rt_next
      )

      if (length(idx_local_min) > 0) {
        df_min_class <- data.frame(
          class            = cls,
          DB               = db,
          C_min            = db_pred$C[idx_local_min],
          RTsec_rtpred_min = db_pred$RTsec_rtpred_prev[idx_local_min],
          stringsAsFactors = FALSE
        )
        # Distance to model median C
        df_min_class$distance <- abs(df_min_class$C_min - model_C_median)
        df_min_class$min <- "Ignore"
        df_min_class$slope <- NA_real_

        # Choose the nearest local min on the LEFT side of the model median (C < median)
        left_mask <- df_min_class$C_min < model_C_median
        if (any(left_mask, na.rm = TRUE)) {
          idx_rel <- which.min(df_min_class$distance[left_mask])
          idx_abs <- which(left_mask)[idx_rel]
          df_min_class$min[idx_abs] <- "min"

          # Compute slope at selected min using the C+1 point from the same db_pred
          c0 <- df_min_class$C_min[idx_abs]
          rt0 <- df_min_class$RTsec_rtpred_min[idx_abs]

          idx_cplus <- which(db_pred$C == (c0 + 1))
          if (length(idx_cplus) >= 1) {
            # If multiple isomers with the same C exist, take the first deterministically
            rt_cplus <- db_pred$RTsec_rtpred_prev[idx_cplus[1]]
            df_min_class$slope[idx_abs] <- rt_cplus - rt0
          }
        }

        df_min <- rbind(df_min, df_min_class)
      }

    }
  }

  # -----------------------------
  # Classes affected by local extrema restrictions  (NUEVO)
  # -----------------------------
  affected_classes <- unique(c(
    df_max$class[df_max$max == "max"],
    df_min$class[df_min$min == "min"]
  ))
  affected_classes <- affected_classes[!is.na(affected_classes)]

  # -----------------------------
  # Apply UPPER CAP from selected maxima (VIEJO)
  # -----------------------------
  if (nrow(df_max) > 0) {
    chosen_max <- df_max[df_max$max == "max", , drop = FALSE]
    if (nrow(chosen_max) > 0) {
      for (cls in unique(chosen_max$class)) {
        chosen_cls <- chosen_max[chosen_max$class == cls, , drop = FALSE]
        for (db in unique(chosen_cls$DB)) {
          chosen_db <- chosen_cls[chosen_cls$DB == db, , drop = FALSE]
          if (nrow(chosen_db) == 0) next

          # If multiple 'max' exist per DB (should not), pick the first deterministically
          chosen_c  <- as.numeric(chosen_db$C_max[1])
          chosen_rt <- as.numeric(chosen_db$RTsec_rtpred_max[1])

          # predictions: cap for C > chosen_c
          idx_pred <- which(predictions$LipidSubclass == cls &
                              predictions$DB == db &
                              !is.na(predictions$C) &
                              predictions$C > chosen_c &
                              !is.na(predictions$RTsec_rtpred_prev))
          if (length(idx_pred) > 0) {
            predictions$RTsec_rtpred[idx_pred] <- chosen_rt
          }

          # data: cap for C > chosen_c
          idx_cand <- which(data$LipidSubclass == cls &
                              data$DB == db &
                              !is.na(data$C) &
                              data$C > chosen_c &
                              !is.na(data$RTsec_rtpred_prev))
          if (length(idx_cand) > 0) {
            data$RTsec_rtpred[idx_cand] <- chosen_rt
          }
        }
      }
    }
  }

  # -----------------------------
  # Apply LOWER DIAGONAL CAP from selected minima
  # -----------------------------
  if (!"slope" %in% names(predictions)) predictions$slope <- NA

  mins <- df_min[df_min$min == "min", , drop = FALSE]

  if (nrow(mins) > 0) {
    if (!"slope" %in% names(mins)) mins$slope <- NA # ensure column presence

    for (i in seq_len(nrow(mins))) {
      cls   <- mins$class[i]
      db    <- mins$DB[i]
      cmin  <- as.numeric(mins$C_min[i])
      rtmin <- as.numeric(mins$RTsec_rtpred_min[i])

      # 1) Prefer slope stored in df_min
      slope <- as.numeric(mins$slope[i])

      # 2) If missing, recompute slope from predictions at C = cmin + 1
      if (is.na(slope)) {
        idx_cplus <- which(predictions$LipidSubclass == cls &
                             predictions$DB == db &
                             !is.na(predictions$C) &
                             predictions$C == (cmin + 1) &
                             !is.na(predictions$RTsec_rtpred_prev))
        if (length(idx_cplus) >= 1) {
          idx_cplus <- idx_cplus[1]
          slope <- as.numeric(predictions$RTsec_rtpred_prev[idx_cplus] - rtmin)
        } else {
          # No C+1 available to estimate slope: skip this min
          next
        }
      }

      # 3) Apply descending diagonal cap to 'predictions' for C <= cmin
      idx_pred <- which(predictions$LipidSubclass == cls &
                          predictions$DB == db &
                          !is.na(predictions$C) &
                          predictions$C <= cmin)

      if (length(idx_pred) > 0) {
        deltaC_pred <- cmin - predictions$C[idx_pred]   # 0,1,2,...
        predictions$RTsec_rtpred[idx_pred] <- rtmin - slope * deltaC_pred
        predictions$slope[idx_pred] <- slope
      }

      # 4) Apply consistently to 'data' for C <= cmin (length-aligned)
      idx_cand <- which(data$LipidSubclass == cls &
                          data$DB == db &
                          !is.na(data$C) &
                          data$C <= cmin &
                          !is.na(data$RTsec_rtpred_prev))

      if (length(idx_cand) > 0) {
        deltaC_data <- cmin - data$C[idx_cand]          # 0,1,2,...
        data$RTsec_rtpred[idx_cand] <- rtmin - slope * deltaC_data
      }

    }
  }


  # -----------------------------
  # Recalculate diffrt for affected classes  (NUEVO)
  # -----------------------------
  if (length(affected_classes) > 0) {

    if (!"diffrt" %in% names(predictions)) {
      predictions$diffrt <- NA_real_
    }

    for (cls in affected_classes) {

      pred_cls <- predictions[
        predictions$LipidSubclass == cls,
        ,
        drop = FALSE
      ]

      if (nrow(pred_cls) == 0) next

      for (c in unique(pred_cls$C)) {

        subset_data <- pred_cls[
          pred_cls$C == c,
          ,
          drop = FALSE
        ]

        subset_data <- subset_data[
          order(subset_data$DB),
          ,
          drop = FALSE
        ]

        if (nrow(subset_data) < 2) {

          predictions$diffrt[
            predictions$LipidSubclass == cls &
              predictions$C == c
          ] <- NA_real_

          next
        }

        distances <- abs(diff(as.numeric(subset_data$RTsec_rtpred))) /
          diff(as.numeric(subset_data$DB))

        predictions$diffrt[
          predictions$LipidSubclass == cls &
            predictions$C == c
        ] <- mean(distances, na.rm = TRUE)
      }
    }
  }

  # -----------------------------
  # Write-back to rtdata          (VIEJO)
  # -----------------------------
  if (typeanalysis == "filtering") {
    rtdata$filtering$data        <- data
    rtdata$filtering$predictions <- predictions
    rtdata$filtering$restrictrtpred$max  <- df_max
    rtdata$filtering$restrictrtpred$min  <- df_min
  } else { # "prediction"
    rtdata$prediction$candidates  <- data
    rtdata$prediction$predictions <- predictions
    rtdata$prediction$restrictrtpred$max  <- df_max
    rtdata$prediction$restrictrtpred$min  <- df_min
  }

  return(rtdata)
}


# addscorertpred
#' Add new column: Score_rtpred
#'
#' Add a new column: Score_rtpred. Initially, calculate the score_p, but since
#' this value can range from -Inf to 1 depending on the predicted RTs, generate
#' a new variable, Score_rtpred, to scale the values to the range [0,1], also
#' considering the distant predictions
#'
#' @param rtdata rtdata
#' @param rtscorethr numeric. Minimum Score_rtpred value from which a prediction or will be considered valid.
#' @param rt_thr numeric. Threshold for the difference (in seconds) between the
#' observed retention time (RT) and the predicted retention time (RTsec_rtpred).
#' @param typeanalysis Character. Specifies whether the function is being used for
#' filtering or RT predictions. Possible values are "filtering" or "prediction".
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
addscorertpred <- function(rtdata, rtscorethr, rt_thr = 30, typeanalysis) {

  typeanalysis <- match.arg(typeanalysis, choices = c("filtering", "prediction"))

  if (typeanalysis == "filtering") { # Filtering ...............................

    data <- rtdata$filtering$data

    # Score_rt based on score_p
    data$Score_rtpred <- data$Score_p
    data$Score_rtpred[is.na(data$Score_p)] <- 0
    data$Score_rtpred[data$Score_p < 0 | is.na(data$Score_p)] <- 0

    if (is.null(data$checkfilter)) {
      data$checkfilter <- data$Filter
    }

    if (is.null(data$checkmodel)) {
      data$checkmodel <- data$checkfilter
    }

    # Recovered with RT modelling
    # data$checkfilter[(
    #   data$Filter %in% c("Confidence level among isobars", "Removed distant data","Removed by RT rules",
    #                      "In-source fragments", "Outlier", "RT coherence between subclasses") &
    #     !(data$checkfilter %in% c("Kept manually by the user", "Removed manually by the user")) &
    #     !is.na(data$Score_p) & data$Score_p > rtscorethr)] <- "Recovered with RT modelling"

    data$checkmodel[(
      data$Filter %in% c("Confidence level among isobars", "Removed distant data","Removed by RT rules",
                         "In-source fragments", "Outlier") &
        !(data$checkfilter %in% c("Kept manually by the user", "Removed manually by the user")) &
        !is.na(data$Score_p) & data$Score_p > rtscorethr &
        !grepl("P[A-Z]-(O|P)|LP[A-Z]-(O|P))$", data$LipidSubclass)
    )] <- "Recovered with RT modelling"

    data$checkmodel[data$checkmodel == "Recovered with RT modelling" &
                      !is.na(data$Score_p) & data$Score_p < rtscorethr] <-
      data$Filter[data$checkmodel == "Recovered with RT modelling" &
                    !is.na(data$Score_p) & data$Score_p < rtscorethr]


    data$checkfilter[(
      data$Filter %in% c("Confidence level among isobars", "Removed distant data","Removed by RT rules",
                         "In-source fragments", "Outlier") &
        !(data$checkfilter %in% c("Kept manually by the user", "Removed manually by the user")) &
        !is.na(data$Score_p) & data$Score_p > rtscorethr &
        !grepl("P[A-Z]-(O|P)|LP[A-Z]-(O|P))$", data$LipidSubclass))] <- "Recovered with RT modelling"

    data$checkfilter[data$checkfilter == "Recovered with RT modelling" &
                      !is.na(data$Score_p) & data$Score_p < rtscorethr] <-
      data$Filter[data$checkfilter == "Recovered with RT modelling" &
                    !is.na(data$Score_p) & data$Score_p < rtscorethr]


    # Removed with RT modelling
    # data$checkmodel[(
    #   !data$Filter %in% c("Confidence level among isobars", "Removed distant data","Removed by RT rules",
    #                       "In-source fragments", "Outlier", "RT coherence between subclasses") &
    #     !(data$checkfilter %in% c("Kept manually by the user", "Removed manually by the user")) &
    #     !is.na(data$Score_p) & data$Score_p <= 0)] <- "Removed with RT modelling"
    #
    # data$checkmodel[data$checkmodel == "Removed with RT modelling" &
    #                    !is.na(data$Score_p) & data$Score_p > 0] <-
    #   data$Filter[data$checkmodel == "Removed with RT modelling" &
    #                 !is.na(data$Score_p) & data$Score_p > 0]
    #
    # data$checkfilter[(
    #   !data$Filter %in% c("Confidence level among isobars", "Removed distant data","Removed by RT rules",
    #                       "In-source fragments", "Outlier", "RT coherence between subclasses") &
    #     !(data$checkfilter %in% c("Kept manually by the user", "Removed manually by the user")) &
    #     !is.na(data$Score_p) & data$Score_p <= 0)] <- "Removed with RT modelling"
    #
    # data$checkfilter[data$checkfilter == "Removed with RT modelling" &
    #                   !is.na(data$Score_p) & data$Score_p > 0] <-
    #   data$Filter[data$checkfilter == "Removed with RT modelling" &
    #                 !is.na(data$Score_p) & data$Score_p > 0]


    rtdata$filtering$data <- data

    return(rtdata)

  } else if (typeanalysis == "prediction") { # Prediction ......................

    data <- rtdata$filtering$data
    candidates <- rtdata$prediction$candidates

    if (nrow(candidates) > 0) {
      # Score_rt based on score_p and range
      candidates$Score_rtpred <- candidates$Score_p
      candidates$Score_rtpred[candidates$Score_p < 0 | is.na(candidates$Score_p)] <- 0
      candidates$Score_rtpred[is.na(candidates$Score_p)] <- 0
      # candidates$Score_rtpred[candidates$range == "out of range"] <- 0

      # order by group
      candidates <- candidates[order(candidates$LipidRTid), ]

      # Score_rtpred for manualRTprediction
      candidates$checkscore_rtpred <- ifelse(candidates$Score_rtpred > rtscorethr,
                                             "Kept", "Removed")

      candidates$checkscore_rtpred <- ifelse((candidates$checkscore_rtpred == "Removed" &
                                                abs(candidates$RTsec - candidates$RTsec_rtpred) <= rt_thr),
                                             "Discarded",
                                             candidates$checkscore_rtpred)

      if (is.null(candidates$checkcandidate)) {
        candidates$checkcandidate <- candidates$checkscore_rtpred
      }

      # removal_types <- c("Confidence level among isobars", "Outlier", #"RT coherence between plasmalogens"
      #                    "RT coherence between subclasses")
      removal_types <- c()
      for (removal_type in removal_types) {
        # Get the groups that contain that type of filtered data
        fragments_group <- unique(data$LipidRTid[data$checkfilter == removal_type])
        for (i in fragments_group) {
          # Get the lipids from the current group marked with the corresponding filtering
          fragments <- data[data$LipidRTid == i & data$checkfilter == removal_type, ]
          for (j in 1:nrow(fragments)) {
            # Mark as "Removed" in candidates if they match in group, LipidClass, C, and DB
            candidates$checkscore_rtpred[
              candidates$LipidRTid == i &
                candidates$LipidSubclass == fragments$LipidSubclass[j] &
                candidates$C == fragments$C[j] &
                candidates$DB == fragments$DB[j]] <- "Removed"
          }
        }
      }

      # generate indices for the candidates
      candidates$index <- seq(max(data$index) + 1, max(data$index) + nrow(candidates))

    } else {
      if (!"Score_rtpred" %in% names(candidates)) {
        candidates[["Score_rtpred"]] <- rep(NA, nrow(candidates))
      }
      if (!"checkscore_rtpred" %in% names(candidates)) {
        candidates[["checkscore_rtpred"]] <- rep(NA, nrow(candidates))
      }
      if (!"checkcandidate" %in% names(candidates)) {
        candidates[["checkcandidate"]] <- rep(NA, nrow(candidates))
      }
      if (!"index" %in% names(candidates)) {
        candidates[["index"]] <- rep(NA, nrow(candidates))
      }
    }

    rtdata$prediction$candidates <- candidates

    return(rtdata)
  }
}


# readfiledbs
#' Allows reading a CSV that contains lipid classes and adducts to perform
#' annotations and retention time predictions.
#'
#' Allows reading a CSV that contains lipid classes and adducts to perform
#' annotations and retention time predictions.
#'
#' @param filedbs csv with classes to predict
#' @param sep Character. Field separator used when reading the CSV file.
#' @param adducts Data.frame. modify adductsTable.
#'
#' @return filedbs data frame
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
readfiledbs <- function(filedbs = NULL, sep = ",", adducts = NULL) {

  # Build default file DBs based on the given polarity
  filedbs_default <- defaultfiledbs(sep = sep)

  if (!is.null(filedbs)) {

    if (is.character(filedbs)) {  # CSV path provided

      # Read CSV
      filedbs_csv <- read.csv(filedbs, header = TRUE, sep = sep,
                              stringsAsFactors = FALSE, check.names = TRUE,
                              na.strings = c("", "NA"))

      # Replace all NAs by empty strings
      filedbs_csv[is.na(filedbs_csv)] <- ""

      # Adducts columns (at least one must be present)
      adducts_cols <- c("adducts_pos", "adducts_neg")
      present_adducts <- adducts_cols[adducts_cols %in% colnames(filedbs_csv)]

      # Error if none of the adducts columns exist
      if (length(present_adducts) == 0) {
        stop("At least one of 'adducts_pos' or 'adducts_neg' must be present in the data frame.")
      }

      # Common columns (mandatory + present adducts)
      required_cols <- c("classes", present_adducts, "fas", "sph", "alkyl", "vinyl")
      cols_trim     <- c(present_adducts, "fas", "sph", "alkyl", "vinyl")

      # Validate required columns
      missing_cols <- setdiff(required_cols, colnames(filedbs_csv))
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
      }

      # Keep only required columns, in order
      filedbs_csv <- filedbs_csv[, required_cols, drop = FALSE]


      # Keep only rows where at least one adduct column has a non-missing, non-empty value
      filedbs_csv <- filedbs_csv[
        apply(filedbs_csv[, present_adducts, drop = FALSE], 1, function(row) {
          any(!is.na(row) & grepl("[^\\s;]", row))
        }),
      ]


      # Trim and normalize whitespace across key string columns
      normalize_spaces <- function(x) gsub("\\s+", "", x)
      filedbs_csv[cols_trim] <- lapply(filedbs_csv[cols_trim],function(x) gsub("\\s+", "", x))

      # Restrict to supported classes based on defaults
      supported_classes <- unique(filedbs_default$classes)
      classes_to_remove <- setdiff(unique(filedbs_csv$classes), supported_classes)
      filedbs_csv <- filedbs_csv[filedbs_csv$classes %in% supported_classes, , drop = FALSE]

      if (length(classes_to_remove)) {
        message(
          "The following classes are not supported and will be ignored: ",
          paste(classes_to_remove, collapse = ", ")
        )
      }

      # Compare chain specs (fas/sph/alkyl/vinyl) between CSV and defaults to detect diffs
      cols_check <- c("fas", "sph", "alkyl", "vinyl")

      # Prepare minimal frames for merging
      left_df  <- filedbs_csv[, c("classes", cols_check), drop = FALSE]
      right_df <- filedbs_default[, c("classes", cols_check), drop = FALSE]

      # Merge by 'classes' (left join) to align CSV vs defaults
      df_cmp <- merge(left_df, right_df, by = "classes",
                      suffixes = c("_csv", "_default"), all.x = TRUE)

      # Define a normalized equality check (ignores whitespace and treats NA/"")
      norm <- function(x) {
        x <- ifelse(is.na(x), "", x)
        gsub("\\s+", "", x)
      }

      # Compute a logical matrix of differences across the chain columns
      diff_matrix <- sapply(cols_check, function(col)
        norm(df_cmp[[paste0(col, "_csv")]]) != norm(df_cmp[[paste0(col, "_default")]])
      )
      # Ensure it's a 2D matrix even if length(cols_check) == 1
      if (is.null(dim(diff_matrix))) diff_matrix <- matrix(diff_matrix, ncol = 1)

      # Any row with at least one differing chain column
      diff_rows <- rowSums(diff_matrix) > 0

      # Subset only rows with any change and keep 'classes' and default dbs key
      df_diff <- df_cmp[diff_rows, c("classes",
                                     paste0(cols_check, "_csv")), drop = FALSE]

      # Start from default DBs but limited to CSV-declared dbs
      dbs <- assignDBs()
      dbs <- dbs[names(dbs) %in% unique(c(filedbs_csv$classes, "adductsTable"))]

      if (!is.null(adducts)) {
        if (is.data.frame(adducts)) {
          if (all(c("adduct", "mdiff", "charge", "n") %in% names(dbs$adductsTable))) {
            dbs$adductsTable <- adducts
          } else {
            warning("Missing columns in adducts data frame. Must include 'adduct', 'mdiff', 'charge', 'n'\n")
          }
        }
      }

      # Helper to split chain strings safely into character vectors
      split_chains <- function(s) {
        s <- trimws(s)
        if (!nzchar(s)) return(character(0))
        # If sep_adducts_chains could be multi-char, avoid fixed=TRUE; otherwise set fixed=TRUE for speed
        out <- strsplit(s, "\\|", fixed = FALSE)[[1]]
        trimws(out[nzchar(out)])
      }

      # Rebuild DBs for classes that have changes in chain definitions
      if (nrow(df_diff) > 0) {
        for (n in seq_len(nrow(df_diff))) {
          lipid_class <- df_diff$classes[[n]]

          fas_v   <- split_chains(df_diff$fas_csv[[n]])
          sph_v   <- split_chains(df_diff$sph_csv[[n]])
          alkyl_v <- split_chains(df_diff$alkyl_csv[[n]])
          vinyl_v <- split_chains(df_diff$vinyl_csv[[n]])

          # Build the class-specific DB with updated chains and replace in global list
          built <- createLipidDBs(
            lipid = lipid_class,
            chains = fas_v,
            sph    = sph_v,
            alkyl  = alkyl_v,
            vinyl  = vinyl_v)

          # Replace only the corresponding element
          if (!is.null(built[[lipid_class]])) {
            dbs[[lipid_class]] <- built[[lipid_class]]
          } else {
            warning("Expected DBs '", lipid_class, "' was not produced.")
          }
        }
      }

      # Verify that the adducts match those defined in dbs$adductsTable
      lipidClasses <- list()
      warning_msg <- character()

      for (i in seq_len(nrow(filedbs_csv))) {
        class_name <- filedbs_csv$classes[i]

        # Extract all adducts present according to the detected columns
        adducts_raw <- unlist(
          lapply(present_adducts, function(col) {
            # Split by sep_adducts_chains and trim spaces
            trimws(strsplit(filedbs_csv[[col]][i], "\\|")[[1]])
          })
        )

        # Filter only valid adducts according to the database
        valid_adducts <- adducts_raw[adducts_raw %in% dbs$adductsTable$adduct]

        # Detect adducts that are not valid
        invalid_adducts <- setdiff(adducts_raw, valid_adducts)

        # Keep only the valid adducts in the final list
        lipidClasses[[class_name]] <- valid_adducts

        # Clean the original columns of the row so that only the valid ones remain
        for (col in present_adducts) {
          if (!is.na(filedbs_csv[[col]][i]) && nchar(filedbs_csv[[col]][i]) > 0) {
            # Split by sep_adducts_chains, remove spaces, filter valid ones, and join again
            adducts_clean <- trimws(strsplit(filedbs_csv[[col]][i], "\\|")[[1]])
            adducts_clean <- adducts_clean[adducts_clean %in% dbs$adductsTable$adduct]
            filedbs_csv[[col]][i] <- paste(adducts_clean, collapse = paste0(gsub("\\\\", "", "\\|"), " "))
          }
        }

        # Generate warning message if invalid adducts were removed
        if (length(invalid_adducts) > 0) {
          msg <- paste0(
            "   - '", class_name, "': ",
            paste(invalid_adducts, collapse = "\\|"))

          warning_msg <- c(warning_msg, msg)
        }
      }


      if (length(warning_msg) > 0) {

        # Valid adducts by polarity
        valid_pos <- dbs$adductsTable$adduct[dbs$adductsTable$charge > 0]
        valid_neg <- dbs$adductsTable$adduct[dbs$adductsTable$charge < 0]

        warning_msg <- c(
          "Warning message:",
          "The following adduct(s) were not valid and have been removed:",
          warning_msg
        )

        # Create complete message
        full_msg <- paste(
          c(warning_msg,
            paste0("\nAllowed adducts for positive polarity: ", paste(unique(valid_pos), collapse = ", ")),
            paste0("Allowed adducts for negative polarity: ", paste(unique(valid_neg), collapse = ", ")),
            "\nNote: Only valid adducts will be used. To add additional valid adducts, use the 'adducts' parameter of the function."
          ),
          collapse = "\n"
        )

        message(full_msg)
      }


      filedbs_final <- filedbs_csv

    } else {
      stop("'filedbs' must be a path (character) to a CSV file or NULL.")
    }

  } else {
    # No CSV provided: use defaults dbs
    dbs <- assignDBs()

    if (!is.null(adducts)) {
      if (is.data.frame(adducts)) {
        if (all(c("adduct", "mdiff", "charge", "n") %in% names(dbs$adductsTable))) {
          dbs$adductsTable <- adducts
        } else {
          warning("Missing columns in adducts data frame. Must include 'adduct', 'mdiff', 'charge', 'n'\n")
        }
      }
    }

    filedbs_final <- filedbs_default
  }

  return(list(filedbs = filedbs_final, dbs = dbs))
}


# createadductsmzdbs
#' Allows reading a CSV that contains lipid classes and adducts to perform
#' annotations and retention time predictions.
#'
#' Allows reading a CSV that contains lipid classes and adducts to perform
#' annotations and retention time predictions.
#'
#' @param filedbs csv with classes to predict
#' @param sep Character. Field separator used when reading the CSV file.
#' @param adducts Data.frame. modify adductsTable.
#' @param classes vector.
#'
#' @return adductsmzdbs list
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
createadductsmzdbs <- function(filedbs = NULL, sep = ",", adducts = NULL, classes = NULL) {

  # Build default file DBs based on the given polarity
  if (!is.null(filedbs)) {
    read_filedbs <- readfiledbs(filedbs = filedbs,
                                sep = sep,
                                adducts = adducts)

    filedbs <- read_filedbs$filedbs
    dbs <- read_filedbs$dbs
  } else {
    # Build default file DBs based on the given polarity
    filedbs <- defaultfiledbs(sep = sep)
    dbs <- assignDBs()
  }


  if (is.data.frame(filedbs) & length(dbs) > 0) {

    # Once the file has been cleaned and the base dbs object has been updated,
    # add the columns containing the m/z values of the adducts
    # Ensure the database object is valid and contains the required 'AdductsTable'
    if (!is.list(dbs) || is.null(names(dbs)) || !("adductsTable" %in% names(dbs))) {
      stop(
        "Missing required resource: 'adductsTable'. ",
        "The database must include the 'adductsTable' in order to propose candidates. ",
        "Please verify that your database is complete and try again."
      )
    }

    # Also ensure the table is not NULL or empty
    if (is.null(dbs[["adductsTable"]]) || nrow(dbs[["adductsTable"]]) == 0) {
      stop(
        "'adductsTable' is present but empty or invalid. ",
        "The database must include a valid 'adductsTable' to propose candidates."
      )
    }
  }


  db_adducts <- dbs[["adductsTable"]]
  adductsmzdbs <- c()
  for (i in 1:nrow(filedbs)) {
    lipidclass <- filedbs$classes[[i]]

    # Combine all available adducts for the row i
    adducts <- unlist(
      lapply(c("adducts_pos", "adducts_neg"), function(col) {
        if (col %in% colnames(filedbs)) {
          trimws(strsplit(filedbs[[col]][i], "\\|")[[1]])
        } else {
          NULL
        }
      })
    )

    db <- dbs[[lipidclass]]
    if (!is.data.frame(db) && !(is.list(db) && length(db) == 0)) {
      db <- db[[1]]
    } else if (is.list(db) && length(db) == 0) {
      next
    }

    mz_adducts <- c()
    for (a in 1:length(adducts)) {
      adduct <- adducts[a]
      db$x <- db$Mass + ((db_adducts[c(db_adducts$adduct==adduct),"mdiff"] *
                            abs(db_adducts[c(db_adducts$adduct==adduct),"n"]))/
                           abs(db_adducts[c(db_adducts$adduct==adduct),"charge"]))
      names(db)[names(db) == "x"] <- adduct
    }
    adductsmzdbs[[lipidclass]] <- db
  }

  adductsmzdbs[["adductsTable"]] <- db_adducts

  return(list(filedbs = filedbs, adductsmzdbs = adductsmzdbs))
}


# proposertcandidates
#' Propose candidates for the features without annotations
#'
#' Propose candidates for the features without annotations from the classes of
#' interest, based on the mz of the features, considering the mz tolerance value
#' provided by the user. This value can be given in ppm (parts per million) or
#' absolute mass error
#'
#' @param rtdata rtdata
#' @param filedbs csv with classes to predict.
#' @param sep Character. Field separator used when reading the CSV file.
#' @param classes A vector containing the lipid classes for which retention time
#' is to be predicted. If classes = NULL, predictions will be made for all classes
#' where there are candidates and data to create a model.
#' @param adducts Data.frame. modify adductsTable.
#' @param dmz m/z tolerance in ppm or absolute mass error.
#' @param ppm logical. If TRUE, m/z tolerance in ppm. If FALSE, m/z tolerance in
#' absolute mass error.
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
proposertcandidates <- function(rtdata,  filedbs = NULL, sep = ",", classes = NULL,
                                adducts = NULL, dmz = 5, ppm = TRUE){

  data_isNA <- rtdata$data[(is.na(rtdata$data$lipidname) |
                              rtdata$data$lipidname == "" |
                              !(rtdata$data$Filter %in% c("Kept", "Ignore duplicated", "Unsupported class"))) &
                             !is.na(rtdata$data$mz) &
                             !is.na(rtdata$data$RTsec),
                           c("LipidRTid", "mz", "RTsec"), drop = FALSE]

  # Identify LipidRTid that must be excluded
  exclude_ids <- unique(
    rtdata$filtering$data$LipidRTid[
      rtdata$filtering$data$checkfilter %in%
        c("Kept", "Kept manually by the user", "Recovered with RT modelling")
    ]
  )

  # Keep ONLY LipidRTid that never appear in those categories
  filter_isNA <- rtdata$filtering$data[
    !(rtdata$filtering$data$LipidRTid %in% exclude_ids),
    c("LipidRTid", "mz", "RTsec"),
    drop = FALSE
  ]

  isNA <- unique(rbind(data_isNA, filter_isNA))


  polarity <- rtdata$metadata$polarity

  if (nrow(isNA) != 0 && polarity != "both" && !is.null(polarity)) {
    # Obtain the m/z LipidMS database for the classes on which you want to predict
    generateadductsmzdbs <- createadductsmzdbs(filedbs = filedbs,
                                               sep = sep,
                                               classes = classes,
                                               adducts = adducts)
    filedbs <- generateadductsmzdbs$filedbs
    adductsmzdbs <- generateadductsmzdbs$adductsmzdbs

    if (is.character(adductsmzdbs)) { # If the provided classes are not valid
      return(adductsmzdbs)
    } else { # If the provided classes are valid
      # Propose candidates
      candidates <- data.frame()


      # Valid adducts according to analysis polarity
      if (polarity == "positive") {
        valid_adducts <- adductsmzdbs$adductsTable$adduct[adductsmzdbs$adductsTable$charge > 0]
      } else if (polarity == "negative") {
        valid_adducts <- adductsmzdbs$adductsTable$adduct[adductsmzdbs$adductsTable$charge < 0]
      }
      # invalid_adducts, adducts in adductsTable but not in valid_adducts
      invalid_adducts <- adductsmzdbs$adductsTable$adduct[!adductsmzdbs$adductsTable$adduct %in% valid_adducts]


      rtdata$metadata$prediction$dbs <- adductsmzdbs
      rtdata$metadata$prediction$filedbs <- filedbs


      # For each feature without an annotation
      for (LipidRTid in unique(as.numeric(isNA$LipidRTid))) {

        mz <- unique(as.numeric(isNA$mz[isNA$LipidRTid == LipidRTid]))
        rt <- unique(as.numeric(isNA$RTsec[isNA$LipidRTid == LipidRTid]))

        # For each lipid class m/z table
        # check if any m/z matches the feature's m/z within the +- tolerance value
        for (LipidSubclass.cand in setdiff(names(adductsmzdbs), "adductsTable")) {
          mz_db <- adductsmzdbs[[LipidSubclass.cand]]

          adduct_columns <- names(mz_db[!(names(mz_db) %in% c("formula", "total", "Mass", "ID", invalid_adducts))])

          if (ppm == TRUE) { # If the tolerance value is given in ppm
            matches <- apply(mz_db[adduct_columns], 2,
                             function(x) (((abs(as.numeric(x) - mz))/x)*1e6 <= as.numeric(dmz)))
          } else { # If the tolerance value is given in absolute mass error
            matches <- apply(mz_db[adduct_columns], 2,
                             function(x) (abs(as.numeric(x) - mz) <= as.numeric(dmz)))
          }

          # If any m/z matches, store the lipid along with its adduct as a candidate
          if (length(which(matches)) > 0) {
            for (cand in which(matches)) {
              CDB <- rep(mz_db$total, length(adduct_columns))[cand]

              if (grepl("_d$", LipidSubclass.cand)) {
                lipidname.cand <- paste0(gsub("_d$", "", LipidSubclass.cand),
                                         "(d", CDB, ")")
              } else if (grepl("_t$", LipidSubclass.cand)) {
                lipidname.cand <- paste0(gsub("_t$", "", LipidSubclass.cand),
                                         "(t", CDB, ")")
              } else if (grepl("-O$", LipidSubclass.cand)) {
                lipidname.cand <- paste0(gsub("-O$", "", LipidSubclass.cand),
                                         "(O-", CDB, ")")
              } else if (grepl("-P$", LipidSubclass.cand)) {
                lipidname.cand <- paste0(gsub("-P$", "", LipidSubclass.cand),
                                         "(P-", CDB, ")")
              } else {
                lipidname.cand <- paste0(LipidSubclass.cand, "(", CDB, ")")
              }

              # remove any spaces just in case
              lipidname.cand <- gsub("\\s+", "", lipidname.cand)
              LipidClass.cand <- sub("_.*", "", LipidSubclass.cand)
              LipidClass.cand <- sub("-.*", "", LipidSubclass.cand)
              LipidClass.cand <- gsub("(_d|_t|-O|-P)", "", sub("-.*", "", LipidSubclass.cand))

              # Esto funciona asi ???
              lipidname.cand <- convertlipidannotation(lipidname.cand, LipidClass.cand, LipidSubclass.cand)

              C.cand <- as.numeric(unlist(strsplit(CDB, ":"))[1])
              DB.cand <- as.numeric(unlist(strsplit(CDB, ":"))[2])
              adducts.cand <- adduct_columns[ceiling(cand/dim(mz_db)[1])]

              new_row <- data.frame("mz" = mz,
                                    "RTsec" = rt,
                                    "LipidRTid" = LipidRTid,
                                    "lipidname" = lipidname.cand,
                                    "LipidClass" = LipidClass.cand,
                                    "LipidSubclass" = LipidSubclass.cand,
                                    "C" = C.cand,
                                    "DB" = DB.cand,
                                    "adduct" = adducts.cand)
              candidates <- rbind(candidates, new_row)
            }
          }
        }
      }

      if (nrow(candidates) != 0) {
        candidates <- candidates[order(candidates$C, candidates$DB), , drop = FALSE]

        # Add a column with the grouped lipid classes to easily separate them when building models
        candidates <- grouplipidclasses(candidates, map_type = "class")

        rtdata$prediction$candidates <- candidates

      } else {
        # empty data frame
        candidates <- data.frame(
          C = integer(0),
          DB = integer(0),
          mz = numeric(0),
          RTsec = numeric(0),
          LipidRTid = integer(0),
          lipidname = character(0),
          LipidClass = character(0),
          LipidSubclass = character(0),
          adduct = character(0),
          groupClass = character(0)
        )

        rtdata$prediction$candidates <- candidates
        warning("No candidate lipids found for features that were not already annotated.\n")
      }
    }

  } else {
    # empty data frame
    candidates <- data.frame(
      C = integer(0),
      DB = integer(0),
      mz = numeric(0),
      RTsec = numeric(0),
      LipidRTid = integer(0),
      lipidname = character(0),
      LipidClass = character(0),
      LipidSubclass = character(0),
      adduct = character(0),
      groupClass = character(0)
    )

    rtdata$prediction$candidates <- candidates

    # Valid adducts according to analysis polarity
    rtdata$metadata$prediction$dbs <- NULL
    rtdata$metadata$prediction$filedbs <- NULL

    if (is.null(polarity) | polarity == "both") { # Target
      warning(
        paste("Adducts must correspond to either positive or negative polarity only.\n",
              "Please ensure all adducts belong to a single polarity group.\n",
              "\n - No valid polarity could be determined. Consequently, no candidate
              lipids will be proposed and retention time predictions will not be performed.\n"))

    }
    if (nrow(isNA) == 0) {
      warning(
        paste("\n - No unannotated features were identified for candidate generation.
              Consequently, no additional lipid candidates will be proposed.\n"))
    }
  }

  return(rtdata)
}

# createpredictgraph
#' Plot the retention time models and predictions for the different lipid classes
#'
#' Plot the retention time models and predictions for the different lipid classes
#'
#' @param data_graph the annotated peak list that will be used for generating the graph.
#' @param min_c Numeric. Minimum carbon number.
#' @param max_c Numeric. Maximum carbon number.
#' @param min_rt Numeric. Minimum retention time.
#' @param max_rt Numeric. Maximum retention time.
#' @param title character. graph title
#' @param subtitle character. graph subtitle
#' @param prediction_layer predictions for each candidate (prediction$predictions)
#' @param shapes_present character. shape of the data present for plotting
#' @param fills_present character. fills of the data present for plotting
#'
#' @return graph
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
createpredictgraph <- function(data_graph, min_c, max_c, min_rt, max_rt, title,
                               subtitle = NULL, prediction_layer = NULL,
                               shapes_present, fills_present) {

  C <- RTsec <- DB <- RTsec_rtpred <- legend_label <- comments <- c()

  graph <- ggplot() +
    prediction_layer +
    geom_point(data = data_graph, aes(x = C, y = RTsec, color = factor(DB), shape = legend_label,
                                      fill = legend_label), size = 2.5, stroke = 1.1) +
    labs(title = title, x = "Carbon number", y = "Retention time (s)") +
    theme_minimal() +
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
    scale_shape_manual(values = shapes_present) +
    scale_fill_manual(values = fills_present) +
    theme(legend.position = "right", legend.box.spacing = unit(1, "lines"),
          legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
    guides(shape = guide_legend(title = "Predictions", title.position = "top",
                                override.aes = list(fill = fills_present, stroke = 1.1)),
           fill = guide_legend(title = "Predictions", title.position = "top",
                               override.aes = list(shape = shapes_present, size = 2.5))) +
    scale_x_continuous(limits = c(min_c, max_c), breaks = seq(min_c, max_c, by = 2))+
    scale_y_continuous(limits = c(min_rt, max_rt))

  if ("comments" %in% colnames(data_graph) &&
      any(!is.na(data_graph$comments) & nzchar(trimws(data_graph$comments)))) {
    graph <- graph +

      geom_point(
        data = subset(data_graph, !is.na(comments) & comments != ""),
        aes(C, RTsec),
        shape = 22,
        size = 5,
        fill = "black",
        alpha = 0.25
      )
  }

  if (!is.null(subtitle)) {
    graph <- graph +
      labs(subtitle = subtitle) +
      theme(plot.subtitle = element_text(color = "#5A8D78"))
  }


  return(graph)
}

# predictgraph
#' Plot the retention time models and predictions for the different lipid classes
#'
#' Plot the retention time models and predictions for the different lipid classes
#'
#' @param rtdata rtdata
#' @param classes vector of classes to plot.
#' @param eval_predict_plot logical. If TRUE, in the prediction plot include a
#' third graph that shows all the proposed candidates regardless of their
#' prediction score. If FALSE, do not add the new graph.
#'
#' @return rtdata
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
predictgraph <- function(rtdata, classes = NULL,
                         eval_predict_plot = TRUE) {

  data <- rtdata$filtering$data
  candidates <- rtdata$prediction$candidates
  predictions <- rtdata$prediction$predictions
  models <- rtdata$prediction$models

  if ("Filter_final" %in% colnames(candidates)) {
    data <- data[!is.na(data$lipidname) & data$lipidname != "" &
                   data$Filter_final %in% c("Kept",
                                            "Kept manually by the user",
                                            "Recovered with RT modelling") &
                   !is.na(data$RTsec), , drop = FALSE]

    classes_candidates <- unique(candidates$LipidSubclass[candidates$Filter_final %in% c("Kept", "Discarded",
                                                                                         "Kept manually by the user",
                                                                                         "Removed manually by the user",
                                                                                         "Discarded - In-source fragments",
                                                                                         "Discarded - RT coherence between subclasses",
                                                                                         "Discarded - Removed by RT rules",
                                                                                         "Discarded - Outlier",
                                                                                         "Discarded - Distant data")])

  } else {
    data <- data[!is.na(data$lipidname) & data$lipidname != "" &
                   data$checkfilter %in% c("Kept",
                                           "Kept manually by the user",
                                           "Recovered with RT modelling") &
                   !is.na(data$RTsec), , drop = FALSE]

    classes_candidates <- unique(candidates$LipidSubclass[candidates$checkcandidate %in% c("Kept", "Discarded",
                                                                                           "Kept manually by the user",
                                                                                           "Removed manually by the user")])
  }

  # All valid classes in the dataset
  classes_validas <- unique(c(na.omit(data$LipidSubclass[data$LipidSubclass != ""]),
                              classes_candidates))
  if (is.null(classes)) {
    classes <- classes_validas
  } else {
    classes <- classes[classes %in% classes_validas]
  }


  plots <- c()

  # Iterate over lipid classes in predictions
  for (lipidsubclass in classes) {
    C <- RT <- DB <- RTsec_rtpred <- legend_label <- Score_rtpred <- c()

    if ("Filter_final" %in% colnames(candidates)) {
      candidates$checkcandidate <- candidates$Filter_final
    }

    # Filter modeled data and candidates by lipid class
    data_model <- data[data$LipidSubclass == lipidsubclass, , drop = FALSE]
    data_candidates <- candidates[candidates$LipidSubclass == lipidsubclass &
                                    candidates$checkcandidate %in% c("Kept",
                                                                     "Discarded",
                                                                     "Discarded - Removed by RT rules",
                                                                     "Discarded - Outlier",
                                                                     "Discarded - In-source fragments",
                                                                     "Discarded - RT coherence between subclasses",
                                                                     "Discarded - Distant data",
                                                                     "Kept manually by the user",
                                                                     "Removed manually by the user"), , drop = FALSE]

    wanted_cols <- c("RTsec", "C", "DB", "shape", "fill", "legend_label", "comments")

    if (nrow(data_model) > 0 && !is.null(data_model)) {
      # Asegurar que data_model tenga la columna comments
      if (!"comments" %in% colnames(data_model)) {
        data_model$comments <- ""
      }

      data_model$shape <- 16
      data_model$fill <- "black"
      data_model$legend_label <- "Data Model"
    } else {
      data_model <- setNames(
        data.frame(matrix(ncol = length(wanted_cols), nrow = 0)),
        wanted_cols
      )
    }


    if (nrow(data_candidates) > 0) {

      # Asegurar que data_candidates tenga la columna comments
      if (!"comments" %in% colnames(data_candidates)) {
        data_candidates$comments <- ""
      }

      data_candidates$shape <- 22
      data_candidates$fill <- ifelse(
        data_candidates$checkcandidate %in% c("Kept", "Kept manually by the user"),
        "white", "grey"
      )
      data_candidates$legend_label <- ifelse(
        data_candidates$checkcandidate %in% c("Kept", "Kept manually by the user"),
        "Accepted Candidates",
        "Discarded Candidates"
      )

      data_class <- rbind(
        data_model[, wanted_cols, drop = FALSE],
        data_candidates[, wanted_cols, drop = FALSE]
      )

    } else {

      data_class <- data_model[, wanted_cols, drop = FALSE]

    }


    data_class$legend_label <- factor(data_class$legend_label,
                                      levels = c("Data Model",
                                                 "Accepted Candidates",
                                                 "Discarded Candidates"))

    # Global shapes and fills
    shapes <- c("Data Model" = 16,
                "Accepted Candidates" = 22,
                "Discarded Candidates" = 22)
    fills  <- c("Data Model" = "black",
                "Accepted Candidates" = "white",
                "Discarded Candidates" = "grey")


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
        subtitle_text <- paste("This model was built using the lipid subclass:", lipid_subclasses)
      } else {
        subtitle_text <- paste("This model was built using the following lipid subclasses:\n   -",
                               paste0(lipid_subclasses, collapse = ", "))
      }
    } else {
      subtitle_text <- "No model for this lipid class"
    }


    min_c <- min(unique(data_class$C)) - 1
    max_c <- max(unique(data_class$C)) + 1
    min_c <- ifelse(min_c %% 2 == 0, min_c - 2, min_c - 1)
    max_c <- ifelse(max_c %% 2 == 0, max_c + 2, max_c + 1)

    min_rt <- floor(min(data_class$RTsec)) - 50
    if (min_rt < 0) {
      min_rt <- 0
    }
    max_rt <- ceiling(max(data_class$RTsec)) + 50

    # check if prediction layer is required
    prediction_layer <- predictionlayergraph(rtdata = rtdata,
                                             lipidsubclass = lipidsubclass,
                                             min_c = min_c,
                                             max_c = max_c,
                                             min_rt = min_rt,
                                             max_rt = max_rt,
                                             alpha = 0.8,
                                             typeanalysis = "prediction")

    # Graph1
    data_graph1 <- data_class[data_class$legend_label %in% c("Data Model"), , drop = FALSE]
    graph1 <- createpredictgraph(
      data_graph = data_graph1,
      min_c = min_c,
      max_c = max_c,
      min_rt = min_rt,
      max_rt = max_rt,
      title = paste0("Modeling of the ", lipidsubclass," class"),
      subtitle = subtitle_text,
      prediction_layer = prediction_layer,
      shapes_present = shapes["Data Model"],
      fills_present = fills["Data Model"]
    )

    # Graph2
    graph2 <- createpredictgraph(
      data_graph = data_class,
      min_c = min_c,
      max_c = max_c,
      min_rt = min_rt,
      max_rt = max_rt,
      title = "Candidates based on \nm/z and preferred adducts",
      prediction_layer = prediction_layer,
      shapes_present = shapes[names(shapes) %in% unique(data_class$legend_label)],
      fills_present = fills[names(fills) %in% unique(data_class$legend_label)]
    )

    # Graph3
    data_graph3 <- data_class[data_class$legend_label %in%
                                c("Data Model", "Accepted Candidates"), , drop = FALSE]
    graph3 <- createpredictgraph(
      data_graph = data_graph3,
      min_c = min_c,
      max_c = max_c,
      min_rt = min_rt,
      max_rt = max_rt,
      title = "New proposed annotations \nbased on RT",
      prediction_layer = prediction_layer,
      shapes_present = shapes[names(shapes) %in% unique(data_graph3$legend_label)],
      fills_present = fills[names(fills) %in% unique(data_graph3$legend_label)]
    )

    # Generate and save graphs to a PDF file
    grDevices::pdf(NULL)
    grDevices::dev.control(displaylist = "enable")

    if (eval_predict_plot) {
      final_graph <- graph1 + theme(legend.position = "none") +
        graph2 +
        graph3 + theme(legend.position = "none") +
        patchwork::plot_layout(ncol = 3, nrow = 1, guides = "collect")
    } else {
      final_graph <- graph1 + graph3 +
        patchwork::plot_layout(ncol = 2, nrow = 1, guides = "collect")
    }
    print(final_graph)

    plots[[lipidsubclass]] <- grDevices::recordPlot() # Save plot in msbatch
    invisible(grDevices::dev.off()) # Close PDF file
  }

  rtdata$prediction$plots <- plots

  return(rtdata)
}

# groupidrtpred
#' Update rtdata$features with the retention time predictions of potential candidates
#'
#' Update rtdata$features  with the retention time predictions of potential candidates
#'
#' @param rtdata rtdata
#'
#' @return featuresPrediction
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
groupidrtpred <- function(rtdata) {

  features   <- rtdata$features
  candidates <- rtdata$prediction$candidates
  source     <- rtdata$metadata$source

  sep.lipid <- "|"

  ## ------------------------------
  ## Column configuration by source
  ## ------------------------------
  if (source == "lipidms") {
    position_cols <- c(
      "mz", "minmz", "maxmz", "RT", "RTsec", "minRT", "maxRT",
      "iniRT", "endRT", "RTminutes", "npeaks", "group",
      "isotope", "isoGroup", "LipidMSid", "Adduct",
      "confidenceLevel", "Score", "ScoreInt", "nsamples",
      "LipidRTid", "GroupID",
      "filter_label", "comments_filter", "lipidname_filter", "expandedlipidname_filter",
      "LipidClass_filter", "LipidSubclass_filter",
      "C_filter", "DB_filter", "adduct_filter",
      "confidenceLevel_filter", "Score_filter"
    )
  } else if (source == "msdial") {
    position_cols <- c(
      "Alignment ID", "Average Rt(min)", "Average Mz",
      "Metabolite name", "Adduct type",
      "LipidRTid", "GroupID",
      "filter_label", "comments_filter", "lipidname_filter", "expandedlipidname_filter",
      "LipidClass_filter", "LipidSubclass_filter",
      "C_filter", "DB_filter", "adduct_filter",
      "confidenceLevel_filter", "Score_filter"
    )
  } else {
    position_cols <- c(
      "mz", "RTsec", "lipidname", "adduct",
      "confidenceLevel", "LipidRTid", "GroupID",
      "filter_label", "comments_filter", "lipidname_filter", "expandedlipidname_filter",
      "LipidClass_filter", "LipidSubclass_filter",
      "C_filter", "DB_filter", "adduct_filter",
      "confidenceLevel_filter", "Score_filter"
    )
  }

  prediction_cols <- c(
    "lipidname_rtpred", "expandedlipidname_rtpred",
    "LipidClass_rtpred", "LipidSubclass_rtpred",
    "C_rtpred", "DB_rtpred", "adduct_rtpred",
    "confidenceLevel_rtpred", "RTsec_rtpred",
    "Score_rtpred", "comments_rtpred",
    "lipidname_final", "expandedlipidname_final",
    "LipidClass_final", "LipidSubclass_final",
    "C_final", "DB_final", "adduct_final",
    "confidenceLevel_final", "Score_final",
    "comments_final"
  )

  ## ------------------------------
  ## Aggregate prediction candidates
  ## ------------------------------
  if (!is.null(candidates) && nrow(candidates) > 0) {

    if ("Filter_final" %in% colnames(candidates)) {
      candidates$checkcandidate <- candidates$Filter_final
    }

    candidates <- candidates[
      candidates$checkcandidate %in% c("Kept", "Kept manually by the user"),
      , drop = FALSE
    ]

    if (nrow(candidates) > 0) {

      if (!"confidenceLevel" %in% colnames(candidates)) {
        candidates$confidenceLevel <- "Predicted"
      }
      candidates$confidenceLevel <- factor(
        candidates$confidenceLevel,
        levels = c("FA position", "FA", "Subclass", "MS-only", "MSMS", "Predicted",
                   "Other", "Unsupported class", "Not pass QC", ""),
        ordered = TRUE)

      candidates <- candidates[order(candidates$Score_rtpred, decreasing = TRUE),, drop = FALSE]

      split_candidates <- split(candidates, candidates$LipidRTid)

      featuresCandidates <- lapply(split_candidates, function(df) {

        ## ---- SAME LOGIC AS FILTERING ----
        same_core <-
          length(unique(df$LipidClass[df$LipidClass != ""])) == 1 &&
          length(unique(df$LipidSubclass[df$LipidSubclass != ""])) == 1 &&
          length(unique(df$C[df$C != ""])) == 1 &&
          length(unique(df$DB[df$DB != ""])) == 1 &&
          length(unique(df$adduct[df$adduct != ""])) == 1

        force_paste <- !same_core

        data.frame(
          LipidRTid = df$LipidRTid[1],
          GroupID_rtpred = collapse_conditional(df$GroupID, sep.lipid, force_paste),
          lipidname_rtpred = collapse_conditional(df$lipidname, sep.lipid, force_paste),
          expandedlipidname_rtpred = collapse_conditional(df$lipidname_expanded, sep.lipid, force_paste),
          LipidClass_rtpred = collapse_conditional(df$LipidClass, sep.lipid, force_paste),
          LipidSubclass_rtpred = collapse_conditional(df$LipidSubclass, sep.lipid, force_paste),
          C_rtpred = collapse_conditional(df$C, sep.lipid, force_paste),
          DB_rtpred = collapse_conditional(df$DB, sep.lipid, force_paste),
          adduct_rtpred = collapse_conditional(df$adduct, sep.lipid, force_paste),
          confidenceLevel_rtpred = collapse_conditional(df$confidenceLevel, sep.lipid, force_paste),
          RTsec_rtpred = collapse_conditional(df$RTsec_rtpred, sep.lipid, force_paste),
          Score_rtpred = collapse_conditional(df$Score_rtpred, sep.lipid, force_paste),
          comments_rtpred = collapse_conditional(df$comments, sep.lipid, force_paste),
          check.names = FALSE,
          stringsAsFactors = FALSE

        )
      })

      featuresCandidates <- do.call(rbind, featuresCandidates)

      ## Remove previous prediction columns
      features <- features[, !names(features) %in% prediction_cols, drop = FALSE]

      ## Merge predictions
      features <- merge(features, featuresCandidates,
                        by = "LipidRTid", all.x = TRUE)
    }
  }

  ## ------------------------------
  ## NA handling
  ## ------------------------------
  for (col in prediction_cols) {
    if (!col %in% names(features)) {
      features[[col]] <- ""
    } else {
      features[[col]][is.na(features[[col]]) | features[[col]] == "NA"] <- ""
    }
  }

  ## ------------------------------
  ## Build FINAL columns
  ## ------------------------------
  if (all(c("GroupID", "GroupID_rtpred") %in% names(features))) {

    gid  <- features$GroupID
    gidp <- features$GroupID_rtpred

    features$GroupID <- ifelse(
      gidp == "" | is.na(gidp),                   # si rtpred no existe
      gid,
      ifelse(
        gid == "" | is.na(gid),                    # si original no existe
        gidp,
        ifelse(
          gid == gidp,                             # si son iguales
          gid,
          paste(gid, gidp, sep = sep.lipid)        # si son distintos
        )
      )
    )

    features$GroupID_rtpred <- NULL
  }

  features$lipidname_final       <- paste0(features$lipidname_filter,
                                           features$lipidname_rtpred,
                                           sep = sep.lipid)
  features$expandedlipidname_final    <- paste0(features$expandedlipidname_filter,
                                                features$expandedlipidname_rtpred,
                                                sep = sep.lipid)
  features$LipidClass_final      <- paste0(features$LipidClass_filter,
                                           features$LipidClass_rtpred,
                                           sep = sep.lipid)
  features$LipidSubclass_final   <- paste0(features$LipidSubclass_filter,
                                           features$LipidSubclass_rtpred,
                                           sep = sep.lipid)
  features$C_final               <- paste0(features$C_filter,
                                           features$C_rtpred,
                                           sep = sep.lipid)
  features$DB_final              <- paste0(features$DB_filter,
                                           features$DB_rtpred,
                                           sep = sep.lipid)
  features$adduct_final          <- paste0(features$adduct_filter,
                                           features$adduct_rtpred,
                                           sep = sep.lipid)
  features$confidenceLevel_final <- paste0(features$confidenceLevel_filter,
                                           features$confidenceLevel_rtpred,
                                           sep = sep.lipid)
  features$Score_final           <- paste0(features$Score_filter,
                                           features$Score_rtpred,
                                           sep = sep.lipid)
  features$comments_final        <- paste0(features$comments_filter,
                                           features$comments_rtpred,
                                           sep = sep.lipid)

  ## Clean separators
  clean_sep <- paste0("^\\", sep.lipid, "|\\", sep.lipid, "$")
  final_cols <- grep("_final$", names(features), value = TRUE)
  features[final_cols] <- lapply(features[final_cols], gsub,
                                 pattern = clean_sep, replacement = "")

  ## Column order
  col_order <- c(
    position_cols,
    prediction_cols,
    setdiff(names(features), c(position_cols, prediction_cols))
  )

  rtdata$features <- features[, col_order[col_order %in% names(features)], drop = FALSE]

  # summarizedFeatures
  cols <- c("LipidRTid", "GroupID", "mz", "RTsec", "RTmin", "lipidname_final", "expandedlipidname_final",
            "LipidClass_final", "LipidSubclass_final", "C_final", "DB_final",
            "adduct_final", "confidenceLevel_final", "Score_final", "comments_final")

  # Encontrar la última columna cuyo nombre termina en "_final"
  final_cols <- which(grepl("_final$", names(features)))

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
    stop("Neither RTsec nor RT column found in features")
  }

  # RTmin
  if (!"RTmin" %in% colnames(summ) &
      "RTsec" %in% colnames(summ)) {
    summ$RTmin <- as.numeric(summ$RTsec) / 60
  }

  # RTsec
  if (!"RTsec" %in% colnames(summ) &
      "RTmin" %in% colnames(summ)) {
    summ$RTsec <- as.numeric(summ$RTmin) * 60
  }

  final_cols <- cols[cols %in% colnames(summ)]

  rtdata$summarizedFeatures <- summ[, final_cols, drop = FALSE]

  return(rtdata)
}




# get_prediction_summary
#' Get a summary of RT-based prediction results to the console.
#'
#' This function create a formatted summary of the RT prediction step,
#' reporting, for each lipid class, the number of newly accepted candidates
#' proposed by the RT model and the total number of final annotated IDs. Results
#' are printed directly to the console, with lipid class names highlighted
#' (if the console supports ANSI formatting).
#'
#' @param rtdata rtdata
#'
#' @return information about prediction
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
get_prediction_summary <- function(rtdata) {

  lines <- c("RT Prediction Summary:\n")

  for (lipid in unique(rtdata$filtering$information$LipidSubclass)) {

    kept <- rtdata$filtering$information$keep[
      rtdata$filtering$information$LipidSubclass == lipid
    ] +
      rtdata$filtering$information$recover[
        rtdata$filtering$information$LipidSubclass == lipid
      ]

    if (lipid %in% rtdata$prediction$information$LipidSubclass) {

      new_annotations <- rtdata$prediction$information$new_rtannotations[
        rtdata$prediction$information$LipidSubclass == lipid
      ]

      if (!is.na(new_annotations) &&
          new_annotations != "No proposed candidates") {

        total_ids <- new_annotations + kept

        lines <- c(
          lines,
          sprintf(
            "  - %-10s |  %4d accepted candidates |  %4d total final IDs",
            lipid, new_annotations, total_ids
          )
        )

      } else {

        lines <- c(
          lines,
          sprintf(
            "  - %-10s |    No proposed candidates |  %4d total final IDs",
            lipid, kept
          )
        )
      }

    } else {

      if (lipid %in% unique(rtdata$prediction$candidates$LipidSubclass)) {

        lines <- c(
          lines,
          sprintf(
            "  - %-10s |           model not built |  %4d total final IDs",
            lipid, kept
          )
        )

      } else {

        lines <- c(
          lines,
          sprintf(
            "  - %-10s |    No proposed candidates |  %4d total final IDs",
            lipid, kept
          )
        )
      }
    }
  }

  return(lines)
}


# printpredictioninfo
#' Print a summary of RT-based prediction results to the console.
#'
#' This function displays a formatted summary of the RT prediction step,
#' reporting, for each lipid class, the number of newly accepted candidates
#' proposed by the RT model and the total number of final annotated IDs. Results
#' are printed directly to the console, with lipid class names highlighted
#' (if the console supports ANSI formatting).
#'
#' @param rtdata rtdata
#'
#' @return information about prediction
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
printpredictioninfo <- function(rtdata) {
  cat(paste(get_prediction_summary(rtdata), collapse = "\n"))
  cat("\n")
}


# grouplipidclasses
#' Grouped lipid classes to easily separate them when building models.
#'
#' Grouped lipid classes to easily separate them when building models.
#'
#' @param data annotated peaklist, with each annotation separated into a row (filtering$data)
#'
#' @return db_mapping, a list with Classes and Adducts
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
grouplipidclasses <- function(data = NULL, return_map = FALSE, map_type = c("class", "model")) {

  # Ensure map_type is either "class" or "model"
  map_type <- match.arg(map_type)

  # Mapping for lipid classification:
  # Each lipid class is assigned to a single category (no overlaps, priority-based)
  class_map <- list(
    PX   = c("PC","PE","PS","PI","PG","PA"),
    LPX  = c("LPC","LPE","LPS","LPI","LPG","LPA"),
    `SM-Cer` = c("SM"),
    PCx  = c("PC-O","PC-P"),
    PEx  = c("PE-O","PE-P"),
    xCer = c("SHexCer_d","HexCer_d","Cer_d","AHexCer_d",
             "SHexCer_t","HexCer_t","Cer_t","AHexCer_t"),
    xDG  = c("DG","MGDG","DGDG")
  )

  # Mapping for modeling purposes:
  # Defines which lipid classes are grouped together when building joint models
  # Overlaps are allowed here (e.g., "PC" appears in both PX and PCx)
  model_map <- list(
    PX   = c("PC","PE","PS","PI","PG","PA"),
    LPX  = c("LPC","LPE","LPS","LPI","LPG","LPA"),
    `SM-Cer` = c("SM","Cer_d","Cer_t"),
    PCx  = c("PC","PC-O","PC-P"),
    PEx  = c("PE","PE-O","PE-P"),
    xCer = c("SHexCer_d","HexCer_d","Cer_d","AHexCer_d",
             "SHexCer_t","HexCer_t","Cer_t","AHexCer_t"),
    xDG  = c("DG","MGDG","DGDG")
  )

  # Select which map to use depending on the requested type
  selected_map <- if (map_type == "class") class_map else model_map

  # If no data is provided and the user only wants the map, return it
  if (is.null(data) && return_map) {
    return(selected_map)
  }

  # If no data is provided and no map is requested, issue a warning
  if (is.null(data)) {
    warning("You must provide 'data' or set return_map = TRUE")
  }

  # Only apply transformation when working with classification
  if (map_type == "class") {

    # Initialize the new column
    data$groupClass <- NA

    # Assign each lipid to a group based on class_map
    # The condition 'is.na' ensures that once assigned, it is not overwritten (priority matters)
    for (grp in names(class_map)) {
      idx <- data$LipidSubclass %in% class_map[[grp]] & is.na(data$groupClass)
      data$groupClass[idx] <- grp
    }

    # For lipid classes not found in any group, keep original class
    data$groupClass[is.na(data$groupClass)] <- data$LipidClass[is.na(data$groupClass)]
  }

  # Return modified dataset (or NULL if no data was provided)
  return(data)
}



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
#' @param distC_predict integer. Minimum difference in number of carbons from which lipids
#' are considered distant from the rest of lipids from the same class.
#' @param minSpecies_predict numeric. minimum number of lipids in a cluster from
#' which the cluster is considered correct from the rest.
#' @param plot logical. If TRUE is plotted. If FALSE, is not plotted
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
filterpredictionsdata <- function(rtdata, distC_predict, minSpecies_predict, plot = TRUE) {

  # save rtscorethr parameter in rtdata
  # rtdata$metadata$prediction$rtscorethr <- as.numeric(rtscorethr)

  if (!"Filter_final" %in% colnames(rtdata$filtering$data)) {
    rtdata$filtering$data$Filter_final <- rtdata$filtering$data$checkfilter
  }
  if (!"Filter_final" %in% colnames(rtdata$prediction$candidates)) {
    rtdata$prediction$candidates$Filter_final <- rtdata$prediction$candidates$checkcandidate
  }

  ##############################################################################
  # Check for in-source fragments to remove incorrect annotations
  rtdata <- checkinsourcefragments(rtdata = rtdata,
                                   typeanalysis = "prediction")

  ##############################################################################
  # Remove distant data points for RT modelling
  rtdata <- removedistantrtdatapoints(rtdata = rtdata,
                                      distC_filter = distC_predict,
                                      minSpecies_filter = minSpecies_predict,
                                      typeanalysis = "prediction")

  ##############################################################################
  # Check for RT coherence between different classes
  rtdata <- checkRTclasses(rtdata = rtdata,
                           typeanalysis = "prediction")

  ##############################################################################
  # Filtering of each lipid class according to its RT pattern based on its
  # composition of C and DB
  rtdata <- filterbyrtrules(rtdata = rtdata,
                            typeanalysis = "prediction")

  ##############################################################################
  # Check outliers in terms of number of carbons
  rtdata <- checkoutliers(rtdata = rtdata,
                          typeanalysis = "prediction")

  #============================================================================#
  # Plot filtering results
  #============================================================================#
  if (plot == TRUE) {
    rtdata <- filteringgraph(rtdata = rtdata)

    rtdata <- predictgraph(rtdata = rtdata,
                           classes = NULL,
                           eval_predict_plot = TRUE)
  }

  # if (is.null(rtdata$prediction$candidates$comments)) {
  #   rtdata$prediction$candidates$comments <- ""
  # }

  #============================================================================#
  # Organize and save results
  #============================================================================#
  ##############################################################################
  rtdata <- organizertdata(rtdata = rtdata,
                           typeanalysis = "filtering")
  rtdata <- organizertdata(rtdata = rtdata,
                           typeanalysis = "prediction")


  return(rtdata)
}


# filterpredictionsbyadducts
#' This function filters lipid predictions based on the level of adduct support
#' observed for each feature within the same lipid class but different molecular
#' compositions. When multiple candidate annotations are available, those supported
#' by more than one adduct are prioritised and retained.
#'
#' If all candidate annotations for a feature are supported by only a single adduct,
#' all annotations are retained, as no additional evidence is available to
#' discriminate between them.
#'
#' @param rtdata rtdata object.
#'
#' @return rtdata object with filtered lipid predictions based on adduct support.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' rtdata <- filterpredictionsbyadducts(rtdata)
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
filterpredictionsbyadducts <- function(rtdata) {

  data <- rtdata$prediction$candidates

  if (!"checkcandidate" %in% names(data)) {
    data$checkcandidate <- NA_character_
  }

  for (rtid in unique(data$LipidRTid)) {

    idx_rtid <- which(data$LipidRTid == rtid)

    if (length(idx_rtid) < 2) {
      next
    }

    subdf <- data[idx_rtid, , drop = FALSE]

    for (subclass in unique(subdf$LipidSubclass)) {

      idx_subclass_local <- which(subdf$LipidSubclass == subclass)

      if (length(idx_subclass_local) < 2) {
        next
      }

      sub_subdf <- subdf[idx_subclass_local, , drop = FALSE]

      # ¿Hay distintas composiciones?
      comp <- paste(sub_subdf$C, sub_subdf$DB, sep = "_")

      if (length(unique(comp)) < 2) {
        next
      }

      n_adducts <- vapply(
        sub_subdf$adducts_expanded,
        count_adductsexpanded,
        integer(1)
      )

      has_multi <- any(n_adducts > 1)
      has_single <- any(n_adducts <= 1)

      if (has_multi && has_single) {

        idx_remove_local <- idx_subclass_local[n_adducts <= 1]
        idx_remove_global <- idx_rtid[idx_remove_local]

        # if (data$checkcandidates[idx_remove_global] %in% c("Kept", "Kept manually by the user")) {
        if (data$checkcandidate[idx_remove_global] %in% c("Kept")) {
          # cat(paste0("idx:", rtid, "\n"))
          data$checkcandidate[idx_remove_global] <- "Discarded"
        }

      }
    }
  }

  rtdata$prediction$candidates <- data

  return(rtdata)
}
# filterpredictionsbyadducts <- function(rtdata) {
#
#   data <- rtdata$prediction$candidates
#
#   if (!"checkcandidate" %in% names(data)) {
#     data$checkcandidate <- NA_character_
#   }
#
#   for (rtid in unique(data$LipidRTid)) {
#
#     idx_rtid <- which(data$LipidRTid == rtid)
#
#     if (length(idx_rtid) < 2) next
#
#     subdf <- data[idx_rtid, , drop = FALSE]
#
#     for (subclass in unique(subdf$LipidSubclass)) {
#
#       idx_subclass_local <- which(subdf$LipidSubclass == subclass)
#
#       if (length(idx_subclass_local) < 2) next
#
#       sub_subdf <- subdf[idx_subclass_local, , drop = FALSE]
#
#       comp <- paste(sub_subdf$C, sub_subdf$DB, sep = "_")
#
#       if (length(unique(comp)) < 2) next
#
#       n_adducts <- vapply(
#         sub_subdf$adducts_expanded,
#         count_adductsexpanded,
#         integer(1)
#       )
#
#       has_multi <- any(n_adducts > 1)
#       has_single <- any(n_adducts <= 1)
#
#       if (has_multi && has_single) {
#
#         idx_remove_local <- idx_subclass_local[n_adducts <= 1]
#         idx_remove_global <- idx_rtid[idx_remove_local]
#
#         data$checkcandidate[idx_remove_global] <- "Discarded"
#       }
#     }
#   }
#
#   rtdata$prediction$candidates <- data
#   return(rtdata)
# }


# count_adductsexpanded
#' Count unique adducts in a string of adduct annotations
#'
#' This function counts the number of unique adducts present in a given string
#' of adduct annotations. Adducts are expected to be separated by the "|" character.
#' The function removes empty values and NA entries before computing the number
#' of unique adduct types.
#'
#' @param x character string containing adduct annotations separated by "|".
#'
#' @return integer. Number of unique adducts found in the input string.
#'
#' @keywords internal
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
count_adductsexpanded <- function(x) {
  if (is.na(x) || x == "") {
    return(0L)
  }

  adducts <- unlist(strsplit(x, "\\|"))
  adducts <- trimws(adducts)
  adducts <- adducts[!is.na(adducts) & adducts != ""]

  length(unique(adducts))
}
