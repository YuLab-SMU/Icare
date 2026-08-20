
#' Universal Prediction Dispatcher for Deployment
#'
#' Internal function that handles model selection and prediction for any
#' model type (base, tuned, or ensemble) stored in the `Train_Model` object.
#'
#' @param object The `Train_Model` S4 object.
#' @param newdata A data frame of new observations.
#' @param preproc A `preProcess` object or `NULL`.
#' @param class_labels Character vector of length 2.
#' @param selected_model Character identifier. Special values: `"auto"`
#'   chooses the best available (tuned > first ensemble > first base model);
#'   `"Ensemble: Strategy"` selects the corresponding ensemble from
#'   `process.info$ensembles`; otherwise matches a model name.
#' @param tuned_label The exact label identifying the fine-tuned model.
#'
#' @return A data frame with two columns (negative, positive) of probabilities.
#' @keywords internal
.predict_deployment <- function(object, newdata, preproc = NULL,
                                class_labels = NULL, selected_model = "auto",
                                tuned_label = NULL) {
  
  # ---- 1. Preprocess new data ----
  newdata <- as.data.frame(newdata)
  if (!is.null(preproc) && inherits(preproc, "preProcess")) {
    newdata <- tryCatch(
      predict(preproc, newdata),
      error = function(e) {
        stop("Preprocessing failed on newdata: ", e$message, call. = FALSE)
      }
    )
  }
  
  # ---- 2. Resolve selected_model ----
  res <- object@best.model.result
  all_models <- object@train.models
  ensembles <- object@process.info$ensembles
  
  if (selected_model == "auto") {
    if (!is.null(tuned_label) && !is.null(res$fine_tuned_model)) {
      selected_model <- tuned_label
    } else if (is.list(ensembles) && length(ensembles) > 0) {
      selected_model <- paste0("Ensemble: ", names(ensembles)[1])
    } else if (!is.null(res$ensemble)) {
      selected_model <- "Ensemble Stacking"
    } else if (length(all_models) > 0) {
      selected_model <- names(all_models)[1]
    } else {
      stop("No model available for prediction.")
    }
  }
  
  # ---- 3. Dispatch ----
  prob_matrix <- NULL
  
  # 3a. Ensemble (multiple strategies)
  if (grepl("^Ensemble: ", selected_model)) {
    strategy <- sub("^Ensemble: ", "", selected_model)
    if (is.list(ensembles) && strategy %in% names(ensembles)) {
      ens_obj <- ensembles[[strategy]]
      predict_fn <- NULL
      if (is(ens_obj, "Train_Model")) {
        predict_fn <- ens_obj@best.model.result$ensemble$predict_fn
      } else if (is.list(ens_obj)) {
        predict_fn <- ens_obj$predict_fn
      }
      if (is.function(predict_fn)) {
        pred <- tryCatch(
          predict_fn(newdata),
          error = function(e) {
            stop("Ensemble '", strategy, "' prediction failed: ", e$message, call. = FALSE)
          }
        )
        prob_matrix <- .ensure_probability_matrix(pred, strategy)
      } else {
        stop("Ensemble '", strategy, "' does not have a valid predict_fn.", call. = FALSE)
      }
    } else {
      warning("Ensemble strategy '", strategy, "' not found. Falling back to 'auto'.",
              call. = FALSE)
      return(.predict_deployment(object, newdata, preproc, class_labels, "auto"))
    }
  }
  
  # 3a-legacy. Single legacy stacking ensemble
  if (is.null(prob_matrix) && identical(selected_model, "Ensemble Stacking") &&
      is.function(object@best.model.result$ensemble$predict_fn)) {
    pred <- tryCatch(
      object@best.model.result$ensemble$predict_fn(newdata),
      error = function(e) {
        stop("Legacy stacking ensemble prediction failed: ", e$message, call. = FALSE)
      }
    )
    prob_matrix <- .ensure_probability_matrix(pred, "Ensemble Stacking")
  }
  
  # 3b. Fine tuned model
  if (is.null(prob_matrix) && !is.null(tuned_label) &&
      identical(selected_model, tuned_label) && !is.null(res$fine_tuned_model)) {
    prob_matrix <- tryCatch(
      as.data.frame(predict(res$fine_tuned_model, newdata, type = "prob")),
      error = function(e) {
        stop("Fine tuned model prediction failed: ", e$message, call. = FALSE)
      }
    )
  }
  
  # 3c. Base models
  if (is.null(prob_matrix) && selected_model %in% names(all_models)) {
    prob_matrix <- tryCatch(
      as.data.frame(predict(all_models[[selected_model]], newdata, type = "prob")),
      error = function(e) {
        stop("Model '", selected_model, "' prediction failed: ", e$message, call. = FALSE)
      }
    )
  }
  
  # 3d. Fallback
  if (is.null(prob_matrix)) {
    if (!is.null(res$fine_tuned_model)) {
      prob_matrix <- as.data.frame(predict(res$fine_tuned_model, newdata, type = "prob"))
    } else if (is.list(ensembles) && length(ensembles) > 0) {
      first_ens <- ensembles[[1]]
      predict_fn <- NULL
      if (is(first_ens, "Train_Model")) {
        predict_fn <- first_ens@best.model.result$ensemble$predict_fn
      } else if (is.list(first_ens)) {
        predict_fn <- first_ens$predict_fn
      }
      if (is.function(predict_fn)) {
        pred <- predict_fn(newdata)
        prob_matrix <- .ensure_probability_matrix(pred, names(ensembles)[1])
      }
    } else if (length(all_models) > 0) {
      prob_matrix <- as.data.frame(predict(all_models[[1]], newdata, type = "prob"))
    } else {
      stop("No fallback model could be used for prediction.", call. = FALSE)
    }
  }
  
  # ---- 4. Ensure correct column names and order ----
  if (is.null(prob_matrix)) {
    stop("Prediction matrix could not be generated.", call. = FALSE)
  }
  
  if (ncol(prob_matrix) == 1) {
    probs_pos <- prob_matrix[, 1]
    probs_neg <- 1 - probs_pos
    prob_matrix <- data.frame(Neg = probs_neg, Pos = probs_pos,
                              stringsAsFactors = FALSE)
  } else if (ncol(prob_matrix) == 2) {
    orig_cols <- colnames(prob_matrix)
    if (!is.null(orig_cols) && length(orig_cols) == 2) {
      pos_idx <- grep("pos|positive|case|1|high|risk", orig_cols, ignore.case = TRUE)
      if (length(pos_idx) == 1) {
        neg_idx <- setdiff(1:2, pos_idx)
        prob_matrix <- prob_matrix[, c(neg_idx, pos_idx), drop = FALSE]
      } else {
        warning("Cannot determine positive class column. Assuming first column = negative, second = positive.",
                call. = FALSE)
      }
    }
    colnames(prob_matrix) <- c("Neg", "Pos")
    if (!is.null(class_labels) && length(class_labels) == 2) {
      colnames(prob_matrix) <- class_labels
    }
  } else {
    stop("Prediction returned ", ncol(prob_matrix),
         " columns; expected 1 or 2.", call. = FALSE)
  }
  
  return(prob_matrix)
}


#' Helper: Convert prediction output to a two column probability matrix
#'
#' @param pred Output from an ensemble predict_fn.
#' @param strategy Name of the ensemble (for error messages).
#' @return A data frame with two columns (negative, positive) of probabilities.
#' @keywords internal
.ensure_probability_matrix <- function(pred, strategy) {
  if (is.numeric(pred)) {
    # If it's a vector, assume it's the positive class probability
    if (is.vector(pred) || (is.matrix(pred) && ncol(pred) == 1)) {
      pos <- as.numeric(pred)
      neg <- 1 - pos
      return(data.frame(Neg = neg, Pos = pos, stringsAsFactors = FALSE))
    } else if (is.matrix(pred) && ncol(pred) == 2) {
      return(as.data.frame(pred))
    } else {
      stop("Ensemble '", strategy, "' returned a numeric object with ",
           if (is.matrix(pred)) ncol(pred) else length(pred),
           " values; cannot determine class probabilities.", call. = FALSE)
    }
  } else if (is.data.frame(pred)) {
    if (ncol(pred) == 1) {
      pos <- pred[, 1]
      neg <- 1 - pos
      return(data.frame(Neg = neg, Pos = pos, stringsAsFactors = FALSE))
    } else if (ncol(pred) == 2) {
      return(pred)
    } else {
      stop("Ensemble '", strategy, "' returned a data frame with ", ncol(pred),
           " columns; expected 1 or 2.", call. = FALSE)
    }
  } else if (is.factor(pred) || is.character(pred)) {
    # For voting or categorical output: convert to probability as fraction of positive class
    # This is a heuristic; better to ensure ensemble returns probabilities.
    warning("Ensemble '", strategy, "' returned class labels, not probabilities. ",
            "Converting to binary probability (0/1) based on the positive class.",
            call. = FALSE)
    # NOTE: unique(pred) returns values in order of FIRST APPEARANCE in the
    # data, which is essentially arbitrary and NOT a reliable way to decide
    # which class is "positive" -- a previous version used uniq[2] from
    # unique(pred), which could silently flip which class counts as
    # positive depending on row order. If `pred` is a factor, its levels
    # encode the package-wide convention (2nd level = positive, see
    # ModelTrainAnalysis()); for plain character input there is no such
    # ordering information at all, so we fall back to alphabetical order
    # and warn loudly that this may not match the intended positive class.
    if (is.factor(pred)) {
      uniq <- levels(pred)
      if (length(uniq) != 2) uniq <- unique(as.character(pred))
    } else {
      uniq <- sort(unique(pred))
      warning("Ensemble '", strategy, "' returned character labels (not a ",
              "factor), so the positive class cannot be determined from ",
              "level order and was inferred alphabetically as '", 
              if (length(uniq) == 2) uniq[2] else NA, "'. Verify this is ",
              "correct, or have the ensemble return a factor/probabilities ",
              "instead.", call. = FALSE)
    }
    if (length(uniq) == 2) {
      pos_class <- uniq[2]
      pos <- as.numeric(pred == pos_class)
      neg <- 1 - pos
      return(data.frame(Neg = neg, Pos = pos, stringsAsFactors = FALSE))
    } else {
      stop("Ensemble '", strategy, "' returned non binary categorical output. ",
           "Cannot convert to probabilities.", call. = FALSE)
    }
  } else {
    stop("Ensemble '", strategy, "' returned an unsupported type: ",
         paste(class(pred), collapse = "/"), call. = FALSE)
  }
}


#' Model Deployment Constructor with Multi‑Ensemble Support
#'
#' Creates a deployment-ready object that encapsulates all trained models,
#' fine tuned models, and all stored ensembles (from `process.info$ensembles`).
#' The returned list includes a unified prediction function that can dispatch
#' to any of these models based on a user selected string.
#'
#' @param object A `Train_Model` S4 object, typically returned by
#'   `ModelTrainAnalysis()` and subsequent `TrainEnsemble()` calls.
#' @param preproc A `caret::preProcess` object (or `NULL`) to apply to new data
#'   before prediction. Must be compatible with `predict(preproc, newdata)`.
#' @param class_labels Character vector of length 2 giving the names of the
#'   negative and positive classes, respectively. Default: `c("Normal", "High Risk")`.
#' @param model_description A short descriptive string about the model platform,
#'   displayed in the Shiny app's metadata card.
#'
#' @return A list of class `"ModelDeployment"` with components:
#'   \describe{
#'     \item{object}{The original `Train_Model` object.}
#'     \item{preproc}{The preprocessing object.}
#'     \item{ref_data}{The training data (used for reference/UI).}
#'     \item{group_col}{The name of the outcome column.}
#'     \item{class_labels}{The class labels.}
#'     \item{model_desc}{The description string.}
#'     \item{model_list}{Character vector of all selectable model identifiers.}
#'     \item{predict_fn}{A function `function(newdata, selected_model = "auto")`
#'       that returns a two column data frame of class probabilities.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' deploy <- ModelDeployment(
#'   object = model_obj,
#'   preproc = preProc,
#'   class_labels = c("Control", "Case"),
#'   model_description = "Omics based diagnostic ensemble"
#' )
#' pred <- deploy$predict_fn(newdata, selected_model = "Ensemble: Stacking")
#' }
ModelDeployment <- function(object,
                            preproc = NULL,
                            class_labels = c("Normal", "High Risk"),
                            model_description = "Clinlabomics based diagnostic ensemble.") {
  
  # ----- 1. Build the list of available model identifiers -----
  available <- c()
  
  # 1a. Base caret models (trained via ModelTrainAnalysis)
  base_models <- names(object@train.models)
  if (length(base_models) > 0) {
    available <- c(available, base_models)
  }
  
  # 1b. Fine tuned model (if present)
  tuned_label <- NULL
  if (!is.null(object@best.model.result$fine_tuned_model)) {
    tuned_method <- object@best.model.result$model_type %||% "tuned"
    tuned_label <- paste0(tuned_method, "_tuned")
    counter <- 0L
    while (tuned_label %in% available) {
      counter <- counter + 1L
      tuned_label <- paste0(tuned_method, "_tuned_", counter)
    }
    available <- c(tuned_label, available)
  }
  
  # 1c. Ensembles stored in process.info$ensembles (multiple strategies)
  ensembles <- object@process.info$ensembles
  if (is.list(ensembles) && length(ensembles) > 0) {
    valid_ensembles <- names(ensembles)[vapply(ensembles, function(e) {
      if (is(e, "Train_Model")) {
        !is.null(e@best.model.result$ensemble) &&
          is.function(e@best.model.result$ensemble$predict_fn)
      } else if (is.list(e)) {
        is.function(e$predict_fn)
      } else {
        FALSE
      }
    }, logical(1))]
    
    if (length(valid_ensembles) > 0) {
      ensemble_labels <- paste0("Ensemble: ", valid_ensembles)
      available <- c(ensemble_labels, available)
    } else {
      warning("No valid ensembles found in process.info$ensembles.")
    }
  } else {
    # Fallback: single legacy ensemble
    if (!is.null(object@best.model.result$ensemble) &&
        is.function(object@best.model.result$ensemble$predict_fn)) {
      available <- c("Ensemble Stacking", available)
    }
  }
  
  if (length(available) == 0) {
    stop("No trainable models or ensembles found.")
  }
  
  # Reference data (with fallbacks)
  ref_data <- object@split.data$training
  if (is.null(ref_data) || nrow(ref_data) == 0) ref_data <- object@filtered.set$training
  if (is.null(ref_data) || nrow(ref_data) == 0) ref_data <- object@clean.df
  if (is.null(ref_data) || nrow(ref_data) == 0) {
    stop("No reference data available. Please populate split.data$training, filtered.set$training, or clean.df.")
  }
  ref_data <- as.data.frame(ref_data)
  
  deployment <- list(
    object       = object,
    preproc      = preproc,
    ref_data     = ref_data,
    group_col    = object@group_col,
    class_labels = class_labels,
    model_desc   = model_description,
    model_list   = available,
    predict_fn   = function(newdata, selected_model = "auto") {
      .predict_deployment(object, newdata, preproc, class_labels, selected_model, tuned_label)
    }
  )
  class(deployment) <- "ModelDeployment"
  return(deployment)
}

#' @title Clinlabomics Terminal: Light-Tech Edition
#' @description Futuristic UI with high-contrast elements on a light grey base for clinical model deployment.
#' 
#' @param deployment A ModelDeployment object.
#' @param title Character string for the application window title.
#'
#' @import shiny bslib
#' @export
deploy_clinlab_app <- function(deployment, title = "Clinlabomics Terminal") {
  
  train_data <- deployment$ref_data
  feat_cols <- setdiff(colnames(train_data), c(".outcome", "group", "Group", deployment$group_col))
  
  # Base Tech Theme Configuration
  tech_theme <- bslib::bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#00dfc0", # Electric Cyan
    bg = "#f7f7f7",      # Light Grey Background
    fg = "#1a1d23",      # Dark Text
    base_font = bslib::font_google("Inter")
  )
  
  ui <- bslib::page_navbar(
    title = title, theme = tech_theme, bg = "#1a1d23", # Dark Tech Header
    
    # Custom CSS for Cyber-Clinical Aesthetics
    header = tags$style(HTML("
      body { background-color: #f7f7f7 !important; }
      .sidebar { background-color: #1a1d23 !important; color: white !important; }
      .sidebar h6 { color: #00dfc0 !important; font-weight: 800; letter-spacing: 1px; }
      .card { border-radius: 12px; border: 1px solid #e0e0e0; box-shadow: 0 4px 15px rgba(0,0,0,0.05); }
      .card-header { background-color: #ffffff !important; border-bottom: 1px solid #eee; font-weight: 700; color: #1a1d23; }
      .btn-primary { background-color: #00dfc0 !important; border: none; color: #1a1d23 !important; font-weight: 800; }
      .irs-bar { background: #00dfc0 !important; border-top: 1px solid #00dfc0 !important; border-bottom: 1px solid #00dfc0 !important; }
      .irs-from, .irs-to, .irs-single { background: #00dfc0 !important; color: #1a1d23 !important; }
    ")),
    
    bslib::nav_panel("Diagnostic Terminal",
                     bslib::layout_sidebar(
                       sidebar = bslib::sidebar(
                         title = "CORE CONFIG", width = 350,
                         h6("ALGORITHM SELECT"),
                         selectInput("model_choice", NULL, choices = deployment$model_list),
                         
                         h6("SENSITIVITY THRESHOLD"),
                         sliderInput("threshold", NULL, min = 0, max = 1, value = 0.5, step = 0.01),
                         
                         hr(style = "border-top: 1px solid #444;"),
                         h6("BIOMARKER INPUTS"),
                         lapply(feat_cols, function(f) {
                           val <- train_data[[f]]
                           if (is.numeric(val)) {
                             numericInput(paste0("in_", f), f, value = signif(stats::median(val, na.rm = TRUE), 4))
                           } else {
                             selectInput(paste0("in_", f), f, choices = levels(as.factor(val)))
                           }
                         }),
                         actionButton("go", "INITIATE ANALYSIS", class = "btn-primary w-100"),
                         downloadButton("download_json", "GENERATE DATA PACK", class = "btn-outline-light w-100 mt-2")
                       ),
                       
                       bslib::layout_column_wrap(
                         width = 1,
                         bslib::card(
                           bslib::card_header("ANALYTICAL ENGINE STATUS"),
                           bslib::layout_column_wrap(
                             width = 1/2,
                             uiOutput("res_ui"),
                             plotly::plotlyOutput("plot", height = "300px")
                           )
                         ),
                         bslib::card(
                           bslib::card_header("METADATA & SPECIFICATIONS"),
                           p(deployment$model_desc, style = "font-style: italic;"),
                           markdown(paste0(
                             "**Target:** ", deployment$group_col, "  \n",
                             "**System:** Clinlabomics-X v2.1  \n",
                             "**Status:** Synchronized with Training Set"
                           ))
                         )
                       )
                     )
    )
  )
  
  server <- function(input, output, session) {
    report_data <- eventReactive(input$go, {
      
      input_list <- lapply(feat_cols, function(f) {
        raw_val  <- input[[paste0("in_", f)]]
        orig_col <- train_data[[f]]
        
        if (is.numeric(orig_col)) {
          return(as.numeric(raw_val))
        } else if (is.factor(orig_col)) {
          return(factor(raw_val, levels = levels(orig_col)))
        } else {
          return(as.character(raw_val))
        }
      })
      
      names(input_list) <- feat_cols
      df_input <- as.data.frame(input_list, stringsAsFactors = FALSE)
      
      probs <- deployment$predict_fn(df_input, selected_model = input$model_choice)
      
      list(
        inputs        = as.list(df_input[1, ]), 
        probabilities = as.list(probs[1, ]), 
        threshold     = input$threshold, 
        model_used    = input$model_choice, 
        time          = Sys.time()
      )
    })
    
    output$res_ui <- renderUI({
      req(res <- report_data())
      pos_prob <- res$probabilities[[2]]
      is_high  <- pos_prob >= res$threshold
      label    <- if (is_high) names(res$probabilities)[2] else names(res$probabilities)[1]
      
      text_color <- if (is_high) "#ff2d55" else "#00dfc0" 
      
      div(style = "text-align:center; padding: 40px;",
          h4("PREDICTION RESULT", style = "color: #888; font-size: 0.9rem; letter-spacing: 2px;"),
          h1(label, style = paste0("color: ", text_color, "; font-weight: 900; font-size: 3.5rem; text-shadow: 0 0 10px ", text_color, "44;")),
          h5(sprintf("PREDICTED RISK: %.2f%%", pos_prob * 100), style = "color: #555; margin-top: 10px;")
      )
    })
    
    output$plot <- plotly::renderPlotly({
      req(res <- report_data())
      p_vals <- unlist(res$probabilities)
      plotly::plot_ly(x = names(p_vals), y = as.numeric(p_vals), type = "bar", 
                      marker = list(color = c('#e0e0e0', '#00dfc0'))) %>%
        plotly::layout(
          yaxis = list(title = "Probability", range = c(0, 1), gridcolor = "#f0f0f0"),
          xaxis = list(title = ""),
          paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor = 'rgba(0,0,0,0)'
        )
    })
    
    output$download_json <- downloadHandler(
      filename = function() { paste0("Clinlab_Export_", format(Sys.time(), "%Y%m%d_%H%M"), ".json") },
      content  = function(file) { writeLines(jsonlite::toJSON(report_data(), auto_unbox = TRUE, pretty = TRUE), file) }
    )
  }
  
  shinyApp(ui, server)
}