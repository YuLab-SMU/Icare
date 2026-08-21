# =============================================================================
# model_hyperparameter_tuning.R
# Flexible Hyperparameter Tuning Module
# =============================================================================

# 0. Package Check -----------------------------------------------------------
.check_tune_pkgs <- function() {
  required <- c("caret", "rBayesianOptimization", "dplyr", "ggplot2", "ggprism")
  missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Install them before using this module.")
  }
  invisible(TRUE)
}


# =============================================================================
# Step 1: Inspect tunable parameters for any caret model
# =============================================================================

#' Inspect Tunable Hyperparameters for a caret Model
#'
#' Returns a data frame listing all tuneable parameters and their default
#' ranges (when available), so the user can review and modify them before
#' passing to \code{BuildTuningBounds}.
#'
#' @param method Character. A caret model name (e.g., "rf", "xgbTree").
#' @return A data frame with columns: \code{parameter}, \code{label}, 
#'   \code{class}, and \code{default_range} (a character hint).
#' @export
#' @examples
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   InspectHyperParams("rf")
#' }
InspectHyperParams <- function(method) {
  .require_pkgs("rBayesianOptimization")
  .check_tune_pkgs()
  
  info <- caret::getModelInfo(method, regex = FALSE)[[1]]
  if (is.null(info)) stop("Model '", method, "' not found in caret.")
  params <- info$parameters
  # Add a suggestive default range based on parameter type and common practice
  params$default_range <- sapply(params$parameter, function(p) {
    switch(p,
           mtry               = "[1, floor(sqrt(n_features))]",
           n.trees            = "[50, 500]",
           nrounds            = "[50, 300]",
           max_depth          = "[2, 10]",
           interaction.depth  = "[1, 9]",
           eta                = "[0.01, 0.3]",
           shrinkage          = "[0.001, 0.1]",
           gamma              = "[0, 5]",
           colsample_bytree   = "[0.4, 1]",
           min_child_weight   = "[1, 10]",
           subsample          = "[0.5, 1]",
           sigma              = "[0.001, 0.1]",
           C                  = "[0.1, 10]",
           alpha              = "[0, 1]",
           lambda             = "[0.0001, 0.5]",
           cp                 = "[0.0001, 0.01]",
           degree             = "[1, 3]",
           nprune             = "[2, 20]",
           paste0("[0, 1]  # please set manually")
    )
  })
  
  cat("\n========================================\n")
  cat("  Hyperparameters for:", method, "\n")
  cat("========================================\n")
  print(params[, c("parameter", "label", "class", "default_range")], row.names = FALSE)
  cat("\nUse these parameters to build bounds with BuildTuningBounds().\n")
  cat("You can modify the default ranges as needed.\n\n")
  
  invisible(params)
}


# =============================================================================
# Step 2: Build a user-defined bounds list from inspected parameters
# =============================================================================

#' Build Bounds List for Bayesian Optimization
#'
#' Converts a user-provided named list of \code{c(lower, upper)} vectors into
#' the format required by \code{rBayesianOptimization::BayesianOptimization}.
#'
#' @param ... Named arguments, each a numeric vector of length 2 
#'   \code{c(lower, upper)}.
#' @return A named list of bounds.
#' @export
#'
#' @examples
#' bounds <- BuildTuningBounds(mtry = c(1, 5), n.trees = c(50, 200))
#' print(bounds)
BuildTuningBounds <- function(...) {
  bounds <- list(...)
  
  # Validate
  for (nm in names(bounds)) {
    if (length(bounds[[nm]]) != 2 || !is.numeric(bounds[[nm]]))
      stop("Each parameter must be a numeric vector of length 2: c(lower, upper). ",
           "Problem with: ", nm)
    if (bounds[[nm]][1] >= bounds[[nm]][2])
      stop("Lower bound must be < upper bound for: ", nm)
  }
  
  cat("Built bounds for", length(bounds), "parameters:", 
      paste(names(bounds), collapse = ", "), "\n")
  
  return(bounds)
}


# =============================================================================
# Step 3: Run Bayesian Optimization with user-defined bounds
# =============================================================================
#' Bayesian Optimization for Model Fine-Tuning
#'
#' @param model_obj An S4 object of class 'Train_Model'.
#' @param method Caret model method string.
#' @param bounds Named list defining parameter search boundaries.
#' @param init_points Initial random exploration points. Default 15.
#' @param n_iter Iteration steps for Bayesian Optimization. Default 30.
#' @param cv_folds Cross-validation folds. Default 5.
#' @param metric Optimization metric. Default "ROC".
#' @param summaryFun Performance evaluation summary function.
#' @param use_scaled Logical; whether to use scaled training dataset.
#' @param sampling Sampling strategy for class imbalance ("smote", "up", "down").
#' @param class_weights Logical; whether to apply class frequency weighting.
#' @param seed Random seed for reproducibility. Default 123.
#' @param verbose Logical; whether to print progress logs.
#'
#' @return Updated \code{model_obj} containing the fine-tuned model and optimization logs.
#' @export
FineTuneModel <- function(model_obj,
                          method,
                          bounds,
                          init_points    = 15,
                          n_iter          = 30,
                          cv_folds        = 5,
                          metric          = "ROC",
                          summaryFun      = caret::twoClassSummary,
                          use_scaled      = FALSE,
                          sampling        = NULL,
                          class_weights   = FALSE,
                          seed            = 123,
                          verbose         = TRUE) {
  .require_pkgs("rBayesianOptimization")
  
  .check_tune_pkgs()
  set.seed(seed)
  if (!inherits(model_obj, "Train_Model")) stop("model_obj must be a 'Train_Model' object.")
  
  # Step 1: Extract dataset
  train_data <- if (use_scaled) model_obj@split.scale.data$training else model_obj@filtered.set$training
  if (is.null(train_data)) train_data <- model_obj@split.data$training
  if (is.null(train_data)) stop("Training dataset is empty.")
  
  # Step 2: Align factor levels with make.names
  gc <- model_obj@group_col
  if (!is.factor(train_data[[gc]])) train_data[[gc]] <- as.factor(train_data[[gc]])
  levels(train_data[[gc]]) <- make.names(levels(train_data[[gc]]))
  
  n_features <- ncol(train_data) - 1
  
  # Step 3: Adjust mtry bounds for all tree-based ensemble methods
  rf_methods <- c("rf", "ranger", "Rborist", "randomForest")
  if (method %in% rf_methods && "mtry" %in% names(bounds)) {
    bounds$mtry[1] <- max(1, bounds$mtry[1])
    bounds$mtry[2] <- min(bounds$mtry[2], n_features)
    if (verbose) cat(">>> Adjusting mtry range to: [", bounds$mtry[1], ",", bounds$mtry[2], "]\n")
  }
  
  # Step 4: Compute class weights if requested
  wts <- NULL
  if (class_weights) {
    tab <- table(train_data[[gc]])
    wts <- as.numeric(1 / tab[as.character(train_data[[gc]])])
  }
  
  # Step 5: Identify integer parameters
  model_info <- tryCatch(
    caret::getModelInfo(method, regex = FALSE)[[1]],
    error = function(e) stop("Model method '", method, "' not found in caret.")
  )
  param_df <- model_info$parameters
  int_params <- param_df$parameter[param_df$class == "integer"]
  extra_int  <- c("n.trees", "nrounds", "n.minobsinnode", "min_child_weight", "mtry", "max_depth")
  int_params <- unique(c(int_params, intersect(extra_int, names(bounds))))
  
  # Step 6: Objective Function with Chance-Level Fallback Penalty
  obj_func <- function(...) {
    params <- list(...)
    for (p in intersect(names(params), int_params)) {
      params[[p]] <- round(params[[p]])
    }
    tune_grid <- do.call(data.frame, params)
    
    ctrl <- caret::trainControl(
      method          = "cv",
      number          = cv_folds,
      classProbs      = TRUE,
      summaryFunction = summaryFun,
      sampling        = sampling,
      verboseIter     = FALSE
    )
    
    res <- tryCatch({
      mod <- caret::train(
        as.formula(paste(gc, "~ .")),
        data      = train_data,
        method    = method,
        tuneGrid  = tune_grid,
        trControl = ctrl,
        weights   = wts,
        metric    = metric
      )
      
      score <- max(mod$results[[metric]], na.rm = TRUE)
      if (is.na(score) || is.infinite(score)) score <- 0.5
      list(Score = score, Pred = 0)
      
    }, error = function(e) {
      if (verbose) message("\n[Training Warning]: ", e$message)
      # Return chance-level score (0.5 for ROC) instead of -1e6 to preserve GP regression stability
      return(list(Score = 0.5, Pred = 0)) 
    })
    
    return(res)
  }
  
  # Step 7: Execute Bayesian Optimization
  if (verbose) cat("\n>>> Starting Bayesian Optimization process...\n")
  opt_res <- rBayesianOptimization::BayesianOptimization(
    FUN         = obj_func,
    bounds      = bounds,
    init_points = init_points,
    n_iter      = n_iter,
    acq         = "ucb", 
    kappa       = 2.576,
    verbose     = verbose
  )
  
  # Step 8: Retrain final optimal model
  best_params <- as.list(opt_res$Best_Par)
  for (p in intersect(names(best_params), int_params)) {
    best_params[[p]] <- round(best_params[[p]])
  }
  final_grid <- do.call(data.frame, best_params)
  
  ctrl_final <- caret::trainControl(
    method          = "cv",
    number          = cv_folds,
    classProbs      = TRUE,
    summaryFunction = summaryFun,
    sampling        = sampling
  )
  
  final_model <- caret::train(
    as.formula(paste(gc, "~ .")),
    data      = train_data,
    method    = method,
    tuneGrid  = final_grid,
    trControl = ctrl_final,
    weights   = wts,
    metric    = metric
  )
  
  model_obj@best.model.result$fine_tuned_model <- final_model
  model_obj@best.model.result$tuning_result      <- opt_res
  
  return(model_obj)
}
# =============================================================================
# Helper: Plot tuning history
# =============================================================================

#' Plot Bayesian Optimization History
#'
#' Shows the best score found so far across iterations.
#'
#' @param model_obj A \code{Train_Model} object that has been fine-tuned.
#' @param save_plot Logical.
#' @param save_dir  Output directory.
#' @return A ggplot.
#' @export
PlotTuningHistory <- function(model_obj, save_plot = FALSE, save_dir = NULL) {
  opt_res <- model_obj@best.model.result$tuning_result
  if (is.null(opt_res)) stop("No tuning result found. Run FineTuneModel first.")
  
  hist_dt <- opt_res$History
  if (!is.data.frame(hist_dt) || nrow(hist_dt) == 0) stop("No tuning history available.")
  
  scores <- hist_dt$Value
  valid_scores <- ifelse(scores <= 0.5 & max(scores, na.rm = TRUE) > 0.5, NA_real_, scores)
  best_curve <- cummax(ifelse(is.na(valid_scores), -Inf, valid_scores))
  best_curve[is.infinite(best_curve)] <- 0.5
  
  hist_df <- data.frame(Iteration = seq_along(best_curve), Best_Score = best_curve)
  
  p <- ggplot2::ggplot(hist_df, ggplot2::aes(x = Iteration, y = Best_Score)) +
    ggplot2::geom_line(color = "#b2e2e2", linewidth = 1) +
    ggplot2::geom_point(color = "#006d2c", size = 2) +
    ggplot2::labs(title = "Bayesian Optimization History", x = "Iteration", y = "Best ROC Score") +
    .pub_theme(13)
  
  print(p)
  if (save_plot && !is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuning_history.pdf"), plot = p, width = 6, height = 4, dpi = 300)
  }
  return(p)
}


# Internal theme helper (shared with viz_functions.R)
.pub_theme <- function(base_size = 13) {
  ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.title    = ggplot2::element_text(face = "bold"),
      legend.title  = ggplot2::element_text(face = "bold"),
      strip.text    = ggplot2::element_text(face = "bold")
    )
}

.save_plot <- function(p, dir, filename, width, height, format = "pdf") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, paste0(tools::file_path_sans_ext(filename), ".", format))
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
  cat("Plot saved to:", path, "\n")
}

# =============================================================================
# Tuned Model Visualisation Functions
# =============================================================================
#' Plot ROC Curve for Tuned Model
#'
#' @description Generates an ROC curve for the fine-tuned model and optionally 
#' compares it with the original best model. Automatically handles factor level 
#' alignment for numeric-origin groups.
#'
#' @param tuned_model A caret \code{train} object (the fine-tuned model).
#' @param original_best A caret \code{train} object for comparison. Default is NULL.
#' @param test_data A data frame containing the test set.
#' @param group_col Character. Name of the grouping/target column.
#' @param palette_name Character. Wesanderson palette name. Default "Darjeeling1".
#' @param save_plot Logical. Whether to save the plot as a PDF.
#' @param save_dir Character. Output directory for the saved plot.
#' @param width,height Numeric. Plot dimensions in inches.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @importFrom pROC roc auc
#' @importFrom ggplot2 ggplot aes geom_line geom_abline scale_color_manual labs coord_equal theme element_text ggsave
#' @importFrom wesanderson wes_palette
#' @importFrom ggprism theme_prism
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   set.seed(123)
#'   idx <- sample(1:nrow(mtcars), 20)
#'   model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#'   trained <- ModelTrainAnalysis(model, methods = c("glm"), 
#'   control = list(method = "cv", number = 3), save_plots = FALSE)
#'   # Assume we have a tuned model (from FineTuneModel)
#'   # For demo, use the trained model as tuned (not actually tuned)
#'   PlotTunedROC(trained@train.models$glm, test_data = trained@filtered.set$testing, 
#'   group_col = "am", save_plot = FALSE)
#'   PlotTunedConfusion(trained@train.models$glm, test_data = trained@filtered.set$testing, 
#'   group_col = "am", save_plot = FALSE)
#'   PlotTunedCalibration(trained@train.models$glm, test_data = trained@filtered.set$testing, 
#'   group_col = "am", save_plot = FALSE)
#' }
#' }
PlotTunedROC <- function(tuned_model,
                         original_best = NULL,
                         test_data,
                         group_col,
                         palette_name = "Darjeeling1",
                         save_plot = FALSE,
                         save_dir = NULL,
                         width = 7,
                         height = 6) {
  
  if (!inherits(tuned_model, "train")) stop("tuned_model must be a caret train object.")
  
  # Align test labels with the model's make.names() transformation
  test_labels <- as.factor(test_data[[group_col]])
  levels(test_labels) <- make.names(levels(test_labels))
  
  # Helper to compute ROC data
  .compute_roc <- function(model, model_label) {
    # Get probabilities for the second class
    probs <- stats::predict(model, newdata = test_data, type = "prob")[, 2]
    
    roc_obj <- pROC::roc(test_labels, probs, 
                         levels = levels(test_labels), 
                         direction = "auto", quiet = TRUE)
    
    auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
    data.frame(
      Sensitivity = roc_obj$sensitivities,
      Specificity = roc_obj$specificities,
      Model = paste0(model_label, " (AUC = ", auc_val, ")")
    )
  }
  
  df_tuned <- .compute_roc(tuned_model, "Tuned")
  
  if (!is.null(original_best)) {
    df_best <- .compute_roc(original_best, "Original Best")
    roc_df <- rbind(df_tuned, df_best)
  } else {
    roc_df <- df_tuned
  }
  
  n_models <- length(unique(roc_df$Model))
  cols <- wesanderson::wes_palette(palette_name, max(2, n_models), type = "discrete")
  
  p <- ggplot2::ggplot(roc_df, 
                       ggplot2::aes(x = 1 - .data$Specificity, y = .data$Sensitivity, color = .data$Model)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(title = "ROC Curve: Tuned vs Original",
                  x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::coord_equal() +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = c(0.75, 0.25))
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuned_ROC.pdf"), plot = p,
                    width = width, height = height, dpi = 300)
  }
  
  return(p)
}

#' Plot Confusion Matrix for Tuned Model
#'
#' @description Generates a heatmap-style confusion matrix. Percentages are 
#' calculated per actual class (column-wise).
#'
#' @param tuned_model A caret \code{train} object.
#' @param test_data A data frame containing the test set.
#' @param group_col Character. Name of the grouping column.
#' @param palette Character vector of two colors for the gradient.
#' @param save_plot Logical.
#' @param save_dir Character.
#' @param width,height Numeric.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient labs
PlotTunedConfusion <- function(tuned_model,
                               test_data,
                               group_col,
                               palette = c("#b2e2e2", "#006d2c"),
                               save_plot = FALSE,
                               save_dir = NULL,
                               width = 5,
                               height = 4.5) {
  
  if (!inherits(tuned_model, "train")) stop("tuned_model must be a caret train object.")
  
  # Align labels
  truth <- as.factor(test_data[[group_col]])
  levels(truth) <- make.names(levels(truth))
  
  pred <- stats::predict(tuned_model, newdata = test_data, type = "raw")
  
  cf <- table(Predicted = factor(pred, levels = levels(truth)),
              Actual = truth)
  
  cf_df <- as.data.frame(cf)
  # Calculate percentages per Actual class
  cf_df <- cf_df %>%
    dplyr::group_by(.data$Actual) %>%
    dplyr::mutate(Pct = round(.data$Freq / sum(.data$Freq) * 100, 1)) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(cf_df, ggplot2::aes(x = .data$Actual, y = .data$Predicted, fill = .data$Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(.data$Freq, "\n(", .data$Pct, "%)")),
                       size = 4.5, fontface = "bold") +
    ggplot2::scale_fill_gradient(low = palette[1], high = palette[2]) +
    ggplot2::labs(title = "Confusion Matrix: Tuned Model",
                  x = "Actual Class", y = "Predicted Class", fill = "Count") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuned_confusion.pdf"), plot = p,
                    width = width, height = height, dpi = 300)
  }
  
  return(p)
}

#' Plot Calibration Curve for a Tuned Caret Model (Enhanced Binning & Smoothing)
#'
#' @param tuned_model A caret \code{train} object for a binary classifier.
#' @param test_data Data frame containing the test set.
#' @param group_col Name of the column in \code{test_data} with true class labels.
#' @param n_bins Number of bins for calibration. If \code{NULL}, determined automatically.
#' @param bin_method Binning strategy: \code{"auto"} (default, chooses best based on data),
#'   \code{"quantile"} (equal sample size per bin), or \code{"equal_width"} (equal width).
#' @param palette Color for calibration points and smooth curve.
#' @param base_size Base font size for plot.
#' @param se Logical; if \code{TRUE}, show confidence band around smooth.
#' @param show_stats Logical; display metrics on plot.
#' @param boot_ci Logical; compute bootstrap confidence intervals for metrics.
#' @param boot_n Number of bootstrap replicates.
#' @param ci_level Confidence level for bootstrap intervals.
#' @param seed Random seed for reproducibility.
#' @param save_plot Logical; save plot to file.
#' @param save_dir Directory to save plot.
#' @param width,height Dimensions of saved plot.
#' @param smooth_method Method for smoothing curve: \code{"lm"} (linear logit),
#'   \code{"loess"} (local regression), or \code{"none"} (no smooth line).
#' @return Invisibly returns a list with plot, metrics, and binned summary.
#' @export
PlotTunedCalibration <- function(tuned_model,
                                 test_data,
                                 group_col,
                                 n_bins = NULL,
                                 bin_method = c("auto", "quantile", "equal_width"),
                                 palette = "#969696",
                                 base_size = 13,
                                 se = FALSE,
                                 show_stats = TRUE,
                                 boot_ci = TRUE,
                                 boot_n = 1000,
                                 ci_level = 0.95,
                                 seed = 1,
                                 save_plot = FALSE,
                                 save_dir = NULL,
                                 width = 6,
                                 height = 5.5,
                                 smooth_method = c("lm", "loess", "none")) {
  
  bin_method <- match.arg(bin_method)
  smooth_method <- match.arg(smooth_method)
  
  # ---- Input validation ----
  if (!inherits(tuned_model, "train")) {
    stop("`tuned_model` must be a caret `train` object.")
  }
  if (is.null(tuned_model$modelType) || tuned_model$modelType != "Classification") {
    stop("`tuned_model` must be a classification model trained with `classProbs = TRUE`.")
  }
  if (is.null(tuned_model$levels) || length(tuned_model$levels) != 2) {
    stop("`tuned_model` must be a binary classifier; ",
         "`tuned_model$levels` must contain exactly two class labels.")
  }
  if (!group_col %in% names(test_data)) {
    stop("`group_col` ('", group_col, "') was not found in `test_data`.")
  }
  
  # Drop missing outcomes
  n_missing_outcome <- sum(is.na(test_data[[group_col]]))
  if (n_missing_outcome > 0) {
    warning(n_missing_outcome, " row(s) with missing `", group_col, "` were dropped.")
    test_data <- test_data[!is.na(test_data[[group_col]]), , drop = FALSE]
  }
  
  # Encode truth
  truth_raw <- test_data[[group_col]]
  truth_labels <- make.names(as.character(truth_raw))
  observed_levels <- unique(truth_labels)
  if (length(observed_levels) != 2) {
    stop("`group_col` must be binary; found ", length(observed_levels),
         " distinct level(s) in `test_data`: ", paste(observed_levels, collapse = ", "))
  }
  model_levels <- tuned_model$levels
  pos_level <- model_levels[2]
  if (!pos_level %in% observed_levels) {
    stop("Positive class '", pos_level, "' (from `tuned_model$levels`) was ",
         "not found among the levels of `test_data[[group_col]]` (",
         paste(observed_levels, collapse = ", "), "). Check that `group_col` ",
         "uses the same coding as the outcome the model was trained on.")
  }
  truth_numeric <- as.integer(truth_labels == pos_level)
  
  # Predict probabilities
  prob_mat <- stats::predict(tuned_model, newdata = test_data, type = "prob")
  if (!pos_level %in% colnames(prob_mat)) {
    stop("Predicted probability matrix does not contain a column named '",
         pos_level, "'. Available columns: ", paste(colnames(prob_mat), collapse = ", "))
  }
  probs <- prob_mat[[pos_level]]
  eps <- 1e-06
  probs_clip <- pmax(pmin(probs, 1 - eps), eps)
  
  cal_df <- data.frame(truth = truth_numeric, prob = probs, prob_clip = probs_clip)
  
  # ---- Intelligent binning (based on PlotCalibration) ----
  n_obs <- nrow(cal_df)
  
  # If n_bins not provided, compute a reasonable default
  if (is.null(n_bins)) {
    base_bins <- min(10, max(4, floor(n_obs / 15)))
    if (bin_method == "auto") {
      n_unique <- length(unique(round(probs, 2)))
      spread <- diff(stats::quantile(probs, probs = c(0.1, 0.9), na.rm = TRUE))
      # Use quantile binning if probabilities are concentrated or have few unique values
      if (spread < 0.15 || n_unique < 10) {
        bin_method <- "quantile"
        n_bins <- min(6, max(3, base_bins - 1))
        message("Auto-switched to quantile binning (spread = ", round(spread, 3),
                ", n_unique = ", n_unique, "), n_bins = ", n_bins)
      } else {
        bin_method <- "equal_width"
        n_bins <- base_bins
        message("Auto-selected equal-width binning, n_bins = ", n_bins)
      }
    } else {
      n_bins <- base_bins
    }
  } else {
    # User provided n_bins; we may still adjust if bin_method == "auto"?
    if (bin_method == "auto") {
      n_unique <- length(unique(round(probs, 2)))
      spread <- diff(stats::quantile(probs, probs = c(0.1, 0.9), na.rm = TRUE))
      if (spread < 0.15 || n_unique < 10) {
        bin_method <- "quantile"
        message("Auto-switched to quantile binning (spread = ", round(spread, 3),
                ", n_unique = ", n_unique, "), n_bins = ", n_bins)
      } else {
        bin_method <- "equal_width"
        message("Auto-selected equal-width binning, n_bins = ", n_bins)
      }
    }
  }
  
  # Prevent too many bins relative to unique probabilities
  n_unique_prob <- length(unique(probs))
  if (n_bins > n_unique_prob) {
    warning("Requested n_bins (", n_bins, ") exceeds number of unique predicted probabilities (",
            n_unique_prob, "). Reducing to ", n_unique_prob, ".")
    n_bins <- min(n_bins, n_unique_prob)
  }
  
  # Create breaks
  if (bin_method == "quantile") {
    breaks <- stats::quantile(probs, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
    breaks <- unique(breaks)
    # Ensure full coverage from 0 to 1
    if (min(breaks) > 0) breaks <- c(0, breaks)
    if (max(breaks) < 1) breaks <- c(breaks, 1)
  } else { # equal_width
    breaks <- seq(0, 1, length.out = n_bins + 1)
  }
  breaks <- unique(sort(breaks))
  
  # Bin the data
  bin <- cut(cal_df$prob, breaks = breaks, include.lowest = TRUE)
  cal_df_bin <- data.frame(cal_df, bin = bin)
  cal_sum <- dplyr::summarise(
    dplyr::group_by(cal_df_bin, .data$bin),
    mean_pred = mean(.data$prob),
    obs_rate = mean(.data$truth),
    n = dplyr::n(),
    .groups = "drop"
  )
  cal_sum <- cal_sum[!is.na(cal_sum$bin), , drop = FALSE]
  
  # ---- Metrics computation (internal binning reused for ECE) ----
  compute_metrics <- function(df, n_bins, bin_method) {
    brier <- mean((df$truth - df$prob)^2)
    cal_glm <- tryCatch(
      suppressWarnings(stats::glm(truth ~ log(prob_clip/(1 - prob_clip)),
                                  family = stats::binomial(), data = df)),
      error = function(e) NULL
    )
    if (is.null(cal_glm)) {
      intercept <- NA_real_
      slope <- NA_real_
    } else {
      intercept <- unname(stats::coef(cal_glm)[1])
      slope <- unname(stats::coef(cal_glm)[2])
    }
    # Re-bin for ECE (using same logic)
    if (bin_method == "quantile") {
      br <- stats::quantile(df$prob, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
      br <- unique(br)
      if (min(br) > 0) br <- c(0, br)
      if (max(br) < 1) br <- c(br, 1)
    } else {
      br <- seq(0, 1, length.out = n_bins + 1)
    }
    br <- unique(sort(br))
    bin_cut <- cut(df$prob, breaks = br, include.lowest = TRUE)
    df_bin <- data.frame(df, bin = bin_cut)
    bin_sum <- dplyr::summarise(
      dplyr::group_by(df_bin, .data$bin),
      mean_pred = mean(.data$prob),
      obs_rate = mean(.data$truth),
      n = dplyr::n(),
      .groups = "drop"
    )
    bin_sum <- bin_sum[!is.na(bin_sum$bin), , drop = FALSE]
    ece <- sum(bin_sum$n * abs(bin_sum$obs_rate - bin_sum$mean_pred)) / sum(bin_sum$n)
    list(brier = brier, intercept = intercept, slope = slope, ece = ece, bin_sum = bin_sum)
  }
  
  # Point estimates
  point_est <- compute_metrics(cal_df, n_bins, bin_method)
  metrics_df <- data.frame(
    statistic = c("Brier", "Intercept", "Slope", "ECE"),
    estimate = c(point_est$brier, point_est$intercept, point_est$slope, point_est$ece),
    stringsAsFactors = FALSE
  )
  
  # Bootstrap CIs (if requested)
  if (boot_ci) {
    if (!is.null(seed)) set.seed(seed)
    alpha <- 1 - ci_level
    boot_mat <- matrix(NA_real_, nrow = boot_n, ncol = 4,
                       dimnames = list(NULL, c("brier", "intercept", "slope", "ece")))
    for (b in seq_len(boot_n)) {
      idx <- sample.int(n_obs, n_obs, replace = TRUE)
      boot_res <- tryCatch(
        compute_metrics(cal_df[idx, , drop = FALSE], n_bins, bin_method),
        error = function(e) NULL
      )
      if (!is.null(boot_res)) {
        boot_mat[b, ] <- c(boot_res$brier, boot_res$intercept,
                           boot_res$slope, boot_res$ece)
      }
    }
    ci <- apply(boot_mat, 2, stats::quantile, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    metrics_df$ci_lower <- ci[1, ]
    metrics_df$ci_upper <- ci[2, ]
  }
  
  # ---- Plotting with flexible smoothing ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey40", linewidth = 0.8)
  
  # Add smooth curve if not "none"
  if (smooth_method != "none") {
    p <- p + ggplot2::geom_smooth(
      data = cal_df,
      ggplot2::aes(x = prob, y = truth),
      method = smooth_method,
      se = se,
      color = palette,
      linewidth = 1,
      fill = paste0(palette, "30"),
      formula = if (smooth_method == "lm") y ~ x else y ~ x   # formula same for both
    )
  }
  
  p <- p + ggplot2::geom_point(
    data = cal_sum,
    ggplot2::aes(x = mean_pred, y = obs_rate, size = n),
    color = palette, alpha = 0.8
  ) +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "Samples (n)") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(
      title = "Calibration Curve: Tuned Model",
      subtitle = paste0("Model: ", tuned_model$method, " | Bins: ", n_bins,
                        " (", bin_method, ") | Smooth: ", smooth_method),
      x = "Mean Predicted Probability",
      y = "Observed Proportion"
    ) +
    .pub_theme(base_size)   # ensure this theme exists in your package
  
  if (show_stats) {
    fmt_line <- function(label, row) {
      if (boot_ci && !is.na(row$ci_lower)) {
        sprintf("%s: %.3f [%.3f, %.3f]", label, row$estimate,
                row$ci_lower, row$ci_upper)
      } else {
        sprintf("%s: %.3f", label, row$estimate)
      }
    }
    stats_text <- paste(
      fmt_line("Brier", metrics_df[metrics_df$statistic == "Brier", ]),
      fmt_line("Intercept", metrics_df[metrics_df$statistic == "Intercept", ]),
      fmt_line("Slope", metrics_df[metrics_df$statistic == "Slope", ]),
      fmt_line("ECE", metrics_df[metrics_df$statistic == "ECE", ]),
      sep = "\n"
    )
    p <- p + ggplot2::annotate("text", x = 0.05, y = 0.95,
                               label = stats_text, hjust = 0, vjust = 1,
                               size = base_size * 0.24, family = "mono", fontface = "bold")
  }
  
  print(p)
  
  if (save_plot) {
    if (is.null(save_dir)) {
      stop("`save_dir` must be provided when `save_plot = TRUE`.")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuned_calibration.pdf"),
                    plot = p, width = width, height = height, dpi = 300)
  }
  
  invisible(list(plot = p, metrics = metrics_df, binned_summary = cal_sum))
}

