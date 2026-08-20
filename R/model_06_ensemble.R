# =============================================================================
# model_final_ensemble.R
# Final Modeling Module: Ensemble -Publication Visualizations
# =============================================================================

# 0. Package Check -----------------------------------------------------------
.check_final_pkgs <- function() {
  required <- c("caret", "caretEnsemble", "ggplot2", "wesanderson", "ggprism",
                "pROC", "rBayesianOptimization", "xgboost", "randomForest",
                "glmnet", "gbm", "dplyr", "tidyr", "reshape2", "gridExtra")
  missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Install them before using this module.")
  }
  invisible(TRUE)
}

# -- Internal helpers (shared with viz_functions) ---------------------------
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

.get_palette <- function(palette_name, n) {
  tryCatch(
    as.character(wesanderson::wes_palette(n = n, name = palette_name,
                                          type = if (n > 5) "continuous" else "discrete")),
    error = function(e) RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)]
  )
}

.save_plot <- function(p, dir, filename, width, height, format = "pdf") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, paste0(tools::file_path_sans_ext(filename), ".", format))
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
  cat("Plot saved to:", path, "\n")
  invisible(path)
}


# =============================================================================
# Part 1: No-Hyperparameter Comparison Visualizations
# =============================================================================

#' Top N Model Performance Dot Plot
#'
#' Shows AUC for the top N models with 95% CI as error bars, sorted by AUC.
#'
#' @param model_obj A trained \code{Train_Model} object.
#' @param top_n Number of top models to show. Default 6.
#' @param save_plot Logical.
#' @param save_dir Output directory.
#' @export
PlotTopModelAUC <- function(model_obj, top_n = 6, save_plot = FALSE, save_dir = NULL) {
  .check_final_pkgs()
  
  perf <- model_obj@all.results
  if (nrow(perf) == 0) stop("No performance data available.")
  
  perf <- head(perf[order(-perf$auc), ], top_n)
  perf$Model <- factor(perf$Model, levels = rev(perf$Model))
  
  has_ci <- all(c("auc_lower", "auc_upper") %in% colnames(perf))
  cols <- .get_palette("Darjeeling1", nrow(perf))
  
  p <- ggplot2::ggplot(perf, ggplot2::aes(x = auc, y = Model, color = Model)) +
    ggplot2::geom_point(size = 4, shape = 18) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(title = "Top Models by AUC", x = "AUC (95% CI)", y = NULL) +
    .pub_theme(13) +
    ggplot2::theme(legend.position = "none")
  
  # Only draw error bars if real confidence intervals exist
  if (has_ci) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = auc_lower, xmax = auc_upper),
      height = 0.2, linewidth = 0.8
    )
  }
  
  if (save_plot && !is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "top_model_auc.pdf"), plot = p, width = 6, height = 4, dpi = 300)
  }
  return(p)
}

#' Multi-Metric Faceted Bar Chart
#'
#' @param model_obj A trained Train_Model object.
#' @param top_n Number of top models. Default 4.
#' @param metrics Metrics to display. Default: c("auc", "Sensitivity", 
#'   "Specificity", "accuracy_score", "f1_score").
#' @param palette_name Wesanderson palette. Default "Darjeeling1".
#' @param base_size Base font size. Default 13.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output directory.
#' @param width,height Plot dimensions.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' set.seed(123)
#' idx <- sample(1:nrow(mtcars), 20)
#' model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#' trained <- ModelTrainAnalysis(model, methods = c("glm", "rf"),
#'                               control = list(method = "cv", number = 3),
#'                               save_plots = FALSE)
#' if (interactive()) {
#'   PlotModelComparison(trained, top_n = 2, save_plot = FALSE)
#'   PlotModelHeatmap(trained, save_plot = FALSE)
#'   PlotModelParallel(trained, top_n = 2, save_plot = FALSE)
#' }
#' }
PlotModelComparison <- function(model_obj, 
                                top_n = 4,
                                metrics = c("auc", "Sensitivity", "Specificity", "accuracy_score", "f1_score"),
                                palette_name = "Darjeeling1",
                                base_size = 13,
                                save_plot = FALSE,
                                save_dir = NULL,
                                width = 10,
                                height = 6) {
  .check_final_pkgs()
  perf <- model_obj@all.results
  if (nrow(perf) == 0) stop("No model performance data available.")
  
  available_metrics <- intersect(metrics, colnames(perf))
  if (length(available_metrics) < 1) stop("No valid metrics found in model results.")
  
  perf <- head(perf[order(-perf$auc), ], top_n)
  
  bar_long <- perf %>%
    dplyr::select(dplyr::all_of(c("Model", available_metrics))) %>%
    tidyr::pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")
  
  # Safe dictionary lookup with fallback to original metric name
  metric_dict <- c(
    auc = "AUC", accuracy_score = "Accuracy", f1_score = "F1 Score",
    Sensitivity = "Sensitivity", Specificity = "Specificity",
    Precision = "Precision", recall_score = "Recall"
  )
  
  bar_long$Metric_Label <- ifelse(
    bar_long$Metric %in% names(metric_dict),
    metric_dict[bar_long$Metric],
    bar_long$Metric
  )
  bar_long$Metric_Label <- factor(bar_long$Metric_Label, levels = unique(bar_long$Metric_Label))
  
  cols <- .get_palette(palette_name, length(unique(bar_long$Model)))
  
  p <- ggplot2::ggplot(bar_long, ggplot2::aes(x = Model, y = Value, fill = Model)) +
    ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(Value), "N/A", sprintf("%.3f", Value))),
      vjust = -0.5, size = 3, fontface = "bold"
    ) +
    ggplot2::facet_wrap(~ Metric_Label, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = "Model Performance Comparison", subtitle = paste("Top", top_n, "models by AUC"), x = NULL, y = NULL) +
    .pub_theme(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  print(p)
  if (save_plot && !is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "model_comparison.pdf"), plot = p, width = width, height = height, dpi = 300)
  }
  return(p)
}

#' Model Performance Heatmap
#'
#' @param model_obj A trained Train_Model object.
#' @param metrics Metrics to display.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output directory.
#' @return A ggplot object.
#' @export
PlotModelHeatmap <- function(model_obj,
                             metrics = c("auc", "Sensitivity", "Specificity", 
                                         "accuracy_score", "f1_score", "Precision"),
                             save_plot = FALSE,
                             save_dir = NULL) {
  
  cat("Generating model performance heatmap...\n")
  
  perf <- model_obj@all.results
  available_metrics <- intersect(metrics, colnames(perf))
  
  perf[, available_metrics] <- lapply(perf[, available_metrics], function(x) {
    ifelse(is.na(x), NA_real_, x)
  })
  
  heat_long <- perf %>%
    dplyr::select(dplyr::all_of(c("Model", available_metrics))) %>%
    tidyr::pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")
  
  p <- ggplot2::ggplot(heat_long, 
                       ggplot2::aes(x = Metric, y = Model, fill = Value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(Value), "NA", sprintf("%.3f", Value))),
      colour = "white", fontface = "bold", size = 4
    ) +
    ggplot2::scale_fill_viridis_c(option = "D", na.value = "grey80", 
                                  limits = c(0, 1), name = "Score") +
    ggplot2::labs(title = "Model Performance Heatmap", x = NULL, y = NULL) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, "model_heatmap.pdf")
    ggplot2::ggsave(path, plot = p, width = 8, height = 5, dpi = 300)
    cat("Plot saved to:", path, "\n")
  }
  
  return(p)
}

#' Model Performance Parallel Coordinates
#'
#' @param model_obj A trained Train_Model object.
#' @param top_n Number of top models. Default 5.
#' @param metrics Metrics to display.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output directory.
#' @return A ggplot object.
#' @export
PlotModelParallel <- function(model_obj,
                              top_n = 5,
                              metrics = c("auc", "Sensitivity", "Specificity", 
                                          "accuracy_score", "f1_score"),
                              save_plot = FALSE,
                              save_dir = NULL) {
  
  cat("Generating parallel coordinates plot...\n")
  
  perf <- model_obj@all.results
  available_metrics <- intersect(metrics, colnames(perf))
  
  perf <- head(perf[order(-perf$auc), ], top_n)
  perf[, available_metrics] <- lapply(perf[, available_metrics], function(x) {
    ifelse(is.na(x), 0, x)
  })
  
  cols <- wesanderson::wes_palette("Darjeeling1", top_n, type = "discrete")
  
  # Normalize to 0-1 for parallel coordinates
  norm_perf <- perf
  norm_perf[, available_metrics] <- apply(norm_perf[, available_metrics], 2, 
                                          function(x) (x - min(x)) / (max(x) - min(x) + 1e-8))
  
  parallel_long <- norm_perf %>%
    dplyr::select(dplyr::all_of(c("Model", available_metrics))) %>%
    tidyr::pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")
  
  p <- ggplot2::ggplot(parallel_long, 
                       ggplot2::aes(x = Metric, y = Value, group = Model, color = Model)) +
    ggplot2::geom_line(linewidth = 1.2, alpha = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      title = "Model Performance Parallel Coordinates",
      subtitle = paste("Top", top_n, "models (metrics normalized 0-1)"),
      x = NULL, y = "Normalized Score"
    ) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, "model_parallel.pdf")
    ggplot2::ggsave(path, plot = p, width = 8, height = 5, dpi = 300)
    cat("Plot saved to:", path, "\n")
  }
  
  return(p)
}



#' Prediction Probability Density Plot (fixed)
#'
#' @param model_obj A trained Train_Model object.
#' @param data Test data (default: split.data$testing).
#' @param save_plot Logical.
#' @param save_dir Output directory.
#' @export
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
#'   PlotProbDensity(trained, data = trained@filtered.set$testing, save_plot = FALSE)
#'   PlotProbStrip(trained, data = trained@filtered.set$testing, save_plot = FALSE)
#' }
#' }
PlotProbDensity <- function(model_obj, data = NULL,
                            save_plot = FALSE, save_dir = NULL) {
  cat("Generating probability density plot...\n")
  
  # ---- Extract best model robustly ----
  best_model <- NULL
  if (!is.null(model_obj@best.model.result$model)) {
    best_model <- model_obj@best.model.result$model
  } else if (!is.null(model_obj@best.model.result$model_type)) {
    best_model <- model_obj@train.models[[model_obj@best.model.result$model_type]]
  }
  if (is.null(best_model)) {
    # fallback: use first trained model
    best_model <- model_obj@train.models[[1]]
    warning("No best model found; using first trained model: ", names(model_obj@train.models)[1])
  }
  
  if (is.null(data)) data <- model_obj@split.data$testing
  if (is.null(data)) stop("No test data available.")
  
  gc <- model_obj@group_col
  
  probs <- predict(best_model, data, type = "prob")[, 2]
  true <- data[[gc]]
  
  prob_df <- data.frame(Probability = probs, Class = true)
  
  cols <- c("#2c3e50", "#e74c3c")
  
  p <- ggplot2::ggplot(prob_df, ggplot2::aes(x = Probability, fill = Class)) +
    ggplot2::geom_density(alpha = 0.6, color = NA) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "Predicted Probability Density",
                  x = "Probability of Positive Class", y = "Density") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, "prob_density.pdf")
    ggplot2::ggsave(path, plot = p, width = 6, height = 4, dpi = 300)
    cat("Plot saved to:", path, "\n")
  }
  
  return(p)
}


#' Prediction Probability Strip Chart (fixed)
#'
#' @param model_obj A trained Train_Model object.
#' @param data Test data.
#' @param save_plot Logical.
#' @param save_dir Output directory.
#' @export
PlotProbStrip <- function(model_obj, data = NULL,
                          save_plot = FALSE, save_dir = NULL) {
  cat("Generating probability strip chart...\n")
  
  # ---- Extract best model robustly ----
  best_model <- NULL
  if (!is.null(model_obj@best.model.result$model)) {
    best_model <- model_obj@best.model.result$model
  } else if (!is.null(model_obj@best.model.result$model_type)) {
    best_model <- model_obj@train.models[[model_obj@best.model.result$model_type]]
  }
  if (is.null(best_model)) {
    best_model <- model_obj@train.models[[1]]
    warning("No best model found; using first trained model.")
  }
  
  if (is.null(data)) data <- model_obj@split.data$testing
  if (is.null(data)) stop("No test data available.")
  
  gc <- model_obj@group_col
  
  probs <- predict(best_model, data, type = "prob")[, 2]
  true <- data[[gc]]
  
  strip_df <- data.frame(Probability = probs, Class = true)
  strip_df <- strip_df[order(strip_df$Probability), ]
  strip_df$Index <- seq_len(nrow(strip_df))
  
  cols <- c("#2c3e50", "#e74c3c")
  
  p <- ggplot2::ggplot(strip_df, ggplot2::aes(x = Index, y = Probability, color = Class)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(title = "Predicted Probabilities (Sorted)",
                  x = "Sample (sorted by probability)", y = "Probability") +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, "prob_strip.pdf")
    ggplot2::ggsave(path, plot = p, width = 8, height = 3.5, dpi = 300)
    cat("Plot saved to:", path, "\n")
  }
  
  return(p)
}



# =============================================================================
# Part 2: caretEnsemble Integration
# =============================================================================

#' Multi-Strategy Model Ensemble (Robust to Single-Row Input)
#'
#' Combines multiple trained caret models using one of four strategies.
#' All predict_fn implementations now handle both single-row and multi-row
#' newdata without dimension errors.
#'
#' @param model_obj A trained Train_Model object.
#' @param strategy Ensemble strategy: "stacking", "average", "weighted", "voting".
#' @param meta_method Meta-learner method for stacking (default "glm").
#' @param top_n Number of top models (by AUC) to include. NULL uses all.
#' @param cv_folds Number of CV folds for stacking meta-features.
#' @param weights Named numeric vector for weighted strategy.
#' @param seed Random seed.
#' @return Updated model_obj with ensemble stored in @best.model.result$ensemble.
#' @export
TrainEnsemble <- function(model_obj,
                          strategy   = c("stacking", "average", "weighted", "voting"),
                          meta_method = "glm",
                          top_n       = NULL,
                          cv_folds    = 5,
                          weights     = NULL,
                          seed        = 123) {
  
  strategy <- match.arg(strategy)
  set.seed(seed)
  
  if (!inherits(model_obj, "Train_Model")) {
    stop("model_obj must be a Train_Model object.")
  }
  
  model_list <- model_obj@train.models
  if (length(model_list) < 2) stop("At least 2 trained models required.")
  
  if (!is.null(top_n)) {
    perf <- model_obj@all.results
    if (nrow(perf) == 0) stop("No performance data in @all.results.")
    cv_metrics <- model_obj@process.info$cv_metrics
    if (is.null(cv_metrics) || nrow(cv_metrics) == 0) {
      stop("No CV metrics found in model_obj@process.info$cv_metrics. ",
           "Run ModelTrainAnalysis() first to obtain cross-validated performance.")
    }
    best <- head(cv_metrics[order(-cv_metrics$CV_Mean), ]$Model, top_n)
    model_list <- model_list[intersect(best, names(model_list))]
    if (length(model_list) < 2) stop("Fewer than 2 models after top_n filtering.")
  }
  
  cat(sprintf("Building %s ensemble from %d models...\n", strategy, length(model_list)))
  
  ensemble <- switch(strategy,
                     
                     # -- Stacking --------------------------------------------------
                     stacking = {
                       train_data <- model_obj@filtered.set$training
                       if (is.null(train_data)) stop("filtered.set$training is empty.")
                       gc <- model_obj@group_col
                       outcome_orig <- train_data[[gc]]
                       if (!is.factor(outcome_orig)) outcome_orig <- as.factor(outcome_orig)
                       valid_levels <- make.names(levels(outcome_orig), unique = TRUE)
                       train_data[[gc]] <- factor(outcome_orig, levels = levels(outcome_orig), labels = valid_levels)
                       folds <- caret::createFolds(train_data[[gc]], k = cv_folds, list = TRUE)
                       
                       # Reliability check
                       reliable_models <- names(model_list)[sapply(names(model_list), function(nm) {
                         all(sapply(folds, function(fold_idx) {
                           train_fold <- train_data[-fold_idx, ]
                           test_fold  <- train_data[fold_idx, ]
                           m <- tryCatch(
                             caret::train(as.formula(paste(gc, "~ .")), data = train_fold,
                                          method = model_list[[nm]]$method,
                                          tuneGrid = model_list[[nm]]$bestTune,
                                          trControl = caret::trainControl(method = "none"),
                                          verbose = FALSE),
                             error = function(e) NULL
                           )
                           if (is.null(m)) return(FALSE)
                           preds <- tryCatch(
                             predict(m, test_fold, type = "prob")[, 2],
                             error = function(e) NULL
                           )
                           if (is.null(preds) || any(is.na(preds)) || length(preds) != nrow(test_fold))
                             return(FALSE)
                           TRUE
                         }))
                       })]
                       
                       if (length(reliable_models) < 2) {
                         stop("Fewer than 2 models passed reliability test. Use 'average' or 'voting'.")
                       }
                       cat("Reliable models for stacking:", paste(reliable_models, collapse = ", "), "\n")
                       
                       # Generate meta-features
                       meta_train <- matrix(NA, nrow = nrow(train_data), ncol = length(reliable_models))
                       colnames(meta_train) <- reliable_models
                       for (i in seq_along(folds)) {
                         fold_idx <- folds[[i]]
                         train_fold <- train_data[-fold_idx, ]
                         test_fold  <- train_data[fold_idx, ]
                         for (nm in reliable_models) {
                           m <- caret::train(as.formula(paste(gc, "~ .")), data = train_fold,
                                             method = model_list[[nm]]$method,
                                             tuneGrid = model_list[[nm]]$bestTune,
                                             trControl = caret::trainControl(method = "none"),
                                             verbose = FALSE)
                           meta_train[fold_idx, nm] <- predict(m, test_fold, type = "prob")[, 2]
                         }
                       }
                       meta_df <- as.data.frame(meta_train)
                       meta_df$outcome <- train_data[[gc]]
                       
                       meta_model <- caret::train(
                         outcome ~ ., data = meta_df,
                         method = meta_method,
                         trControl = caret::trainControl(
                           method = "cv", number = 5,
                           classProbs = TRUE, summaryFunction = caret::twoClassSummary
                         ),
                         metric = "ROC",
                         family = if(meta_method == "glm") "binomial"
                       )
                       
                       list(
                         predict_fn = function(newdata) {
                           base_preds <- sapply(reliable_models, function(nm) {
                             predict(model_list[[nm]], newdata, type = "prob")[, 2]
                           })
                           # Ensure matrix even for single row
                           if (is.vector(base_preds)) {
                             base_preds <- matrix(base_preds, nrow = 1,
                                                  dimnames = list(NULL, reliable_models))
                           }
                           base_preds <- as.data.frame(base_preds)
                           predict(meta_model, base_preds, type = "prob")[, 2]
                         },
                         object = list(base_models = model_list[reliable_models], meta_model = meta_model),
                         strategy = "Stacking",
                         method   = meta_method,
                         weights  = NULL
                       )
                     },
                     
                     # -- Simple Average ---------------------------------------------
                     average = {
                       list(
                         predict_fn = function(newdata) {
                           preds <- sapply(model_list, function(m) {
                             predict(m, newdata, type = "prob")[, 2]
                           })
                           if (is.vector(preds)) {
                             mean(preds)
                           } else {
                             rowMeans(preds)
                           }
                         },
                         object   = model_list,
                         strategy = "Average",
                         method   = "Simple Average",
                         weights  = rep(1/length(model_list), length(model_list))
                       )
                     },
                     
                     # -- Weighted Average ------------------------------------------
                     weighted = {
                       cv_metrics <- model_obj@process.info$cv_metrics
                       if (is.null(cv_metrics) || nrow(cv_metrics) == 0) {
                         stop("No CV metrics found. Run ModelTrainAnalysis() first.")
                       }
                       aucs <- sapply(names(model_list), function(nm) {
                         cv_metrics$CV_Mean[cv_metrics$Model == nm]
                       })
                       names(aucs) <- names(model_list)
                       if (is.null(weights)) {
                         weights <- aucs / sum(aucs, na.rm = TRUE)
                         weights[is.na(weights)] <- 0
                       }
                       names(weights) <- names(model_list)
                       
                       list(
                         predict_fn = function(newdata) {
                           preds <- sapply(names(model_list), function(nm) {
                             predict(model_list[[nm]], newdata, type = "prob")[, 2] * weights[nm]
                           })
                           if (is.vector(preds)) {
                             sum(preds)
                           } else {
                             rowSums(preds)
                           }
                         },
                         object   = model_list,
                         strategy = "Weighted",
                         method   = "Weighted Average",
                         weights  = weights
                       )
                     },
                     
                     # -- Majority Voting --------------------------------------------
                     voting = {
                       list(
                         predict_fn = function(newdata) {
                           votes <- sapply(model_list, function(m) {
                             as.character(predict(m, newdata, type = "raw"))
                           })
                           if (is.vector(votes)) {
                             # Single row: votes is a character vector
                             names(which.max(table(votes)))
                           } else {
                             apply(votes, 1, function(x) names(which.max(table(x))))
                           }
                         },
                         object   = model_list,
                         strategy = "Voting",
                         method   = "Majority Voting",
                         weights  = NULL
                       )
                     }
  )
  
  if (is.null(model_obj@best.model.result)) {
    model_obj@best.model.result <- list()
  }
  model_obj@best.model.result$ensemble <- ensemble
  
  cat(sprintf("Ensemble '%s' created.\n", ensemble$strategy))
  return(model_obj)
}


#' Predict Using a Stored Ensemble
#'
#' Calls the \code{predict_fn} stored inside the ensemble to generate
#' predictions on new data.
#'
#' @param model_obj A \code{Train_Model} object that has been processed by
#'   \code{TrainEnsemble}.
#' @param newdata    A data frame with the same predictor columns as the
#'   original training data.
#'
#' @return For classification: a numeric vector of class probabilities
#'   (stacking, average, weighted) or a character vector of predicted
#'   classes (voting).
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'ens' from TrainEnsemble
#' # newdata <- mtcars[1:5, ]
#' # probs <- PredictEnsemble(ens, newdata)
#' # print(probs)
#' }
PredictEnsemble <- function(model_obj, newdata) {
  ens <- model_obj@best.model.result$ensemble
  if (is.null(ens)) stop("No ensemble found. Run TrainEnsemble first.")
  ens$predict_fn(newdata)
}

