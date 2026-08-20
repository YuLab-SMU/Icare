# =============================================================================
# Downstream tooling for NON-BINARY ModelTrainAnalysis() results
# =============================================================================
# ModelTrainAnalysis() (see model_05_train.R) stores results differently
# depending on task_type:
#   - binary:      object@all.results is a data.frame (Model, auc, ...)
#   - multiclass/  object@all.results is a LIST:
#     regression:    list(task_type, cv_selection, train_final, test_final,
#                         preProcess, imbalance_handling)
#                   and object@best.model.result is a list with $model,
#                   $model_type, $task_type, $cv_metric, $test_performance,
#                   $train_performance, $selection_basis.
#
# The functions below are written specifically against that LIST shape, so
# that multiclass/regression pipelines have the same kind of "summarize /
# compare / visualize" tooling that PlotModelComparison(), PlotConfusionMatrix(),
# SelectBestModel(), etc. provide for binary results -- WITHOUT requiring any
# change to what ModelTrainAnalysis() itself returns.
#
# Do NOT call SelectBestModel() on a non-binary object -- it assumes
# object@all.results is a data.frame and will stop with an informative error
# if you try (see model_05_train.R). Model selection for multiclass/
# regression is already done for you inside ModelTrainAnalysis(); these
# functions just report on / visualize that selection.
# =============================================================================

.assert_nonbinary_result <- function(object) {
  if (!inherits(object, "Train_Model")) {
    stop("object must be a 'Train_Model' instance produced by ModelTrainAnalysis().")
  }
  res <- slot(object, "all.results")
  if (!is.list(res) || is.data.frame(res) || is.null(res$task_type)) {
    stop("object@all.results does not look like a multiclass/regression result ",
         "from ModelTrainAnalysis(). For binary results use SelectBestModel(), ",
         "PlotModelComparison(), PlotConfusionMatrix(), etc. instead.")
  }
  res
}

#' Summarize a Multiclass or Regression ModelTrainAnalysis Result
#'
#' Prints (and returns invisibly) a tidy overview of the cross-validation
#' ranking and the final train/test performance of the selected best model,
#' for a `Train_Model` object produced by `ModelTrainAnalysis()` with
#' `task_type = "multiclass"` or `"regression"`.
#'
#' @param object A `Train_Model` object (non-binary result).
#' @return Invisibly, a list with `cv_ranking`, `best_model`, `train_performance`,
#'   `test_performance`.
#' @export
#' @examples
#' \dontrun{
#' data(iris)
#' iris$Species <- as.factor(iris$Species)
#' model <- CreateModelObject(data = iris, group_col = "Species")
#' set.seed(123)
#' idx <- sample(1:nrow(iris), 100)
#' model@filtered.set <- list(training = iris[idx, ], testing = iris[-idx, ])
#' trained <- ModelTrainAnalysis(model,
#'                               methods = c("multinom", "rf"),
#'                               control = list(method = "cv", number = 3),
#'                               task_type = "multiclass",
#'                               save_plots = FALSE)
#' SummarizeNonBinaryResults(trained)
#' }
SummarizeNonBinaryResults <- function(object) {
  res <- .assert_nonbinary_result(object)

  cat("========== Task type:", res$task_type, "==========\n\n")
  cat("--- Cross-validation ranking (", object@process.info$task_type %||% res$task_type, ") ---\n")
  print(res$cv_selection)

  best <- object@best.model.result
  cat("\n--- Selected model:", best$model_type, "---\n")
  cat("Selection basis:", best$selection_basis, "\n\n")

  cat("Train-set performance:\n")
  print(best$train_performance)
  cat("\nTest-set performance:\n")
  print(best$test_performance)

  if (res$task_type == "regression") {
    gap <- best$train_performance$RMSE - best$test_performance$RMSE
    cat("\nTrain vs test RMSE gap:", round(gap, 4),
        if (abs(gap) > 0.1 * best$test_performance$RMSE) " (large -- check for overfitting)\n" else " (acceptable)\n")
  } else {
    gap <- best$train_performance$accuracy - best$test_performance$accuracy
    cat("\nTrain vs test accuracy gap:", round(gap, 4),
        if (gap > 0.1) " (large -- check for overfitting)\n" else " (acceptable)\n")
  }

  invisible(list(
    cv_ranking = res$cv_selection,
    best_model = best$model_type,
    train_performance = best$train_performance,
    test_performance = best$test_performance
  ))
}

# small null-coalescing helper (base R has no built-in %||%)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Plot Cross-Validation Comparison Across Models (any task_type)
#'
#' Bar chart (mean +/- SD of the selection metric) comparing every trained
#' model's cross-validation performance. Works for `Train_Model` objects
#' produced by `ModelTrainAnalysis()` regardless of `task_type` -- it reads
#' `object@process.info$cv_metrics` (binary) or
#' `object@all.results$cv_selection` (multiclass/regression), whichever is
#' present.
#'
#' @param object A `Train_Model` object.
#' @param palette_name Wes Anderson palette name. Default `"AsteroidCity1"`.
#' @param save_plot Logical. Save the plot as PDF. Default `FALSE`.
#' @param save_dir Directory to save to. Required if `save_plot = TRUE`.
#' @param plot_width,plot_height Plot dimensions in inches.
#' @return Invisibly, the ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar scale_fill_manual labs theme element_text
#' @importFrom ggprism theme_prism
#' @export
PlotCVComparison <- function(object, palette_name = "AsteroidCity1",
                             save_plot = FALSE, save_dir = NULL,
                             plot_width = 6, plot_height = 5) {
  if (!inherits(object, "Train_Model")) stop("object must be a 'Train_Model' instance.")

  cv_metrics <- object@process.info$cv_metrics
  metric_name <- object@process.info$task_type
  if (is.null(cv_metrics)) {
    res <- slot(object, "all.results")
    if (is.list(res) && !is.data.frame(res) && !is.null(res$cv_selection)) {
      cv_metrics <- res$cv_selection
    }
  }
  if (is.null(cv_metrics) || nrow(cv_metrics) == 0) {
    stop("No CV metrics found on this object. Was it produced by ModelTrainAnalysis()?")
  }

  # (Re-uses the same internal helper defined in model_05_train.R so the
  # figure style -- palette, theme, error bars -- is identical everywhere
  # it's used, whether called from inside ModelTrainAnalysis(save_plots=TRUE)
  # or here, standalone, after the fact.)
  .plot_cv_bar(cv_metrics, metric_selection = "CV metric", palette_name = palette_name,
               save_plots = save_plot, save_dir = save_dir,
               plot_width = plot_width, plot_height = plot_height,
               filename = "cv_comparison_standalone.pdf")
}

#' Confusion Matrix Heatmap for the Best Multiclass Model
#'
#' @param object A `Train_Model` object with `task_type = "multiclass"`.
#' @param dataset One of `"test"` (default) or `"train"`.
#' @param palette_name Viridis option letter, e.g. `"D"`. Default `"D"`.
#' @param save_plot Logical. Save as PDF. Default `FALSE`.
#' @param save_dir Directory to save to. Required if `save_plot = TRUE`.
#' @param plot_width,plot_height Plot dimensions in inches.
#' @return Invisibly, the confusion-matrix table used for the plot.
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_viridis_c labs theme_minimal
#' @export
PlotMultiClassConfusion <- function(object, dataset = c("test", "train"),
                                    palette_name = "D", save_plot = FALSE,
                                    save_dir = NULL, plot_width = 6, plot_height = 5) {
  res <- .assert_nonbinary_result(object)
  if (res$task_type != "multiclass") {
    stop("PlotMultiClassConfusion() only applies to task_type = 'multiclass' (found '",
         res$task_type, "').")
  }
  dataset <- match.arg(dataset)

  best <- object@best.model.result
  data <- if (dataset == "test") object@filtered.set$testing else object@filtered.set$training
  group_col <- object@group_col

  truth <- factor(data[[group_col]])
  pred  <- stats::predict(best$model, newdata = data)
  cm <- caret::confusionMatrix(pred, truth)

  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Predicted", "Actual", "Freq")

  p <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Actual, y = Predicted, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), color = "black") +
    ggplot2::scale_fill_viridis_c(option = palette_name) +
    ggplot2::labs(title = paste("Confusion Matrix --", best$model_type, "(", dataset, "set )"),
                  x = "Actual class", y = "Predicted class") +
    ggplot2::theme_minimal(base_size = 13)

  print(p)

  if (isTRUE(save_plot)) {
    if (is.null(save_dir)) stop("save_dir must be provided when save_plot = TRUE.")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("multiclass_confusion_", dataset, ".pdf")),
                     plot = p, width = plot_width, height = plot_height)
    write.csv(cm_df, file.path(save_dir, paste0("multiclass_confusion_", dataset, ".csv")),
              row.names = FALSE)
  }

  invisible(cm$table)
}

#' One-vs-Rest ROC Curves for the Best Multiclass Model
#'
#' @param object A `Train_Model` object with `task_type = "multiclass"`.
#' @param dataset One of `"test"` (default) or `"train"`.
#' @param palette_name Wes Anderson palette name. Default `"Darjeeling1"`.
#' @param save_plot Logical. Save as PDF. Default `FALSE`.
#' @param save_dir Directory to save to. Required if `save_plot = TRUE`.
#' @param plot_width,plot_height Plot dimensions in inches.
#' @return Invisibly, a list of per-class `pROC::roc` objects and the macro AUC.
#' @importFrom ggplot2 ggplot aes geom_line geom_abline scale_color_manual labs
#' @importFrom ggprism theme_prism
#' @export
PlotMultiClassROC <- function(object, dataset = c("test", "train"),
                              palette_name = "Darjeeling1", save_plot = FALSE,
                              save_dir = NULL, plot_width = 6, plot_height = 5) {
  res <- .assert_nonbinary_result(object)
  if (res$task_type != "multiclass") {
    stop("PlotMultiClassROC() only applies to task_type = 'multiclass' (found '",
         res$task_type, "').")
  }
  dataset <- match.arg(dataset)
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required.")

  best <- object@best.model.result
  data <- if (dataset == "test") object@filtered.set$testing else object@filtered.set$training
  group_col <- object@group_col
  truth <- factor(data[[group_col]])

  probs <- tryCatch(stats::predict(best$model, newdata = data, type = "prob"),
                    error = function(e) NULL)
  if (is.null(probs)) stop("Selected model '", best$model_type, "' does not support predict(..., type = 'prob').")

  classes <- levels(truth)
  roc_list <- list()
  plot_data <- data.frame()
  for (cls in classes) {
    binary_truth <- factor(ifelse(truth == cls, cls, "rest"), levels = c("rest", cls))
    roc_obj <- pROC::roc(binary_truth, probs[[cls]], levels = c("rest", cls),
                         direction = "auto", quiet = TRUE)
    roc_list[[cls]] <- roc_obj
    auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
    plot_data <- rbind(plot_data, data.frame(
      Specificity = 1 - roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      Class = paste0(cls, " (AUC=", auc_val, ")")
    ))
  }

  macro_auc <- tryCatch(as.numeric(pROC::multiclass.roc(truth, probs)$auc), error = function(e) NA_real_)

  n_colors <- length(classes)
  palette_colors <- tryCatch({
    cols <- wesanderson::wes_palette(palette_name, type = "discrete")
    if (length(cols) < n_colors) rep(cols, length.out = n_colors) else cols[1:n_colors]
  }, error = function(e) viridis::viridis(n_colors))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Specificity, y = Sensitivity, color = Class)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = palette_colors) +
    ggplot2::labs(title = paste0("One-vs-Rest ROC -- ", best$model_type,
                                 " (macro AUC = ", round(macro_auc, 3), ")"),
                  x = "1 - Specificity", y = "Sensitivity", color = "Class") +
    ggprism::theme_prism(base_size = 13)

  print(p)

  if (isTRUE(save_plot)) {
    if (is.null(save_dir)) stop("save_dir must be provided when save_plot = TRUE.")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("multiclass_roc_", dataset, ".pdf")),
                     plot = p, width = plot_width, height = plot_height)
    write.csv(plot_data, file.path(save_dir, paste0("multiclass_roc_", dataset, ".csv")),
              row.names = FALSE)
  }

  invisible(list(roc_objects = roc_list, macro_auc = macro_auc))
}

#' Predicted-vs-Actual and Residual Diagnostics for the Best Regression Model
#'
#' @param object A `Train_Model` object with `task_type = "regression"`.
#' @param dataset One of `"test"` (default) or `"train"`.
#' @param save_plot Logical. Save as PDF. Default `FALSE`.
#' @param save_dir Directory to save to. Required if `save_plot = TRUE`.
#' @param plot_width,plot_height Plot dimensions in inches.
#' @return Invisibly, a data.frame with `actual`, `predicted`, `residual`.
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_hline labs
#' @importFrom ggprism theme_prism
#' @importFrom stats predict
#' @importFrom utils write.csv
#' @export
PlotRegressionFit <- function(object, dataset = c("test", "train"),
                              save_plot = FALSE, save_dir = NULL,
                              plot_width = 10, plot_height = 5) {
  
  res <- .assert_nonbinary_result(object)
  if (res$task_type != "regression") {
    stop("PlotRegressionFit() only applies to task_type = 'regression' (found '",
         res$task_type, "').")
  }
  
  dataset <- match.arg(dataset)
  best <- object@best.model.result
  data <- if (dataset == "test") object@filtered.set$testing else object@filtered.set$training
  group_col <- object@group_col
  
  actual <- as.numeric(data[[group_col]])
  predicted <- as.numeric(stats::predict(best$model, newdata = data))
  fit_df <- data.frame(actual = actual, predicted = predicted, residual = actual - predicted)
  
  perf <- if (dataset == "test") best$test_performance else best$train_performance
  
  p1 <- ggplot2::ggplot(fit_df, ggplot2::aes(x = .data$actual, y = .data$predicted)) +
    ggplot2::geom_point(alpha = 0.6, color = "#3B7EA1") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = paste0(best$model_type, " -- Predicted vs Actual (", dataset, ")"),
                  subtitle = paste0("RMSE = ", round(perf$RMSE, 3),
                                    " | R\u00b2 = ", round(perf$Rsquared, 3)),
                  x = "Actual", y = "Predicted") +
    ggprism::theme_prism(base_size = 12)
  
  p2 <- ggplot2::ggplot(fit_df, ggplot2::aes(x = .data$predicted, y = .data$residual)) +
    ggplot2::geom_point(alpha = 0.6, color = "#C0392B") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = "Residuals vs Predicted", x = "Predicted", y = "Residual") +
    ggprism::theme_prism(base_size = 12)
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    p <- gridExtra::grid.arrange(p1, p2, ncol = 2)
  } else {
    print(p1); print(p2)
    p <- list(fit_plot = p1, residual_plot = p2)
  }
  
  if (isTRUE(save_plot)) {
    if (is.null(save_dir)) stop("save_dir must be provided when save_plot = TRUE.")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    ggplot2::ggsave(file.path(save_dir, paste0("regression_fit_", dataset, ".pdf")),
                    plot = p1, width = plot_width / 2, height = plot_height)
    ggplot2::ggsave(file.path(save_dir, paste0("regression_residuals_", dataset, ".pdf")),
                    plot = p2, width = plot_width / 2, height = plot_height)
    utils::write.csv(fit_df, file.path(save_dir, paste0("regression_fit_", dataset, ".csv")),
                     row.names = FALSE)
  }
  
  invisible(fit_df)
}
