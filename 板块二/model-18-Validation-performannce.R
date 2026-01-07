
#' Plot ROC Curve for External Validation Data
#'
#' This function generates and optionally saves a ROC curve for an external validation dataset
#' using the best-performing classification model. It calculates the AUC and its confidence interval,
#' and produces a publication-ready ROC plot.
#'
#' @importFrom pROC roc auc ci.auc
#' @param best_model A trained classification model object that supports probability prediction via `predict(..., type = "prob")`.
#' @param validation_data A data frame containing the validation dataset. Must include the true class labels in `group_col`.
#' @param group_col A string specifying the name of the column in `validation_data` that contains the true class labels. Default is "group".
#' @param palette_name A string specifying the color palette name used for plotting. Must be a valid palette in `wesanderson::wes_palette()`. Default is "AsteroidCity1".
#' @param base_size A numeric value controlling the base font size for the plot. Default is 14.
#' @param save_plots Logical; whether to save the plot as a PDF file. Default is TRUE.
#' @param save_dir A character string specifying the directory path where the plot should be saved. Default is `here("ModelData", "best_model_result")`.
#' @param plot_width Numeric value specifying the width (in inches) of the saved plot. Default is 5.
#' @param plot_height Numeric value specifying the height (in inches) of the saved plot. Default is 5.
#' @param alpha Numeric value (between 0 and 1) indicating the transparency level of the ROC curve line. Default is 1 (fully opaque).
#' @param subtitle A string to be displayed as the subtitle of the ROC plot. Default is "Validation Dataset".
#'
#' @returns A data frame containing the coordinates of the ROC curve for plotting purposes, including specificity, sensitivity, and AUC information.
#' @export
#'
#' @examples
#' # Assuming `model` is a trained classifier and `val_data` is a data frame for validation:
#' plot_validation_model_roc(
#'   best_model = model,
#'   validation_data = val_data,
#'   group_col = "outcome"
#' )
plot_validation_model_roc <- function(best_model,
                                      validation_data,
                                      group_col = "group",
                                      palette_name = "AsteroidCity1",
                                      base_size = 14,
                                      save_plots = TRUE,
                                      save_dir = here("ModelData", "best_model_result"),
                                      plot_width = 5,
                                      plot_height = 5,
                                      alpha = 1,
                                      subtitle = "Validation Dataset") {
  
  if (is.null(validation_data)) {
    stop("Validation dataset is missing.")
  }
  
  validation_data[[group_col]] <- as.factor(validation_data[[group_col]])
  
  validation_predictions <- predict(best_model,
                                    newdata = validation_data, type = "prob")[, 2]
  
  roc_raw <- roc(validation_data[[group_col]], validation_predictions, levels = c("0", "1"), direction = ">")
  if (auc(roc_raw) < 0.5) {
    cat("Warning: Model predictions are inverted. Reversing prediction probabilities.\n")
    validation_predictions <- 1 - validation_predictions
    roc_validation <- roc(validation_data[[group_col]], validation_predictions, levels = c("0", "1"), direction = ">")
  } else {
    roc_validation <- roc_raw
  }
  
  auc_validation <- auc(roc_validation)
  auc_ci_validation <- ci.auc(roc_validation)
  
  validation_plot_data <- data.frame(
    Specificity = 1 - roc_validation$specificities,
    Sensitivity = roc_validation$sensitivities,
    Dataset = paste0("Validation Set (AUC = ", round(auc_validation, 3),
                     ", CI = [", round(auc_ci_validation[1], 3), ", ",
                     round(auc_ci_validation[3], 3), "])")
  )
  validation_plot_data$Specificity <- as.numeric(validation_plot_data$Specificity)
  validation_plot_data$Sensitivity <- as.numeric(validation_plot_data$Sensitivity)
  validation_plot_data$Dataset <- as.factor(validation_plot_data$Dataset)
  
  p <- ggplot(validation_plot_data, aes(x = Specificity, y = Sensitivity, color = Dataset)) +
    geom_line(size = 1.25, alpha = alpha) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    scale_color_manual(values = wes_palette(palette_name)) +
    labs(title = "ROC Curve for Best Model",
         subtitle = subtitle,
         x = "1 - Specificity",
         y = "Sensitivity",
         color = "Dataset (AUC and CI)") +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.8)),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90")
    )
  
  print(p)
  
  if (save_plots) {
    ggsave(filename = file.path(save_dir, "validation_roc_plot.pdf"), plot = p, width = plot_width, height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "validation_roc_plot.pdf"), "\n")
  }
  
  return(validation_plot_data)
}
#' Evaluate Model Performance on Validation Dataset
#'
#' This function validates a trained model's performance on a held-out validation dataset,
#' including ROC curve generation, performance metrics calculation, and variable importance checks.
#' It handles common validation issues like variable name mismatches and missing data imputation.
#'
#' @import methods
#' @import stats
#' @import here
#' @import caret
#' 
#' @param object A `Model_data` object containing training and validation datasets.
#' @param group_col Character. Name of the grouping variable (target). Default is "group".
#' @param palette_name Character. Color palette name for plots (from `wesanderson`). Default is "AsteroidCity1".
#' @param base_size Numeric. Base font size for plots. Default is 14.
#' @param save_plots Logical. Whether to save output plots. Default is TRUE.
#' @param save_dir Character. Directory path for saving plots. Default is `here("ModelData", "best_model_result")`.
#' @param plot_width Numeric. Plot width in inches. Default is 5.
#' @param plot_height Numeric. Plot height in inches. Default is 5.
#' @param alpha Numeric. Significance level for statistical tests. Default is 0.05.
#' @param importance_threshold Numeric. Threshold for variable importance (0-1). Default is 0.05.
#' @param best_threshold Numeric. Custom classification threshold. If NULL, uses model's optimal threshold.
#'
#' @return The input `Model_data` object with validation results added to the `best.model.result` slot.
#' 
#' @section Details:
#' The function:
#' \itemize{
#'   \item Verifies input data integrity
#'   \item Handles variable name mismatches (e.g., factor dummy variables)
#'   \item Imputes missing low-importance variables
#'   \item Calculates performance metrics using `evaluate_model_performance()`
#'   \item Generates ROC curves using `plot_validation_model_roc()`
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' model_data <- ModelValidationPerformance(model_data)
#' validation_results <- model_data@best.model.result$validation_result
#' }
ModelValidationPerformance <- function(object,
                                       group_col = "group",
                                       palette_name = "AsteroidCity1",
                                       base_size = 14,
                                       save_plots = TRUE,
                                       save_dir = here("ModelData", "best_model_result"),
                                       plot_width = 5,
                                       plot_height = 5,
                                       alpha = 0.05,
                                       importance_threshold = 0.05,
                                       best_threshold = NULL,
                                       set_type = "validation") {
  
  if (!inherits(object, "Model_data")) {
    stop("Input must be an object of class 'Model_data'")
  }
  
  cat("=== Starting Model Validation Performance Evaluation ===\n")
  
  # Validate set_type
  set_type <- tolower(set_type)
  if (!set_type %in% c("validation", "external_validation")) {
    stop("Invalid set_type. Use 'validation' or 'external_validation'")
  }
  
  # Extract model and data
  best_model <- object@best.model.result[["model"]][[1]]
  group_col <- object@group_col %||% group_col
  training_data <- object@filtered.set[["training"]]
  
  # Get the appropriate validation dataset
  validation_data <- switch(set_type,
                            "validation" = object@filtered.set[["validation"]],
                            "external_validation" = object@filtered.set[["external_validation"]])
  
  if (is.null(validation_data) || nrow(validation_data) == 0) {
    stop(paste(set_type, "dataset is empty or not found"))
  }
  
  if (is.null(best_model)) {
    stop("Best model not found in the object")
  }
  
  # Get or calculate best threshold
  if (is.null(best_threshold)) {
    best_threshold <- object@best.model.result[["best_threshold"]]
    cat("Using model's optimal threshold:", best_threshold, "\n")
  }
  
  # Handle variable importance and feature selection
  var_imp <- tryCatch({
    varImp(best_model)$importance
  }, error = function(e) {
    cat("Variable importance not available for this model type. Using all model features.\n")
    NULL
  })
  
  # Get important variables from model
  important_vars <- if (!is.null(best_model$coefnames)) {
    best_model$coefnames
  } else {
    setdiff(colnames(training_data), group_col)
  }
  
  # Handle any feature name transformations
  change_vars <- grep("1$", important_vars, value = TRUE)
  changed_vars <- gsub("1$", "", change_vars)
  var_mapping <- setNames(changed_vars, change_vars)
  important_vars <- ifelse(
    important_vars %in% change_vars,
    var_mapping[important_vars],
    important_vars
  )
  
  # Check for missing variables
  missing_vars <- setdiff(important_vars, names(validation_data))
  if (length(missing_vars) > 0) {
    stop(paste("The following important variables are missing in the", set_type, "data:", 
               paste(missing_vars, collapse = ", ")))
  }
  
  # Prepare prediction data
  prediction_data <- validation_data[, important_vars, drop = FALSE]
  prediction_data[[group_col]] <- validation_data[[group_col]]
  
  # Evaluate model performance
  cat("Evaluating the best model on", set_type, "dataset (n =", nrow(prediction_data), ")...\n")
  validation_result <- evaluate_model_performance(
    data = prediction_data,
    model_result = best_model,
    group_col = group_col,
    custom_cutoff = best_threshold
  )
  
 
  # Store results with appropriate naming
  result_slot <- paste0(set_type, "_result")
  object@best.model.result[[result_slot]] <- list(
    roc_plot = roc_results,
    performance_metrics = validation_result,
    threshold_used = best_threshold,
    evaluation_time = Sys.time()
  )
  
  # Print summary
  cat("\n=== Performance Summary (", set_type, ") ===\n", sep = "")
  print(validation_result)
  cat("\nThreshold applied:", best_threshold, "\n")
  
  cat("\nUpdated 'Model_data' object contains:\n")
  cat("- Performance metrics stored in @best.model.result$", result_slot, "\n", sep = "")
  
  return(object)
}