#' Enhanced Model Comparison Function with Improved ROC Visualization
#'
#' @param model_list List of model objects (named list)
#' @param model_names Vector of model names (optional)
#' @param palette_name Color palette name (default: "Set1")
#' @param base_size Base font size for plots (default: 14)
#' @param save_plots Whether to save plots (default: FALSE)
#' @param save_dir Directory to save plots (default: current directory)
#' @param plot_width Plot width in inches (default: 7)
#' @param plot_height Plot height in inches (default: 5)
#' @param alpha Line transparency (default: 1)
#' @return List containing performance metrics and enhanced ROC curve visualization
#' @importFrom pROC roc auc ci.auc
#' @import caret
#' @export
compare_models <- function(model_list,
                           model_names = NULL,
                           palette_name = "AsteroidCity1",
                           base_size = 14,
                           save_plots = TRUE,
                           save_dir =  here("ModelData", "model_comparison"),
                           plot_width = 5,
                           plot_height = 5,
                           alpha = 1) {

  cat("Starting model comparison analysis...\n")

  if (!is.list(model_list)) stop("model_list must be a list")
  if (is.null(model_names)) {
    model_names <- names(model_list)
    if (is.null(model_names)) {
      model_names <- paste0("Model_", seq_along(model_list))
      cat("  - Generated default model names\n")
    }
  }

  results <- list()
  roc_list <- list()
  plot_data_list <- list()
  auc_results <- numeric()

  cat("- Processing", length(model_list), "models...\n")
  model_data <- lapply(seq_along(model_list), function(i) {
    model <- model_list[[i]]
    current_model_name <- model_names[i]

    cat("  > Analyzing model:", current_model_name, "\n")

    re <- ExtractModel(model)
    best_model <- re$best_model
    train_data <- re$filtered.set[["training"]]
    group_col <- re$group_col

    levels <- validate_binary_classification(train_data, group_col)

    pred_prob <- predict(best_model, newdata = train_data, type = "prob")[,2]

    roc_obj <- pROC::roc(train_data[[group_col]], pred_prob, levels = c("0", "1"), direction = "<")
    auc_value <- pROC::auc(roc_obj)
    auc_ci <- pROC::ci.auc(roc_obj)

    roc_list[[current_model_name]] <<- roc_obj
    auc_results[current_model_name] <<- auc_value

    plot_data_list[[current_model_name]] <<- data.frame(
      Specificity = 1 - roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      Model = paste0(current_model_name, " (AUC = ", round(auc_value, 3),
                     ", CI = [", round(auc_ci[1], 3), ", ", round(auc_ci[3], 3), "])")
    )

    pred_factor <- factor(ifelse(pred_prob > 0.5, 1, 0), levels = c(0, 1))
    true_factor <- factor(train_data[[group_col]], levels = c(0, 1))

    list(
      name = current_model_name,
      model = best_model,
      roc = roc_obj,
      metrics = list(
        AUC = auc_value,
        AUC_CI_lower = auc_ci[1],
        AUC_CI_upper = auc_ci[3],
        Accuracy = mean(true_factor == pred_factor),
        Sensitivity = caret::sensitivity(pred_factor, true_factor, positive = "1"),
        Specificity = caret::specificity(pred_factor, true_factor, negative = "0"),
        F1 = caret::F_meas(pred_factor, true_factor, relevant = "1")
      ),
      train_data = train_data,
      group_col = group_col
    )
  })

  cat("- Compiling performance metrics...\n")
  performance_table <- do.call(rbind, lapply(model_data, function(x) {
    data.frame(
      Model = x$name,
      as.data.frame(x$metrics)
    )
  }))

  combined_plot_data <- do.call(rbind, plot_data_list)

  roc_plot <- ggplot(combined_plot_data, aes(x = Specificity, y = Sensitivity, color = Model)) +
    geom_line(size = 1.25, alpha = alpha) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    scale_color_manual(values = wes_palette(palette_name))  +
    labs(title = "ROC Curves Comparison",
         subtitle = "Including AUC and 95% Confidence Intervals",
         x = "1 - Specificity",
         y = "Sensitivity",
         color = "Model (AUC and CI)") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.8)),
      legend.title = element_text(face = "bold", size = base_size * 0.7),
      legend.text = element_text(size = base_size * 0.6),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40")
    )
  print(roc_plot)
  if (save_plots) {
    cat("- Saving output plots...\n")
    if (!dir.exists(save_dir)) {
      cat("  - Creating output directory:", save_dir, "\n")
      dir.create(save_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(save_dir, "model_comparison_roc.pdf"),
           plot = roc_plot,
           width = plot_width,
           height = plot_height,
           device = "pdf")
    cat("Plot saved to: ", file.path(save_dir, "model_comparison_roc.pdf"))
  }

  cat("- Model comparison completed successfully!\n")

  list(
    performance = performance_table,
    roc_plot = roc_plot,
    auc_values = auc_results,
    roc_objects = roc_list,
    model_data = model_data
  )
}


#' Validate Binary Classification Data
#'
#' Checks and validates that a specified column in a dataframe contains proper binary classification labels.
#' Automatically converts to factor if needed and verifies exactly two levels exist.
#'
#' @param data Input dataframe containing the classification labels
#' @param group_col Name of the column containing classification labels (character)
#'
#' @return Vector containing the two factor levels
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data <- data.frame(outcome = c(0, 1, 0, 1))
#' levels <- validate_binary_classification(data, "outcome")
#' }
#'
#' @export
validate_binary_classification <- function(data, group_col) {
  if (!is.factor(data[[group_col]])) {
    data[[group_col]] <- factor(data[[group_col]])
    warning(paste("Converted", group_col, "to factor"))
  }

  group_levels <- levels(data[[group_col]])
  if (length(group_levels) != 2) {
    stop(paste(group_col, "must have exactly 2 levels. Found:",
               paste(group_levels, collapse = ", ")))
  }
  return(group_levels)
}


#' Extract Model Components from Train_Model Object
#'
#' Extracts key components (filtered dataset, group column, best model, and associated metrics)
#' from a `Train_Model` S4 object. Validates input and returns an S3 object of class `ExtractedModel`.
#'
#' @param object A `Train_Model` S4 object containing modeling results. Must contain at minimum
#'   `filtered.set`, `group_col`, and `best.model.result` slots.
#'
#' @return An S3 object of class `ExtractedModel` containing:
#' \itemize{
#'   \item \code{filtered.set} - The filtered dataset from the input object
#'   \item \code{group_col} - The name of the group/outcome column
#'   \item \code{best_model} - The best performing model (NULL if not found)
#'   \item \code{model_type} - Type/name of the best model
#'   \item \code{model_metrics} - Performance metrics (if available)
#'   \item \code{feature_importance} - Feature importance scores (if available)
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'trained_model' is a valid Train_Model object
#' extracted <- ExtractModel(trained_model)
#' print(extracted$best_model)
#' }
#'
#' @export
#' @seealso \code{\link{validate_binary_classification}} for validating the extracted dataset
ExtractModel <- function(object) {
  if (is.null(object)) {
    stop("Invalid input: 'object' should be provided")
  }

  if (!inherits(object, "Train_Model")) {
    stop("Input object must be of class 'Train_Model'")
  }
  result <- list()

  result$filtered.set <- object@filtered.set
  result$group_col <- object@group_col

  if (length(object@best.model.result) > 0) {
    result$model_type <- object@best.model.result[["model_type"]]
    result$best_model <- object@best.model.result[["model"]][[1]]  # Using [[1]] to extract the first element

    if ("metrics" %in% names(object@best.model.result)) {
      result$model_metrics <- object@best.model.result[["metrics"]]
    }

    if ("feature_importance" %in% names(object@best.model.result)) {
      result$feature_importance <- object@best.model.result[["feature_importance"]]
    }

  } else {
    warning("No best model result found in the Train_Model object")
    result$best_model <- NULL
    result$model_type <- NULL
  }

  class(result) <- "ExtractedModel"

  return(result)
}
