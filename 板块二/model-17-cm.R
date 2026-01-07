#' Generate and Plot Confusion Matrix for the Best Model
#'
#' This function evaluates a classification model on a given testing dataset,
#' computes the confusion matrix, and visualizes it using a heatmap-style plot.
#' It optionally saves the plot to a specified directory.
#'
#'
#' @import caret
#' @import stats
#' @param best_model A trained classification model that supports the `predict()` method. Typically an object returned from training (e.g., via `caret`, `glm`, `randomForest`, etc.).
#' @param testing_data A data frame containing testing observations. Must include all predictor variables required by the model and a column for true labels.
#' @param group_col A character string specifying the name of the column in `testing_data` that contains the true class labels. Default is `"group"`.
#' @param save_plots Logical. If `TRUE`, the confusion matrix plot will be saved as a PDF file. Default is `TRUE`.
#' @param save_dir A character string indicating the path to the directory where the plot should be saved. Default is `here("ModelData", "best_cutoff")`.
#' @param plot_width Numeric value specifying the width of the output PDF plot in inches. Default is `6`.
#' @param plot_height Numeric value specifying the height of the output PDF plot in inches. Default is `6`.
#' @param palette_name A character string specifying the color palette name to be used from the `wesanderson` package for the heatmap. Default is `"AsteroidCity1"`.
#' @param seed An integer used to set the random seed for reproducibility. Default is `123`.
#'
#' @returns A list containing:
#' \itemize{
#'   \item `cm_results`: The result of `caret::confusionMatrix()`, containing classification metrics.
#'   \item `cm_plot`: A `ggplot2` object representing the confusion matrix plot.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
#' model <- readRDS("trained_model.rds")
#' test_data <- read.csv("test_data.csv")
#' result <- generate_best_model_confusion_matrix(model, test_data, group_col = "group")
#' print(result$cm_results)
#' }

generate_best_model_confusion_matrix <- function(best_model,
                                                 testing_data,
                                                 group_col = "group",
                                                 save_plots = TRUE,
                                                 save_dir = here("ModelData", "best_cutoff"),
                                                 plot_width = 6,
                                                 plot_height = 6,
                                                 base_size =14,
                                                 palette_name = "AsteroidCity1",
                                                 seed =123,
                                                 best_threshold=0.5,
                                                 set_type = "train") {
  set.seed(seed)
  testing_predictions <- predict(best_model, newdata = testing_data, type = "prob")[, 2]
  
  predicted_labels <- ifelse(testing_predictions > best_threshold, 1, 0)
  

  true_labels <- testing_data[[group_col]]

  true_labels <- factor(true_labels)
  label_check <- function(pred, true) {
    baseline_acc <- mean((true == 0) == (pred == 0))  
    flipped_acc <- mean((true == 0) == (pred == 1))   
    
    if (flipped_acc > baseline_acc) {
      warning(paste("Label direction may be reversed!",
                    "Original accuracy:", round(baseline_acc, 4),
                    "Flipped accuracy:", round(flipped_acc, 4)))
      return(1 - pred)  
    }
    return(pred)
  }
  predicted_labels <- label_check(predicted_labels, true_labels)
  predicted_labels <- factor(predicted_labels)
  true_labels <- factor(true_labels, levels = c(0, 1))
  predicted_labels <- factor(predicted_labels, levels = c(0, 1))
  if (!all(levels(testing_predictions) == levels(true_labels))) {
    stop("The levels of predictions and true labels are not consistent.")
  }

  cm <- caret::confusionMatrix(predicted_labels, true_labels)
  selected_colors <-as.vector(wes_palette(palette_name))

  print(cm)

  cm_plot <- ggplot(as.data.frame(cm$table), aes(Reference, Prediction)) +
    geom_tile(aes(fill = Freq), color = "white") +
    geom_text(aes(label = Freq), vjust = 1, size = 5, color = "black", fontface = "bold") +
    scale_fill_gradient(low = selected_colors[2], high = selected_colors[1]) +
    labs(title = "Confusion Matrix",
         x = "True Label",
         y = "Predicted Label",
         subtitle = paste("Accuracy: ", round(cm$overall['Accuracy'], 4),
                          " | Kappa: ", round(cm$overall['Kappa'], 4))) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95")
    )

  print(cm_plot)

  if (save_plots) {
    ggsave(
      filename = file.path(save_dir, paste0("confusion_matrix_", set_type, ".pdf")),  # 明确添加.pdf扩展名
      plot = cm_plot, 
      width = plot_width, 
      height = plot_height,
      device = "pdf"  
    )
    cat("Plot saved to:", file.path(save_dir, "confusion_matrix_plot.pdf"), "\n")
  }

  return(list(cm_results = cm, cm_plot = cm_plot))
}

#' Generate Confusion Matrix for Best Model on Specified Dataset
#'
#' @param object A Best_Model object containing trained models and data
#' @param group_col Name of the grouping variable column (default="group")
#' @param palette_name Color palette name for plots (default="AsteroidCity1")
#' @param save_plots Logical indicating whether to save plots (default=TRUE)
#' @param save_dir Directory to save plots (default=here("ModelData","best_model_result"))
#' @param plot_width Plot width in inches (default=6)
#' @param plot_height Plot height in inches (default=6)
#' @param best_threshold Decision threshold (default=NULL will use model's default)
#' @param set_type Dataset to use: "train", "test" or "validation" (default="train")
#'
#' @return Updated Best_Model object with confusion matrix results
#' @export
ModelBestCM <- function(object,
                        group_col = "group",
                        palette_name = "AsteroidCity1",
                        save_plots = TRUE,
                        save_dir = here("ModelData", "best_cutoff"),
                        plot_width = 5,
                        plot_height = 5,
                        best_threshold = NULL,
                        set_type = "train") {
  
  if (!inherits(object, "Best_Model")) {
    stop("Input must be an object of class 'Best_Model'.")
  }
  
  cat("Input is of class 'Best_Model'. Extracting datasets...\n")
  
  best_model<-object@best.model[[1]]
  best_model_type<-object@best.model.type
  group_col  <-  object@group_col
  data_sets <- object@filtered.set
  train_data<-data_sets$training
  if (is.null(best_threshold)) {
    best_threshold <- object@best_threshold[["best_threshold"]]
  }
  
  dataset <- switch(tolower(set_type),
                    "train" = data_sets$training,
                    "test" = data_sets$testing,
                    "validation" = data_sets$validation %||% stop("Validation set not found"),
                    "external_validation" = data_sets$external_validation %||% stop("external_validation set not found"),
                    stop("Invalid set_type. Use 'train', 'test', 'validation' or 'external_validation'"))
  dataset <- na.omit(dataset)
  
  dataset <- match_factor_levels(dataset, train_data)
  
  confusion_matrix_results <- generate_best_model_confusion_matrix(
    best_model = best_model,
    testing_data = dataset,
    group_col = group_col,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    palette_name = palette_name,
    best_threshold = best_threshold,
    set_type = set_type)
  
  slot_name <- paste0("confusion_matrix_", set_type)
  object@cm.result[["confusion_matrix_results"]][[slot_name]] <- confusion_matrix_results
  
  cat("\n=== Confusion Matrix Summary ===\n")
  cat("Threshold applied:", best_threshold, "\n")
  cat("Updating 'Best_Model' object...\n")
  cat("The 'Best_Model' object has been updated with the following slots:\n")
  cat("- 'cm.result' slot updated with confusion matrix results.\n")
  
  return(object)
}


#' Match Factor Levels Between External and Reference Datasets
#'
#' This function ensures that factor variables in an external dataset have the same levels 
#' as those in a reference dataset (e.g., training set). This is particularly useful 
#' before applying a trained model to new data to avoid factor level mismatches.
#'
#' @param external_data A data frame containing new data (e.g., test or validation set) 
#'   to be aligned in terms of factor levels.
#' @param reference_data A data frame used as the reference for factor levels 
#'   (usually the training dataset).
#'
#' @return A modified version of \code{external_data} where factor variables 
#'   that also appear in \code{reference_data} are converted to factors 
#'   with matched levels.
#'
#' @examples
#' # Assume train_data and test_data share some columns
#' test_data <- match_factor_levels(test_data, train_data)
#'
#' @export
#' 
match_factor_levels <- function(external_data, reference_data) {
  common_vars <- intersect(names(external_data), names(reference_data))
    for (var in common_vars) {
    if (is.factor(reference_data[[var]])) {
      external_data[[var]] <- factor(external_data[[var]], 
                                     levels = levels(reference_data[[var]]))
    }
  }
  
  return(external_data)
}
