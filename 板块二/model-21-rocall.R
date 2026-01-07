#' Generate and Save ROC Plots for the Best Model in a 'Best_Model' Object
#'
#' This function extracts the best model and relevant datasets from a 'Best_Model' object,
#' computes ROC curves across different validation cohorts (training, testing, validation, external validation),
#' and visualizes the ROC curves. It also updates the 'Best_Model' object with the resulting ROC plot data.
#'
#' @param object An object of class 'Best_Model' that contains the best model and associated datasets.
#' @param group_col The name of the grouping/label column. Default is "group". If the object is of class 'Best_Model', this will be overwritten.
#' @param palette_name The name of the color palette to use from the wesanderson package. Default is "AsteroidCity1".
#' @param base_size The base font size for plots. Default is 14.
#' @param save_plots Logical; whether to save the ROC plot as a PDF. Default is TRUE.
#' @param save_dir Directory where the ROC plot will be saved. Default is `here("ModelData", "best_model_result")`.
#' @param plot_width Width of the saved plot in inches. Default is 5.
#' @param plot_height Height of the saved plot in inches. Default is 5.
#' @param alpha Transparency level for plot elements. Default is 1 (opaque).
#'
#' @return An updated 'Best_Model' object with ROC plot data stored in the `best.model.result[["roc_plots"]]` slot.
#' @export
#'
#' @examples
#' \dontrun{
#'   object_best <- ModelBestRoc(object = object_best)
#' }

ModelBestRoc <- function(object,
                         group_col = "group",
                         palette_name = "AsteroidCity1",
                         base_size = 14,
                         save_plots = TRUE,
                         save_dir = here("ModelData", "best_model_result"),
                         plot_width = 5,
                         plot_height = 5,
                         alpha = 1) {

  if (inherits(object, "Best_Model")) {
    cat("Input is of class 'Best_Model'. Extracting datasets...\n")
    best_model<-object@best.model[[1]]
    best_model_type<-object@best.model.type
    group_col  <-  object@group_col

    train_data <- object@filtered.set[["training"]]
    test_data <- object@filtered.set[["testing"]]
    validation_data <- if ("validation" %in% names(object@filtered.set)) object@filtered.set$validation else NULL
    external_validation <- if ("external_validation" %in% names(object@filtered.set)) object@filtered.set$external_validation else NULL
  }

  cat("Plotting ROC curve for the best model...\n")
  set.seed(1234)
  roc_results <- plot_best_model_roc(
    best_model = best_model,
    train_data = train_data,
    test_data = test_data,
    validation_data = validation_data,
    external_validation = external_validation,
    group_col = group_col,
    palette_name = palette_name,
    base_size = base_size,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    alpha = alpha
  )

  if (inherits(object, "Best_Model")) {
    object@performance.result[["roc_plots"]] <- roc_results
    cat("Updating 'Best_Model' object...\n")
    cat("The 'Best_Model' object has been updated with the following slots:\n")
    cat("- 'best.model.result' slot updated.\n")

    return(object)
  } else {
    stop("Input must be an object of class 'Best_Model'")
  }
}
#' Plot ROC Curves for the Best Model and Compare AUCs of Multiple Datasets
#'
#' This function generates and visualizes ROC curves for a given trained classification model
#' across multiple datasets (training, testing, validation, and external validation). AUC values
#' and 95% confidence intervals are computed for each dataset. The ROC curves are optionally saved as a PDF.
#'
#' @import stats
#' @import ggplot2
#' @import here
#'
#' @param best_model The best trained classification model.
#' @param train_data (Optional) A data frame of the training set.
#' @param test_data (Optional) A data frame of the testing set.
#' @param validation_data (Optional) A data frame of the internal validation set.
#' @param external_validation (Optional) A data frame of the external validation set.
#' @param group_col Character; the name of the column containing binary class labels (e.g., "0" and "1"). Default is "group".
#' @param palette_name Character; the name of the color palette from the wesanderson package. Default is "AsteroidCity1".
#' @param base_size Numeric; base font size for plot text. Default is 14.
#' @param save_plots Logical; whether to save the ROC plot to file. Default is TRUE.
#' @param save_dir Character; directory to save the ROC plot. Default is `here("ModelData", "best_model_result")`.
#' @param plot_width Numeric; width of the saved plot in inches. Default is 5.
#' @param plot_height Numeric; height of the saved plot in inches. Default is 5.
#' @param alpha Numeric; transparency of ROC lines. Default is 1.
#' @param subtitle Character; subtitle text for the plot. Default is "Training and Testing Datasets".
#'
#' @return A named list of data frames, each containing ROC curve plotting data (specificity, sensitivity, and dataset label).
#' @export
#'
#' @examples
#' \dontrun{
#'   plot_best_model_roc(best_model = trained_model,
#'                       train_data = train_df,
#'                       test_data = test_df)
#' }

plot_best_model_roc <- function(best_model,
                                train_data ,
                                test_data ,
                                validation_data ,
                                external_validation,
                                group_col = "group",
                                palette_name = "AsteroidCity1",
                                base_size = 14,
                                save_plots = TRUE,
                                save_dir = here("ModelData", "best_model_result"),
                                plot_width = 5,
                                plot_height = 5,
                                alpha = 1,
                                subtitle = "Training and Testing Datasets") {

  plot_data_list <- list()

  if (!is.null(train_data)) {
    training_predictions <- predict(best_model, newdata = train_data, type = "prob")[, 2]
    roc_training <- roc(train_data[[group_col]], training_predictions, levels = c("0", "1"), direction = "<")

    auc_training <- auc(roc_training)
    auc_ci_training <- ci.auc(roc_training)
    training_plot_data <- data.frame(
      Specificity = 1 - roc_training$specificities,
      Sensitivity = roc_training$sensitivities,
      Dataset = paste0("Training Set (AUC = ", sprintf("%.3f", auc_training),
                       " ± ", sprintf("%.3f", (auc_ci_training[3]-auc_ci_training[1])/2), ")")
    )

    plot_data_list$training <- training_plot_data
  }

  if (!is.null(test_data)) {
    testing_predictions <- predict(best_model, newdata = test_data, type = "prob")[, 2]
    roc_testing <- roc(test_data[[group_col]], testing_predictions, levels = c("0", "1"), direction = "<")

    auc_testing <- auc(roc_testing)
    auc_ci_testing <- ci.auc(roc_testing)
    testing_plot_data <- data.frame(
      Specificity = 1 - roc_testing$specificities,
      Sensitivity = roc_testing$sensitivities,
      Dataset = paste0("Testing Set (AUC = ", sprintf("%.3f", auc_testing),
                       " ± ", sprintf("%.3f", (auc_ci_testing[3]-auc_ci_testing[1])/2), ")")
    )

    plot_data_list$testing <- testing_plot_data
  }



  if (!is.null(validation_data)) {
    validation_data <- match_factor_levels(validation_data, train_data)
    validation_data[[group_col]] <- factor(validation_data[[group_col]], levels = c("0", "1"))
    validation_predictions <- predict(best_model, newdata = validation_data, type = "prob")[, 2]
    roc_validation <- roc(validation_data[[group_col]], validation_predictions, levels = c("0", "1"), direction = ">")

    if (auc(roc_validation) < 0.5) {
      cat("Warning: Model predictions are inverted. Reversing prediction probabilities.\n")
      validation_predictions <- 1 - validation_predictions
      roc_validation <- roc(validation_data[[group_col]], validation_predictions, levels = c("0", "1"), direction = ">")
    }

    auc_validation <- auc(roc_validation)
    auc_ci_validation <- ci.auc(roc_validation)
    validation_plot_data <- data.frame(
      Specificity = 1 - roc_validation$specificities,
      Sensitivity = roc_validation$sensitivities,
      Dataset = paste0("Validation Set (AUC = ", sprintf("%.3f", auc_validation),
                       " ± ", sprintf("%.3f", (auc_ci_validation[3] - auc_ci_validation[1]) / 2), ")")
    )
    plot_data_list$validation <- validation_plot_data
  }

  if (!is.null(external_validation)) {
    external_validation <- na.omit(external_validation)

    external_validation <- match_factor_levels(external_validation, train_data)

    external_validation[[group_col]] <- factor(external_validation[[group_col]], levels = c("0", "1"))
    external_validation_predictions <- predict(best_model, newdata = external_validation, type = "prob")[, 2]
    roc_external_validation <- roc(external_validation[[group_col]], external_validation_predictions, levels = c("0", "1"), direction = ">")

    if (auc(roc_external_validation) < 0.5) {
      cat("Warning: Model predictions are inverted. Reversing prediction probabilities.\n")
      external_validation_predictions <- 1 - external_validation_predictions
      roc_external_validation <- roc(external_validation[[group_col]], external_validation_predictions, levels = c("0", "1"), direction = ">")
    }

    auc_external_validation <- auc(roc_external_validation)
    auc_ci_external_validation <- ci.auc(roc_external_validation)
    external_validation_plot_data <- data.frame(
      Specificity = 1 - roc_external_validation$specificities,
      Sensitivity = roc_external_validation$sensitivities,
      Dataset = paste0("External Validation Set (AUC = ", sprintf("%.3f", auc_external_validation),
                       " ± ", sprintf("%.3f", (auc_ci_external_validation[3] - auc_ci_external_validation[1]) / 2), ")")
    )
    plot_data_list$external_validation <- external_validation_plot_data
  }


  combined_plot_data <- do.call(rbind, plot_data_list)
  combined_plot_data$Dataset <- factor(combined_plot_data$Dataset,
                                       levels = unlist(lapply(plot_data_list, function(x) unique(x$Dataset))))


  palette_colors <- wes_palette(name = palette_name, n = length(unique(combined_plot_data$Dataset)), type = "discrete")

  p <- ggplot(combined_plot_data, aes(x = Specificity, y = Sensitivity, color = Dataset)) +
    geom_line(size = 1.25, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = palette_colors) +
    labs(
      title = "ROC Curves Validation Comparison",
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Validation Cohort"
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      expand = expansion(mult = 0.01)
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = scales::alpha("white", 0.8)),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  print(p)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(save_dir, "best_model_roc_plot.pdf"),
           plot = p,
           width = plot_width,
           height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "best_model_roc_plot.pdf"), "\n")
  }

  return(plot_data_list)
}
