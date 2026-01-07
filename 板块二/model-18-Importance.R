#' Extract and Visualize Feature Importance from a Model
#'
#' This function calculates feature importance from a trained model, visualizes the top N features,
#' and saves both the plot and importance data to files.
#'
#' @param best_model A trained model object that works with `varImp()` from the caret package
#' @param top_n Number of top features to display (default: 15). If NULL, shows all features.
#' @param palette_name Name of color palette from wesanderson package (default: "AsteroidCity1")
#' @param save_plots Logical indicating whether to save output files (default: TRUE)
#' @param save_dir Directory path where outputs should be saved (default: here::here("ModelData", "best_model_result"))
#' @param plot_width Width of output plot in inches (default: 5)
#' @param plot_height Height of output plot in inches (default: 5)
#' @param base_size Base font size for plot text (default: 14)
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @return A list containing:
#' \itemize{
#'   \item{plot: ggplot object of the feature importance plot}
#'   \item{importance_df: Data frame of all features sorted by importance}
#' }
#' When save_plots=TRUE, also saves:
#' \itemize{
#'   \item{Feature_Importance.pdf: Bar plot of top features}
#'   \item{feature_importance.csv: CSV file with all feature importance values}
#' }
#'
#' @examples
#' \dontrun{
#' # After training a model with caret
#' model <- train(Species ~ ., data = iris, method = "rf")
#'
#' # Get and save feature importance
#' results <- get_feature_importance(
#'   best_model = model,
#'   top_n = 10,
#'   save_dir = "model_results"
#' )
#'
#' # Access results
#' importance_plot <- results$plot
#' importance_data <- results$importance_df
#' }
#'
#' @export
get_feature_importance <- function(
    best_model,
    top_n = 15,
    palette_name = "AsteroidCity1",
    save_plots = TRUE,
    save_dir = here::here("ModelData", "best_model_result"),
    plot_width = 5,
    plot_height = 5,
    base_size = 14,
    seed = 123,
    save_data = TRUE,
    csv_filename = "feature_importance.csv"
) {
  set.seed(seed)

  feature_importance <- tryCatch(
    varImp(best_model),
    error = function(e) stop("Failed to calculate feature importance: ", e$message)
  )

  importance_df <- as.data.frame(feature_importance$importance)

  if (is.null(rownames(importance_df))) {
    stop("Row names are missing in the importance dataframe. Ensure that feature names are set.")
  }

  importance_df <- data.frame(
    Feature = rownames(feature_importance$importance),
    Importance = importance_df[, 1],
    row.names = NULL
  )

  sorted_importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

  if (is.null(top_n)) {
    top_n <- nrow(sorted_importance_df)
  } else if (top_n > nrow(sorted_importance_df)) {
    warning("top_n is greater than the number of features, displaying all features")
    top_n <- nrow(sorted_importance_df)
  }

  # Get top features
  top_features <- head(sorted_importance_df, top_n)

  bar_colors <- tryCatch(
    as.vector(wes_palette(palette_name)),
    error = function(e) {
      warning("Failed to use specified palette, using default colors")
      scales::hue_pal()(top_n)
    }
  )

  p <- ggplot(
    top_features,
    aes(x = reorder(Feature, Importance), y = Importance)
  ) +
    geom_bar(
      stat = "identity",
      fill = rep(bar_colors, length.out = top_n)
    ) +
    labs(
      x = "Features",
      y = "Importance",
      title = paste("Top", top_n, "Feature Importance")
    ) +
    ggprism::theme_prism(base_size = base_size) +
    coord_flip()

  print(p)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("\nCreated directory: ", save_dir)
    }

    plot_path <- file.path(save_dir, "Feature_Importance.pdf")
    ggsave(
      filename = plot_path,
      plot = p,
      width = plot_width,
      height = plot_height,
      device = "pdf"
    )
    cat("\nPlot saved to: ", plot_path)
  }
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("Created directory:", save_dir, "\n")
    }

    data_path <- file.path(save_dir, csv_filename)
    write.csv(sorted_importance_df, data_path, row.names = FALSE)

  }
  result <- list(
    plot = p,
    importance_df = sorted_importance_df
  )

  return(result)
}

#' Generate and Visualize Feature Importance for the Best Model
#'
#' This function calculates and visualizes the feature importance of the best model
#' stored within a `Best_Model` object or from a model object directly. It supports
#' saving the plot and storing the result in the object for further use.
#'
#' @import methods
#' @import stats
#' @import here
#' @import caret
#' @param object An object of class `Best_Model`, or a trained model object directly (e.g., from caret or randomForest).
#' @param top_n Integer specifying the number of top features to display. Default is 15.
#' @param palette_name Character string specifying the color palette name from `wesanderson`. Default is `"AsteroidCity1"`.
#' @param save_plots Logical indicating whether to save the plot to file. Default is `TRUE`.
#' @param save_dir Directory path where plots will be saved. Default is `here("ModelData", "best_model_result")`.
#' @param plot_width Width of the saved plot in inches. Default is 5.
#' @param plot_height Height of the saved plot in inches. Default is 5.
#' @param base_size Base font size for the plot theme. Default is 14.
#'
#' @returns If `object` is of class `Best_Model`, returns the updated object with a feature importance plot stored in the `@best.model.result$top_features_plot` slot.
#' If a model object is provided directly, returns a `ggplot` object of the importance plot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # If `model_data` is a Best_Model object
#' object_model <- FeatureImportance(object_model)}

FeatureImportance <- function(object,
                              top_n = 15,
                              palette_name = "AsteroidCity1",
                              save_plots = TRUE,
                              save_dir =  here("ModelData", "best_model_result"),
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 14) {
  if (is.null(object)) {
    stop("Invalid input: 'object' should be provided")
  }

  if (inherits(object, "Best_Model")) {
    cat("\nInput is of class 'Best_Model'. Extracting the best model...\n")

    best_model<-object@best.model[[1]]
    best_model_type<-object@best.model.type
    group_col  <-  object@group_col
    if (is.null(best_model)) {
      stop("Invalid input: 'object' should contain a valid model")
    }

    cat("\nCalculating feature importance for the best model extracted from 'Best_Model'...\n")
    result <- get_feature_importance(best_model=best_model,
                                     top_n = top_n,
                                     palette_name = palette_name,
                                     save_plots = save_plots,
                                     save_dir =  save_dir,
                                     plot_width = plot_width,
                                     plot_height = plot_height,
                                     base_size = base_size)


    object@feature.importance <- result
    cat("\nUpdating 'Best_Model' object...\n")
    cat("The 'Best_Model' object has been updated with the following slots:\n")
    cat("- 'feature.importance' slot updated.\n")
    return(object)

  } else {
    cat("Object is provided directly. Calculating feature importance...\n")
    result <- FeatureImportance(object, top_n = top_n)
    cat("Feature importance calculated for the provided model.\n")
    return(result)
  }
}
