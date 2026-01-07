#' Clinical Prediction Using a Trained Model
#'
#' Applies a trained model from a Model_data object to new clinical data, generating 
#' predictions and diagnostic visualizations. Handles variable name alignment and 
#' missing data checks.
#'
#' @details Key features:
#' \itemize{
#'   \item Automatically matches variable names between training and prediction data
#'   \item Validates presence of all required predictor variables
#'   \item Generates boxplot-jitter visualization of prediction probabilities
#'   \item Saves PDF visualization to specified directory
#' }
#'
#' @param object A Model_data object containing:
#' \itemize{
#'   \item `@best.model.result$model`: The trained model
#'   \item `@filtered.set$training`: Training data reference
#'   \item `@best.model.result$best_threshold`: (Optional) Default classification threshold
#' }
#' @param new_data Data frame containing samples for prediction. Must include:
#' \itemize{
#'   \item All variables used in model training
#'   \item Row names used as sample identifiers
#' }
#' @param group_col Name of grouping variable (unused in current implementation). 
#'        Retained for backward compatibility.
#' @param palette_name Name of wesanderson color palette for visualization. 
#'        Default "Royal1".
#' @param save_dir Output directory for saved plots. Created if nonexistent.
#' @param plot_width Plot width in inches. Default 6.
#' @param plot_height Plot height in inches. Default 6.
#' @param alpha Transparency level (unused in current implementation).
#' @param base_size Base font size for plot text. Default 14.
#' @param best_threshold Classification probability threshold (0-1). If NULL, 
#'        uses the threshold stored in the Model_data object.
#'
#' @return A list containing:
#' \itemize{
#'   \item `predictions`: Data frame with:
#'     \itemize{
#'       \item `Sample`: Sample IDs (from row names)
#'       \item `Classification`: Predicted class ("Positive"/"Negative")
#'       \item `Probability`: Prediction probabilities
#'     }
#'   \item `best_threshold`: The probability threshold used for classification
#' }
#'
#' @section Warning:
#' Stops execution if:
#' \itemize{
#'   \item Input is not a Model_data object
#'   \item Required variables are missing in new_data
#'   \item Model or new_data are empty
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with model's default threshold
#' results <- ModelClinicalPrediction(model_obj, new_patient_data)
#' 
#' # Custom threshold and output directory
#' object_best<-ModelClinicalPrediction(
#'   object_best, 
#'   new_patient_data,
#'   best_threshold = 0.7,
#'   save_dir = "results/clinical_predictions"
#' )
#' }
#' 
#' @export
ModelClinicalPrediction <- function(
    object,
    new_data,
    group_col = "group",
    palette_name = "Royal1",
    save_dir = here("ModelData", "clinical_predictions"),
    plot_width = 6,
    plot_height = 6,
    alpha = 1,
    base_size = 14,
    best_threshold = NULL
) {
  if (inherits(object, "Model_data")) {
    cat("Input is of class 'Model_data'. Extracting datasets...\n")
    
    best_model <- object@best.model.result[["model"]][[1]]
    training_data <- object@filtered.set[["training"]]
    
    if (is.null(best_threshold)) {
      best_threshold <- object@best.model.result[["best_threshold"]]
    }
    
    if (is.null(new_data) || nrow(new_data) == 0) {
      stop("new_data dataset is empty or not found.")
    }
    if (is.null(best_model)) {
      stop("Best model not found in the object.")
    }
    
    var_imp <- tryCatch(
      {
        varImp(best_model)$importance
      },
      error = function(e) {
        cat("Variable importance not available. Treating all variables as important.\n")
        NULL
      }
    )
    
    important_vars <- best_model$coefnames
    train_vars <- colnames(training_data)
    change_vars <- setdiff(important_vars, train_vars)
    changed_vars <- gsub("1$", "", change_vars)
    var_mapping <- setNames(changed_vars, change_vars)
    
    important_vars <- ifelse(
      important_vars %in% change_vars,
      var_mapping[important_vars],
      important_vars
    )
    
    missing_vars <- setdiff(important_vars, names(new_data))
    if (length(missing_vars) > 0) {
      stop(paste(
        "The following important variables are missing in the data:",
        paste(missing_vars, collapse = ", ")
      ))
    }
    
    prediction_data <- new_data[, important_vars]
    
    predictions <- predict(
      best_model,
      newdata = prediction_data,
      type = "prob"
    )[, 2]
    
    predicted_labels <- ifelse(predictions > best_threshold, 1, 0)
    
    final_result <- data.frame(
      Sample = rownames(prediction_data),
      Classification = as.factor(ifelse(
        predicted_labels == 1,
        "Positive",
        "Negative"
      )),
      Probability = predictions
    )
    
    Classification <- unique(final_result$Classification)
    colors <- wes_palette(
      n = length(Classification),
      name = palette_name,
      type = "discrete"
    )
    names(colors) <- Classification
    
    p <- ggplot(final_result, aes(x = Classification, y = Probability, fill = Classification)) +
      geom_boxplot(
        outlier.shape = 19,
        outlier.colour = colors[1],
        outlier.size = 1
      ) +
      geom_jitter(
        width = 0.2,
        size = 2,
        aes(color = Classification),
        alpha = 0.6
      ) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(
        title = "Visualization of Predicted Group and Probabilities",
        x = "Group",
        y = "Predicted Probability"
      ) +
      theme_minimal(base_size = base_size) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = base_size + 2
        ),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = base_size - 2
        ),
        axis.text.y = element_text(size = base_size - 2),
        axis.title.x = element_text(size = base_size),
        axis.title.y = element_text(size = base_size),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = base_size - 2)
      )
    
    print(p)
    
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    
    ggsave(
      filename = file.path(save_dir, "prediction_visualization.pdf"),
      plot = p,
      width = plot_width,
      height = plot_height,
      device = "pdf"
    )
    
    cat(
      "Prediction visualization saved to:",
      file.path(save_dir, "prediction_visualization.pdf"), "\n"
    )
    
    return(list(
      predictions = final_result,
      best_threshold = best_threshold
    ))
  } else {
    stop("Input must be an object of class 'Model_data'")
  }
}
