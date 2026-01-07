#' Process New Data for Model Prediction
#'
#' Prepares new data for prediction by applying the same preprocessing steps used during model training.
#' Handles missing values, outliers, and normalization consistently with the training procedure.
#' Includes robust error handling to ensure graceful degradation when components are missing.
#'
#' @param object A `Best_Model` object containing preprocessing information and trained model.
#' Must have slots: `group_col`, `process.info`, and `filtered.set`.
#' @param new_data Data frame containing new observations to process. Must contain
#' at least the same features as the training data (excluding the group column).
#' @param miss_threshold Numeric. Maximum percentage of missing values allowed per sample.
#' Default is 25. Will be overridden by value stored in `object` if available.
#' @param save_plots Logical. Should diagnostic plots be saved? Default is FALSE.
#' @param save_dir Character. Directory path for saving outputs. Default is "ModelData/clinical_predictions".
#' @param plot_width Numeric. Plot width in inches. Default is 5.
#' @param plot_height Numeric. Plot height in inches. Default is 5.
#' @param base_size Numeric. Base font size for plots. Default is 14.
#' @param plot_display_num Numeric. Number of plots to display if not saving. Default is 1.
#' @param palette_name Character. Color palette name for plots. Default is "Royal1".
#' @param data_type Character. Type of data to return: "clean" (processed only) or
#' "scale" (with normalization applied). Default is "clean".
#' @param save_data Logical. Should processed data be saved to CSV? Default is TRUE.
#' @param csv_filename Character. Filename for saved data. Default is "new_data.csv".
#'
#' @return A processed data frame ready for model prediction with:
#' \itemize{
#'   \item Consistent variable selection as training data
#'   \item Handled missing values (imputed or removed)
#'   \item Processed outliers (according to training method)
#'   \item Factor levels matching training data
#'   \item Optional normalization applied
#' }
#' Returns original data with warning if critical preprocessing information is missing.
#'
#' @details The function performs the following processing steps:
#' \enumerate{
#'   \item Removes group column if present in new data
#'   \item Selects and reorders variables to match training set
#'   \item Removes samples with excessive missing values
#'   \item Imputes missing values using training method
#'   \item Processes outliers based on training parameters
#'   \item Converts factors to match training data levels
#'   \item Optionally applies normalization
#' }
#' Each step includes error handling and will skip with warning if required information is missing.
#'
#' @examples
#' \dontrun{
#' processed_data <- process_new_data(object_best, new_data)
#' }
#'
#' @export
process_new_data <- function(
    object,
    new_data,
    miss_threshold = 25,
    save_plots = FALSE,
    save_dir = here::here("ModelData", "clinical_predictions"),
    plot_width = 5,
    plot_height = 5,
    base_size = 14,
    plot_display_num = 1,
    palette_name = "Royal1",
    data_type = "clean",
    save_data = TRUE,
    csv_filename = "new_data.csv") {

  if (!inherits(object, "Best_Model")) {
    stop("Input object must be a Best_Model object.")
  }
  if (is.null(new_data) || nrow(new_data) == 0) {
    stop("new_data must be a non-empty data frame.")
  }

  group_col_in_object <- tryCatch(
    object@group_col,
    error = function(e) NULL
  )

  if (!is.null(group_col_in_object) && group_col_in_object %in% colnames(new_data)) {
    warning(paste("Removing group column", group_col_in_object,
                  "from new_data as it should not be present in prediction data"))
    new_data <- new_data[, !colnames(new_data) %in% group_col_in_object, drop = FALSE]
  }

  original_rownames <- rownames(new_data)

  mapping_info <- tryCatch(
    object@process.info[["mapping_info"]],
    error = function(e) NULL
  )
  missing_removal <- tryCatch(
    object@process.info[["missing_removal"]],
    error = function(e) NULL
  )
  missing_info <- tryCatch(
    object@process.info[["missing_info"]],
    error = function(e) NULL
  )
  outlier_info <- tryCatch(
    object@process.info[["outlier_info"]],
    error = function(e) NULL
  )
  retained_vars <- tryCatch(
    colnames(object@filtered.set$training),
    error = function(e) NULL
  )

  if (!is.null(retained_vars) && !is.null(group_col_in_object)) {
    retained_vars <- retained_vars[retained_vars != group_col_in_object]
  }

  if (is.null(retained_vars)) {
    warning("No retained variables found in model object. Returning original data.")
    return(new_data)
  }

  normalize_info <- tryCatch(
    object@process.info[["normalization_info"]],
    error = function(e) NULL
  )

  clean_new_data <- as.data.frame(lapply(new_data, function(x) {
    x[x %in% c("<NA>", "NA", "", "NULL")] <- NA
    return(x)
  }), stringsAsFactors = FALSE)

  clean_new_data <- clean_new_data[, retained_vars, drop = FALSE]
  rownames(clean_new_data) <- original_rownames[1:nrow(clean_new_data)]
  cat("- After variable selection:", nrow(clean_new_data), "rows x", ncol(clean_new_data), "columns\n")

  if (!is.null(missing_removal)) {
    miss_threshold <- missing_removal$miss_threshold
  } else {
    miss_threshold <- 25
    warning("Using default missing threshold (25%)")
  }

  data_rm <- which(rowMeans(is.na(clean_new_data)) * 100 >= miss_threshold)
  if (length(data_rm) > 0) {
    cat("- Detected high-missing samples:", length(data_rm), "samples (threshold =", miss_threshold, "%)\n")
    clean_new_data <- clean_new_data[-data_rm, , drop = FALSE]
  } else {
    cat("- No high-missing samples found (threshold =", miss_threshold, "%)\n")
  }

  if (!is.null(missing_info)) {
    cat("- Starting missing value imputation...\n")
    testing_miss <- suppressWarnings(
      tryCatch(
        apply_imputation(new_data = clean_new_data, missing_info$imputation_info),
        error = function(e) {
          warning("Imputation failed: ", e$message)
          clean_new_data
        }
      )
    )
    rownames(testing_miss) <- rownames(clean_new_data)
    stopifnot(nrow(testing_miss) == nrow(clean_new_data))
    cat("- Imputation method:", missing_info[["impute_method"]], "\n")
  } else {
    warning("No missing info found. Skipping imputation.")
    testing_miss <- clean_new_data
  }

  if (!is.null(outlier_info)) {
    cat("- Starting outlier processing...\n")
    handle_outliers_train <- tryCatch(
      outlier_info$handle_outliers_train,
      error = function(e) NULL
    )
    normal_ranges_list <- tryCatch(
      outlier_info$outlier_result_train$normal_ranges,
      error = function(e) NULL
    )

    if (!is.null(handle_outliers_train) && !is.null(normal_ranges_list)) {
      outlier_result_test <- tryCatch(
        detect_and_mark_outliers(
          data = testing_miss,
          group_col = NULL,
          palette_name = palette_name,
          save_plots = save_plots,
          save_dir = save_dir,
          plot_display_num = plot_display_num,
          plot_width = plot_width,
          plot_height = plot_height,
          base_size = base_size,
          custom_ranges = normal_ranges_list,
          max_unique_values = max_unique_values
        ),
        error = function(e) {
          warning("Outlier detection failed: ", e$message)
          list(data_marked = testing_miss)
        }
      )
      data_marked <- outlier_result_test$data_marked
      handle_outliers_test <- tryCatch(
        apply_outlier_handling(
          data_marked = data_marked,
          handling_result = handle_outliers_train
        ),
        error = function(e) {
          warning("Outlier handling failed: ", e$message)
          list(cleaned_data = data_marked)
        }
      )
      cleaned_test <- handle_outliers_test$cleaned_data
      cat("- Outlier processing method:", outlier_info$method, "\n")
    } else {
      warning("Skipping outlier processing due to missing information")
      cleaned_test <- testing_miss
    }
  } else {
    warning("No outlier info found. Skipping outlier processing.")
    cleaned_test <- testing_miss
  }

  cleaned_test <- tryCatch(
    convert_factors_to_binary(cleaned_test),
    error = function(e) {
      warning("Factor conversion failed: ", e$message)
      cleaned_test
    }
  )

  factor_vars <- tryCatch(
    sapply(object@filtered.set$training, is.factor),
    error = function(e) NULL
  )

  if (!is.null(factor_vars)) {
    for (var in names(factor_vars)[factor_vars]) {
      if (var %in% colnames(cleaned_test)) {
        if (is.factor(object@filtered.set$training[[var]])) {
          cleaned_test[[var]] <- tryCatch(
            factor(cleaned_test[[var]],
                   levels = levels(object@filtered.set$training[[var]])),
            error = function(e) {
              warning("Factor level matching failed for ", var, ": ", e$message)
              cleaned_test[[var]]
            }
          )
        } else {
          cleaned_test[[var]] <- tryCatch(
            as(cleaned_test[[var]], class(object@filtered.set$training[[var]])),
            error = function(e) {
              warning("Type conversion failed for ", var, ": ", e$message)
              cleaned_test[[var]]
            }
          )
        }
      }
    }
  }

  if (data_type == "scale") {
    if (!is.null(normalize_info)) {
      cat("- Applying data normalization...\n")
      normalization_test <- tryCatch(
        apply_previous_normalization(
          new_data = cleaned_test,
          prev_normalization_info = normalize_info,
          normalize_method = normalize_info[["normalize_method"]],
          group_col = NULL),
        error = function(e) {
          warning("Normalization failed: ", e$message)
          cleaned_test
        }
      )
    } else {
      warning("No normalization info found, returning unscaled data")
      normalization_test <- cleaned_test
    }
  } else {
    normalization_test <- cleaned_test
    cat("- Returning raw (unscaled) data\n")
  }

  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("Created directory: ", save_dir)
    }

    data_path <- file.path(save_dir, csv_filename)
    tryCatch(
      {
        write.csv(normalization_test, data_path, row.names = TRUE)
        cat("Saved processed data to: ", data_path)
      },
      error = function(e) {
        warning("Failed to save data: ", e$message)
      }
    )
  }

  return(normalization_test)
}


#' Clinical Prediction Visualization for Best Model
#'
#' Generates predictions on new data using the best model from a Best_Model object,
#' applies the optimal classification threshold, and creates visualizations of the prediction results.
#'
#' @param object A Best_Model object containing a trained model
#' @param new_data Data frame containing new observations for prediction
#' @param group_col Name of the column containing class labels (default: "group")
#' @param palette_name Name of color palette for plots (default: "Royal1")
#' @param save_dir Directory to save plots (default: "ModelData/clinical_predictions")
#' @param plot_width Plot width in inches (default: 5)
#' @param plot_height Plot height in inches (default: 5)
#' @param alpha Transparency level for points (default: 1)
#' @param base_size Base font size for plots (default: 14)
#' @param best_threshold Optional threshold value to override stored threshold (default: NULL)
#'
#' @return A list containing:
#' \itemize{
#'   \item predictions - Data frame with sample IDs, classifications and probabilities
#'   \item best_threshold - The threshold value used for classification
#' }
#'
#' @examples
#' \dontrun{
#' # After creating a Best_Model object
#' prediction_results <- ModelClinicalPrediction(
#'   object = best_model,
#'   new_data = test_data,
#'   palette_name = "Royal1"
#' )
#'
#' # Access prediction results
#' head(prediction_results$predictions)
#' print(prediction_results$best_threshold)
#' }
#'
#' @importFrom wesanderson wes_palette
#' @importFrom here here
#' @import ggplot2
#' @export
ModelClinicalPrediction <- function(
    object,
    new_data,
    group_col = "group",
    palette_name = "Royal1",
    save_dir = here("ModelData", "clinical_predictions"),
    plot_width = 5,
    plot_height = 5,
    alpha = 1,
    base_size = 14,
    best_threshold = NULL,
    save_data = TRUE,
    csv_filename = "clinical_predictions.csv"
) {

  if (inherits(object, "Best_Model")) {

    best_model <- object@best.model[[1]]
    best_model_type <- object@best.model.type
    group_col <- object@group_col
    data_sets <- object@filtered.set
    train_data <- data_sets$training

    mapping_info <- object@process.info[["mapping_info"]]

    cat("Extracted model components:\n")
    cat(" - Best model type:", best_model_type, "\n")

    if (is.null(best_threshold)) {
      best_threshold <- object@best_threshold[["best_threshold"]]
      if (is.null(best_threshold)) {
        stop("No threshold found in the object and none provided")
      }
      cat("Using default best threshold from object:", best_threshold, "\n")
    } else {
      cat("Using provided best threshold:", best_threshold, "\n")
    }

    if (is.null(new_data) || nrow(new_data) == 0) {
      stop("new_data dataset is empty or not found.")
    }

    if (is.null(best_model)) {
      stop("Best model not found in the object.")
    }

    predictions <- predict(
      best_model,
      newdata = new_data,
      type = "prob"
    )[, 2]

    predicted_labels <- ifelse(predictions > best_threshold, 1, 0)

    if (!is.null(mapping_info) && !is.null(mapping_info$binary)) {
      reverse_mapping <- setNames(names(mapping_info$binary), mapping_info$binary)
      original_labels <- reverse_mapping[as.character(predicted_labels)]
    } else {
      original_labels <- ifelse(predicted_labels == 1, "Positive", "Negative")
      warning("No label mapping information found in the object, using default labels")
    }

    final_result <- data.frame(
      Sample = rownames(new_data),
      Classification = original_labels,
      Probability = predictions,
      stringsAsFactors = FALSE
    )

    Classification <- unique(final_result$Classification)

    colors <- wes_palette(
      n = length(Classification),
      name = palette_name,
      type = "discrete"
    )

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

    output_file <- file.path(save_dir, "prediction_visualization.pdf")
    cat("Saving plot to:", output_file, "\n")
    ggsave(
      filename = output_file,
      plot = p,
      width = plot_width,
      height = plot_height,
      device = "pdf"
    )
    if (save_data) {
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("Created directory: ", save_dir)
      }

      # Save prediction results
      results_path <- file.path(save_dir, csv_filename)
      write.csv(final_result, results_path, row.names = FALSE)
      cat("Saved prediction results to: ", results_path)
      }

    return(list(
      predictions = final_result,
      best_threshold = best_threshold
    ))
  } else {
    stop("Input must be an object of class 'Best_Model'")
  }
}
