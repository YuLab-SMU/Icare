#' Filter and Subset Data Features
#'
#' Filters datasets (training, testing, validation) based on selected feature subsets and optionally 
#' exports the filtered data to CSV files. Works with objects of class 'Train_Model'.
#'
#' @param object A `Train_Model` object containing the datasets to be filtered
#' @param feature_subset_name Character name of the feature subset to use (default: "best_features_subset"). 
#'        Must correspond to a named list element in object@feature.result
#' @param group_col Character name of the grouping/response variable column (default: "group")
#' @param data_type Character specifying data type to use: "clean" for raw data or "scale" for 
#'        standardized data (default: "clean")
#' @param use_feature_subset Logical indicating whether to filter using feature subset (default: TRUE)
#' @param save_data Logical indicating whether to save filtered data to CSV files (default: TRUE)
#' @param train_filename Character filename for training data CSV output (default: "train_data.csv")
#' @param test_filename Character filename for test data CSV output (default: "test_data.csv")
#' @param val_filename Character filename for validation data CSV output (default: "val_data.csv")
#' @param ext_val_filename Character filename for external validation data CSV output (default: "ext_val_data.csv")
#' @param save_dir Character directory path where CSV files will be saved 
#'        (default: here::here("ModelData", "Data"))
#'
#' @return Returns the modified `Train_Model` object with filtered datasets stored in the 
#'         `@filtered.set` slot. When `save_data=TRUE`, also saves CSV files to disk.
#'
#' @details This function performs the following operations:
#' \enumerate{
#'   \item Validates input object and data structure
#'   \item Extracts datasets based on specified data_type ("clean" or "scale")
#'   \item Optionally filters datasets using specified feature subset
#'   \item Ensures grouping variable is included in filtered data
#'   \optionally saves filtered datasets as CSV files
#'   \item Updates object@filtered.set with filtered datasets
#' }
#'
#' @section File Output:
#' When `save_data=TRUE`, the function saves these files (if datasets exist):
#' \itemize{
#'   \item Training data (train_data.csv)
#'   \item Testing data (test_data.csv)
#'   \item Validation data (val_data.csv)
#'   \item External validation data (ext_val_data.csv)
#' }
#'
#' @examples
#' \dontrun{
#' # Custom parameters example
#' object_model <- FilterDataFeatures(
#'   object = object_model,
#'   feature_subset_name = "selected_features",
#'   data_type = "scale",
#'   save_dir = "output/filtered_data"
#' )
#' }
#'
#' @export
FilterDataFeatures <- function(
    object,
    feature_subset_name = "best_features_subset",
    group_col = "group",
    data_type = "clean",
    use_feature_subset = TRUE,
    save_data = TRUE,
    train_filename = "train_data.csv",
    test_filename = "test_data.csv",
    val_filename = "val_data.csv",
    ext_val_filename = "ext_val_data.csv",
    save_dir = here::here("ModelData", "Data")
) {
  if (!inherits(object, 'Train_Model')) {
    stop("Input must be an object of class 'Train_Model'.")
  }
  
  group_col <- object@group_col %||% group_col
  
  if (!is.list(object@split.data) ||
      !all(c("training", "testing") %in% names(object@split.data))) {
    stop("The 'split.data' slot must contain 'training' and 'testing' datasets.")
  }
  
  cat("Data type:", data_type, "\n")
  cat("Use feature subset:", use_feature_subset, "\n")
  if (use_feature_subset) {
    cat("Feature subset name:", feature_subset_name, "\n")
  }
  cat("Save to CSV:", save_data, "\n")
  if (save_data) {
    cat("Save directory:", save_dir, "\n")
  }
  cat("\n")
  
  if (data_type == "scale") {
    train_data <- slot(object, "split.scale.data")[["training"]]
    test_data <- slot(object, "split.scale.data")[["testing"]]
    validation <- if ("validation" %in% names(slot(object, "split.scale.data"))) 
      slot(object, "split.scale.data")[["validation"]] else NULL
    external_validation <- if ("external_validation" %in% names(slot(object, "split.scale.data"))) 
      slot(object, "split.scale.data")[["external_validation"]] else NULL
  } else if (data_type == "clean") {
    train_data <- slot(object, "split.data")[["training"]]
    test_data <- slot(object, "split.data")[["testing"]]
    validation <- if ("validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["validation"]] else NULL
    external_validation <- if ("external_validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["external_validation"]] else NULL
  } else {
    stop("Invalid 'data_type'. Use 'clean' or 'scale'.")
  }
  
  if (is.null(train_data) || !is.data.frame(train_data) || nrow(train_data) == 0) {
    stop("No valid training data found.")
  }
  if (is.null(test_data) || !is.data.frame(test_data) || nrow(test_data) == 0) {
    stop("No valid test data found.")
  }
  
  if (!use_feature_subset) {
    cat("Skipping feature subset filtering - using full dataset\n")
    filtered_data <- list(
      training = train_data,
      testing = test_data,
      validation = if (!is.null(validation)) validation else NULL,
      external_validation = if (!is.null(external_validation)) external_validation else NULL
    )
  } else {
    best_features <- object@feature.result[[feature_subset_name]]
    
    if (is.null(best_features) || length(best_features) == 0) {
      cat("Warning: No feature subset found - using full dataset\n")
      filtered_data <- list(
        training = train_data,
        testing = test_data,
        validation = if (!is.null(validation)) validation else NULL,
        external_validation = if (!is.null(external_validation)) external_validation else NULL
      )
    } else {
      if (!is.null(group_col) && group_col %in% names(train_data)) {
        best_features <- unique(c(best_features, group_col))
      }
      
      available_train <- intersect(best_features, names(train_data))
      missing_train <- setdiff(best_features, names(train_data))
      if (length(missing_train) > 0) {
        cat("Note: The following features are missing from training data and will be skipped:\n")
        cat(paste(missing_train, collapse = ", "), "\n")
      }
      
      available_test <- intersect(best_features, names(test_data))
      
      filtered_validation <- if (!is.null(validation)) {
        available_valid <- intersect(best_features, names(validation))
        if (length(setdiff(best_features, names(validation))) > 0) {
          cat("Note: Some features missing from validation data\n")
        }
        validation[, available_valid, drop = FALSE]
      } else NULL
      
      filtered_external <- if (!is.null(external_validation)) {
        available_external <- intersect(best_features, names(external_validation))
        if (length(setdiff(best_features, names(external_validation))) > 0) {
          cat("Note: Some features missing from external validation data\n")
        }
        external_validation[, available_external, drop = FALSE]
      } else NULL
      
      filtered_data <- list(
        training = train_data[, available_train, drop = FALSE],
        testing = test_data[, available_test, drop = FALSE],
        validation = filtered_validation,
        external_validation = filtered_external
      )
    }
  }
  
  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("Created directory:", save_dir, "\n")
      }
      
      write.csv(filtered_data$training, 
                file.path(save_dir, train_filename), 
                row.names = FALSE)
      
      write.csv(filtered_data$testing, 
                file.path(save_dir, test_filename), 
                row.names = FALSE)
      
      if (!is.null(filtered_data$validation)) {
        write.csv(filtered_data$validation, 
                  file.path(save_dir, val_filename), 
                  row.names = FALSE)
      }
      
      if (!is.null(filtered_data$external_validation)) {
        write.csv(filtered_data$external_validation, 
                  file.path(save_dir, ext_val_filename), 
                  row.names = FALSE)
      }
      
      cat("\nSaved CSV files to:\n")
      cat("-", file.path(save_dir, train_filename), "\n")
      cat("-", file.path(save_dir, test_filename), "\n")
      if (!is.null(filtered_data$validation)) {
        cat("-", file.path(save_dir, val_filename), "\n")
      }
      if (!is.null(filtered_data$external_validation)) {
        cat("-", file.path(save_dir, ext_val_filename), "\n")
      }
    }, error = function(e) {
      warning("Failed to save CSV files: ", e$message)
    })
  }
  
  object@filtered.set <- filtered_data
  cat("\nUpdated 'Train_Model' object contains:\n")
  cat("- Filtered datasets in filtered.set slot:\n")
  cat("  * Training (", ncol(filtered_data$training)-1, " features)\n", sep="")
  cat("  * Testing (", ncol(filtered_data$testing)-1, " features)\n", sep="")
  if (!is.null(filtered_data$validation)) {
    cat("  * Validation (", ncol(filtered_data$validation)-1, " features)\n", sep="")
  }
  if (!is.null(filtered_data$external_validation)) {
    cat("  * External validation (", ncol(filtered_data$external_validation)-1, " features)\n", sep="")
  }
  
  return(object)
}

