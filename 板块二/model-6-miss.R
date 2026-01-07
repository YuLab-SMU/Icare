#' Apply Imputation to New Data Based on Training Set Imputation
#'
#' This function applies the same imputation strategy used on the training set
#' to new data (typically a test set), ensuring consistent handling of missing values.
#' It handles both numeric and categorical variables, maintaining factor levels
#' from the training set.
#'
#' @import methods 
#' @import stats 
#' @param new_data A data frame containing the new data to be imputed
#' @param imputation_info A list containing imputation information from impute_missing_values(),
#'                        including imputation values and variable types
#' @return A data frame with missing values imputed consistently with the training set
#' @export
#'
#' @examples
#' # Example 1: Apply imputation to test set
#' imputed_test <- apply_imputation(test_data, train_imputation_info)
#' 
#' # Example 2: Apply to new dataset with warnings for missing variables
#' new_data_imputed <- apply_imputation(new_data, imputation_info = model_imputation)
apply_imputation <- function(new_data, imputation_info) {
  set.seed(123) 
  if (!is.data.frame(new_data)) stop("new_data must be a data frame")
  if (!is.list(imputation_info) || is.null(imputation_info$imputation_values)) {
    stop("imputation_info must be a valid list from impute_missing_values")
  }
  
  imputed_data <- as.data.frame(new_data)
  original_vars <- names(imputation_info$imputation_values)
  missing_vars <- setdiff(original_vars, colnames(imputed_data))
  
  if (length(missing_vars) > 0) {
    warning(paste("Missing variables in new_data:", paste(missing_vars, collapse = ", ")))
  }
  
  numeric_vars <- imputation_info$variable_types$numeric_vars
  categorical_vars <- imputation_info$variable_types$categorical_vars
  
  for (col in intersect(original_vars, colnames(imputed_data))) {
    tryCatch({
      col_info <- imputation_info$imputation_values[[col]]
      imp_value <- if (!is.null(col_info$used_value)) col_info$used_value else col_info$value
      
      if (col %in% numeric_vars) {
        if (!is.numeric(imputed_data[[col]])) {
          warning(sprintf("Converting variable %s to numeric", col))
          imputed_data[[col]] <- as.numeric(imputed_data[[col]])
        }
        imputed_data[[col]][is.na(imputed_data[[col]])] <- imp_value
      } else if (col %in% categorical_vars) {
        imputed_data[[col]] <- as.character(imputed_data[[col]])
        imputed_data[[col]][is.na(imputed_data[[col]])] <- as.character(imp_value)
        
        if (col %in% names(imputation_info$imputed_data)) {
          original_levels <- levels(imputation_info$imputed_data[[col]])
          imputed_data[[col]] <- factor(imputed_data[[col]], levels = original_levels)
        } else {
          imputed_data[[col]] <- as.factor(imputed_data[[col]])
        }
      } else {
        imputed_data[[col]][is.na(imputed_data[[col]])] <- imp_value
      }
    }, error = function(e) {
      warning(sprintf("Failed to impute variable %s: %s", col, e$message))
    })
  }
  
  for (col in intersect(categorical_vars, colnames(imputed_data))) {
    if (!is.factor(imputed_data[[col]])) {
      imputed_data[[col]] <- as.factor(imputed_data[[col]])
    }
  }
  
  return(imputed_data)
}

#' Apply Missing Value Imputation to Model Data
#'
#' This function handles missing values in training, testing, validation, and external validation datasets
#' using specified imputation methods. It works with both Train_Model objects and raw data lists, with
#' optional CSV export of imputed datasets.
#'
#' @param object Either a Train_Model object or a list containing data frames (training, testing, etc.)
#' @param group_col The name of the column containing group/class information (default: "group")
#' @param impute_method The imputation method to use ("mice", "mean", "median", etc.) (default: "mice")
#' @param m Number of multiple imputations (only relevant for "mice" method) (default: 5)
#' @param max_unique_values Maximum unique values for a variable to be considered categorical (default: 5)
#' @param return_imputation_info Whether to return imputation parameters for later application (default: TRUE)
#' @param save_data Logical indicating whether to save imputed datasets as CSV files (default: FALSE)
#' @param train_filename Name for imputed training set CSV file (default: "train_data.csv")
#' @param test_filename Name for imputed test set CSV file (default: "test_data.csv")
#' @param val_filename Name for imputed validation set CSV file (default: "val_data.csv")
#' @param ext_val_filename Name for imputed external validation CSV file (default: "ext_val_data.csv")
#' @param save_dir Directory path for saving CSV files (default: here::here("ModelData","Data"))
#'
#' @return Depending on input:
#'   - For Train_Model objects: Returns updated object with imputed data in split.data slot
#'   - For lists: Returns a list containing imputed datasets and imputation information
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # With Train_Model object and CSV export
#' object_model <- ModelApplyMiss(object_model, 
#'                        impute_method = "median",
#'                        save_data = TRUE)
#'
#' }
ModelApplyMiss <- function(object,
                           group_col = "group",
                           impute_method = "mice",
                           m = 5,
                           max_unique_values = 5,
                           return_imputation_info = TRUE,
                           save_data = TRUE,
                           train_filename = "train_data.csv",
                           test_filename = "test_data.csv",
                           val_filename = "val_data.csv",
                           ext_val_filename = "ext_val_data.csv",
                           save_dir = here::here("ModelData","Data")) {
  
  set.seed(123) 
  cat("\n=== Starting Missing Value Imputation Process ===\n")
  
  if (inherits(object, 'Train_Model')) {
    group_col <- object@group_col
    training <- slot(object, "split.data")[["training"]]
    testing <- slot(object, "split.data")[["testing"]]
    validation <- if ("validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["validation"]] else NULL
    external_validation <- if ("external_validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["external_validation"]] else NULL
    return_object_modelect <- TRUE
  } else if (is.list(object)) {
    training <- object$training
    testing <- object$testing
    validation <- if ("validation" %in% names(object)) object$validation else NULL
    external_validation <- if ("external_validation" %in% names(object)) object$external_validation else NULL
    return_object_modelect <- FALSE
  } else {
    stop("Input must be an object of class 'Train_Model' or a list")
  }
  
  if (is.null(training) || nrow(training) == 0) {
    stop("No valid data found in the input")
  }
  
  if (length(group_col) == 0 || is.null(group_col) || !group_col %in% colnames(training)) {
    cat("Group column is not valid, setting to NULL.\n")
    group_col <- NULL
  }
  
  cat("\n[Step 1] Imputing training data using '", impute_method, "' method...\n", sep = "")
  impute_result <- impute_missing_values(
    data = training,
    group_col = group_col,
    impute_method = impute_method,
    m = m,
    return_imputation_info = return_imputation_info
  )
  
  cat("\n[Step 2] Applying imputation to test data...\n")
  testing_miss <- apply_imputation(new_data = testing, impute_result$imputation_info)
  
  if (!is.null(validation)) {
    cat("\n[Step 3] Applying imputation to validation data...\n")
    validation_miss <- apply_imputation(new_data = validation, impute_result$imputation_info)
  } else {
    validation_miss <- NULL
  }
  
  if (!is.null(external_validation)) {
    cat("\n[Step 4] Applying imputation to external validation data...\n")
    external_validation_miss <- apply_imputation(new_data = external_validation, impute_result$imputation_info)
  } else {
    external_validation_miss <- NULL
  }
  
  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("\nCreated directory:", save_dir, "\n")
      }
      
      save_imputed_data <- function(data, filename) {
        if (!is.null(data)) {
          full_path <- file.path(save_dir, filename)
          write.csv(data, file = full_path, row.names = FALSE)
          cat("- Saved", filename, "to", save_dir, "\n")
        }
      }
      
      cat("\n[Step 5] Saving imputed datasets to CSV...\n")
      save_imputed_data(impute_result[["imputed_data"]], train_filename)
      save_imputed_data(testing_miss, test_filename)
      save_imputed_data(validation_miss, val_filename)
      save_imputed_data(external_validation_miss, ext_val_filename)
      
    }, error = function(e) {
      warning("Failed to save CSV files: ", e$message)
    })
  }
  
  if (return_object_modelect) {
    if ("split.data" %in% slotNames(object)) {
      object@split.data[["training"]] <- impute_result[["imputed_data"]]
      object@split.data[["testing"]] <- testing_miss
      if (!is.null(validation_miss)) object@split.data[["validation"]] <- validation_miss
      if (!is.null(external_validation_miss)) object@split.data[["external_validation"]] <- external_validation_miss
      
      object@process.info[["missing_info"]] <- list(
        impute_method = impute_method,
        imputation_info = impute_result$imputation_info
      )
      
      cat("\nThe 'Train_Model' object has been updated with the following slots:\n")
      cat("- 'split.data' slot updated with imputed data.\n")
      cat("- 'process.info' slot updated with imputation details.\n")
      return(object)
    } else {
      stop("The 'Train_Model' object does not have a 'split.data' slot.")
    }
  }
  
  result <- list(
    training = impute_result[["imputed_data"]],
    testing = testing_miss,
    validation = validation_miss,
    external_validation = external_validation_miss,
    impute_method = impute_method,
    imputation_info = impute_result$imputation_info
  )
  
  cat("\n=== Imputation Completed ===\n")
  if (save_data) cat("Imputed datasets saved to:", save_dir, "\n")
  
  return(result)
}