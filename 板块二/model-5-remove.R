#' Remove Variables and Samples with High Missing Rates
#'
#' This function removes variables and samples with missing rates above a specified threshold
#' from training, testing, validation, and external validation datasets. It works with both
#' Train_Model objects and raw data lists, with optional CSV export of cleaned datasets.
#'
#' @param object Either a Train_Model object or a list containing data frames (training, testing, etc.)
#' @param miss_threshold The maximum allowed percentage of missing values (0-100). 
#'        Variables/samples with missing rates >= this threshold will be removed (default: 25)
#' @param save_data Logical indicating whether to save cleaned datasets as CSV files (default: FALSE)
#' @param train_filename Name for cleaned training set CSV file (default: "train_data.csv")
#' @param test_filename Name for cleaned test set CSV file (default: "test_data.csv")
#' @param val_filename Name for cleaned validation set CSV file (default: "val_data.csv")
#' @param ext_val_filename Name for cleaned external validation CSV file (default: "ext_val_data.csv")
#' @param save_dir Directory path for saving CSV files (default: here::here("ModelData","Data"))
#'
#' @return Depending on input:
#'   - For Train_Model objects: Returns updated object with cleaned data in split.data slot
#'   - For lists: Returns a list containing cleaned datasets and removal information
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # With Train_Model object and CSV export
#' object_model <- ModelRemoveMiss(object_model, 
#'                         miss_threshold = 30,
#'                         save_data = TRUE)
#'
#' }
ModelRemoveMiss <- function(object, 
                            miss_threshold = 25,
                            save_data = TRUE,
                            train_filename = "train_data.csv",
                            test_filename = "test_data.csv",
                            val_filename = "val_data.csv",
                            ext_val_filename = "ext_val_data.csv",
                            save_dir = here::here("ModelData","Data")) {
  
  cat("\n=== Removing Variables and Samples with High Missing Rates ===\n")
  
  if (inherits(object, 'Train_Model')) {
    training <- object@split.data[["training"]]
    testing <- object@split.data[["testing"]]
    validation <- if ("validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["validation"]] else NULL
    external_validation <- if ("external_validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["external_validation"]] else NULL
    return_object_modelect <- TRUE
  } else if (is.list(object)) {
    training <- object$training
    testing <- object$testing
    validation <- object$validation
    external_validation <- object$external_validation
    return_object_modelect <- FALSE
  } else {
    stop("Input should be a Train_Model object or a list containing datasets.")
  }
  
  if (is.null(training) || nrow(training) == 0) {
    stop("Valid training data not found.")
  }
  
  standardize_missing <- function(df) {
    as.data.frame(lapply(df, function(x) {
      x[x %in% c("<NA>", "NA", "", "NULL")] <- NA
      return(x)
    }), stringsAsFactors = FALSE)
  }
  
  cat("\n[Step 1] Processing training data (remove vars + samples)...\n")
  training <- standardize_missing(training)
  
  var_missing_percentage <- colMeans(is.na(training)) * 100
  vars_to_remove <- names(var_missing_percentage[var_missing_percentage >= miss_threshold])
  clean_training <- training[, !colnames(training) %in% vars_to_remove, drop = FALSE]
  
  row_missing_percentage <- rowMeans(is.na(clean_training)) * 100
  samples_to_remove <- which(row_missing_percentage >= miss_threshold)
  if (length(samples_to_remove) > 0) {
    clean_training <- clean_training[-samples_to_remove, , drop = FALSE]
  }
  
  cat("- Removed", length(vars_to_remove), "variables from training.\n")
  cat("- Removed", length(samples_to_remove), "samples from training.\n")
  
  retained_vars <- colnames(clean_training)
  
  clean_testing <- NULL
  if (!is.null(testing)) {
    cat("[Step 2] Processing testing data (remove same vars + high-miss samples)...\n")
    testing <- standardize_missing(testing)
    test_sub <- testing[, retained_vars, drop = FALSE]
    test_rm <- which(rowMeans(is.na(test_sub)) * 100 >= miss_threshold)
    if (length(test_rm) > 0) {
      test_sub <- test_sub[-test_rm, , drop = FALSE]
    }
    clean_testing <- test_sub
    cat("- Removed", length(test_rm), "samples from testing.\n")
  }
  
  clean_validation <- NULL
  if (!is.null(validation)) {
    cat("[Step 3] Processing validation data...\n")
    validation <- standardize_missing(validation)
    val_sub <- validation[, retained_vars, drop = FALSE]
    val_rm <- which(rowMeans(is.na(val_sub)) * 100 >= miss_threshold)
    if (length(val_rm) > 0) {
      val_sub <- val_sub[-val_rm, , drop = FALSE]
    }
    clean_validation <- val_sub
    cat("- Removed", length(val_rm), "samples from validation.\n")
  }
  
  # Process external validation data
  clean_external_validation <- NULL
  if (!is.null(external_validation)) {
    cat("[Step 4] Processing external validation data...\n")
    external_validation <- standardize_missing(external_validation)
    ext_sub <- external_validation[, retained_vars, drop = FALSE]
    ext_rm <- which(rowMeans(is.na(ext_sub)) * 100 >= miss_threshold)
    if (length(ext_rm) > 0) {
      ext_sub <- ext_sub[-ext_rm, , drop = FALSE]
    }
    clean_external_validation <- ext_sub
    cat("- Removed", length(ext_rm), "samples from external validation.\n")
  }
  
  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("\nCreated directory:", save_dir, "\n")
      }
      
      save_clean_data <- function(data, filename) {
        if (!is.null(data)) {
          full_path <- file.path(save_dir, filename)
          write.csv(data, file = full_path, row.names = FALSE)
          cat("- Saved", filename, "to", save_dir, "\n")
        }
      }
      
      cat("\n[Step 5] Saving cleaned datasets to CSV...\n")
      save_clean_data(clean_training, train_filename)
      save_clean_data(clean_testing, test_filename)
      save_clean_data(clean_validation, val_filename)
      save_clean_data(clean_external_validation, ext_val_filename)
      
    }, error = function(e) {
      warning("Failed to save CSV files: ", e$message)
    })
  }
  
  if (return_object_modelect) {
    object@split.data[["training"]] <- clean_training
    object@split.data[["testing"]] <- clean_testing
    if (!is.null(clean_validation)) object@split.data[["validation"]] <- clean_validation
    if (!is.null(clean_external_validation)) object@split.data[["external_validation"]] <- clean_external_validation
    
    object@process.info[["missing_removal"]] <- list(
      miss_threshold = miss_threshold,
      removed_variables = vars_to_remove,
      removed_samples = list(
        training = length(samples_to_remove),
        testing = if (!is.null(testing)) length(test_rm) else NULL,
        validation = if (!is.null(validation)) length(val_rm) else NULL,
        external_validation = if (!is.null(external_validation)) length(ext_rm) else NULL
      )
    )
    
    cat("\nThe 'Train_Model' object has been updated with the following slots:\n")
    cat("- 'split.data' slot updated with the cleaned data.\n")
    cat("- 'process.info' slot updated with missing removal details.\n")
    cat("\n=== Removal Completed. Updated Train_Model returned ===\n")
    return(object)
  }
  
  result <- list(
    training = clean_training,
    testing = clean_testing,
    validation = clean_validation,
    external_validation = clean_external_validation,
    removed_variables = vars_to_remove,
    removed_samples = list(
      training = samples_to_remove,
      testing = if (!is.null(testing)) test_rm else NULL,
      validation = if (!is.null(validation)) val_rm else NULL,
      external_validation = if (!is.null(external_validation)) ext_rm else NULL
    ),
    miss_threshold = miss_threshold
  )
  
  cat("\n=== Removal Completed. Returning cleaned list ===\n")
  return(result)
}