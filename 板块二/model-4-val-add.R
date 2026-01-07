#' Extract and Prepare Validation Data for Model Evaluation
#'
#' Processes validation data from either raw data frame or Stat object, preparing it for
#' model evaluation by:
#' 1. Cleaning special characters (>, <) from numeric values
#' 2. Ensuring unique row/column names
#' 3. Automatic variable type conversion
#' 4. One-hot encoding for categorical variables
#' 5. Binary conversion for factors
#' Finally updates the Train_Model object with validation data.
#'
#' @param data Optional raw data frame containing validation data. 
#'   Must be provided if object_stats is NULL.
#' @param object_stats Optional Stat class object containing pre-processed data.
#'   Must be provided if data is NULL.
#' @param object_model Required Train_Model object to be updated with validation data.
#' @param group_col Optional name of the grouping/outcome variable column.
#' @param max_unique_values Maximum number of unique values for a numeric variable to be
#'   considered for conversion to factor (default: 5).
#' @param ... Additional arguments passed to diagnostic functions.
#'
#' @return The updated Train_Model object with validation data stored in the 
#'   split.data$validation slot.
#'
#' @export
#'
#' @examples
#' # Using raw data frame
#' object_model <- new("Train_Model")
#' valid_df <- data.frame(
#'   group = factor(c("A","B","A")),
#'   var1 = c(1, 2, 1),
#'   var2 = c("X", "Y", "X")
#' )
#' updated_model <- Extract_validata(
#'   data = valid_df,
#'   object_model = object_model,
#'   group_col = "group"
#' )
#'
#' # Using Stat object
#' object_model <- new("Stat", clean.df = valid_df)
#' updated_model <- Extract_validata(
#'   object_stats = object_model,
#'   object_model = object_model,
#'   group_col = "group"
#' )
#'
Extract_validata <- function(
    data = NULL,
    object_stats = NULL,
    object_model = NULL,
    group_col = NULL,
    max_unique_values = 5,
    ...
) {
  if (!is.null(data) && !is.null(object_stats)) {
    stop("Only one of 'data' and 'object_stats' should be provided.")
  }
  if (is.null(data) && is.null(object_stats)) {
    stop("At least one of 'data' and 'object_stats' must be provided.")
  }
  if (is.null(object_model)) {
    stop("'object_model' must be provided.")
  }
  
  clean_symbol_values <- function(data) {
    for(col in colnames(data)) {
      if(is.character(data[[col]])) {
        if(any(grepl("^[<>]", data[[col]]))) {
          data[[col]] <- gsub("[<>]", "", data[[col]])
          data[[col]] <- as.numeric(data[[col]])
          warning(paste("Removed >/< symbols from column", col, "and converted to numeric"))
        }
      }
    }
    return(data)
  }
  
  prepare_data <- function(data, data_name) {
    if (nrow(data) == 0) {
      stop(paste(data_name, "is empty."))
    }
    
    data <- as.data.frame(data)
    data <- clean_symbol_values(data)
    
    if (anyDuplicated(rownames(data))) {
      warning(paste("Duplicate row names found in", data_name, "; they have been made unique."))
      rownames(data) <- make.unique(rownames(data))
    }
    
    if (anyDuplicated(colnames(data))) {
      warning(paste("Duplicate column names found in", data_name, "; they have been made unique."))
      colnames(data) <- make.unique(colnames(data))
    }
    
    if (is.null(colnames(data))) {
      stop(paste(data_name, "is missing column names."))
    }
    
    cat(paste("Data prepared for", data_name))
    return(data)
  }
  
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("The 'data' parameter must be a data frame.")
    }
    data.df <- prepare_data(data, "validation data")
  } else {
    if (!inherits(object_stats, "Stat")) {
      stop("The 'object_stats' parameter must be an instance of class 'Stat'.")
    }
    data.df <- ExtractCleanData(object = object_stats)
    if (is.null(data.df)) {
      stop("Failed to extract clean data from the provided 'Stat' object.")
    }
    data.df <- prepare_data(data.df, "clean data from Stat object")
  }
  
  variable_types <- diagnose_variable_type(data.df, 
                                           group_col = group_col, 
                                           max_unique_values = max_unique_values)
  
  converted_data <- convert_variables(data.df, variable_types)
  
  onehot_data <- one_hot_encode(converted_data, 
                                group_col = group_col,
                                max_unique_values = max_unique_values)
  
  final_data <- convert_factors_to_binary(onehot_data)
  
  object_model@split.data[["validation"]] <- final_data
  
  cat("The independent validation set has been added.\n")
  cat("Updating 'Train_Model' object...\n")
  cat("The 'Train_Model' object has been updated with:\n")
  cat("- 'split.data' slot updated with validation data\n")
  
  return(object_model)
}


#' Extract and Prepare External Validation Data for Model Evaluation
#'
#' This function processes external validation data (either provided directly or extracted from a Stat object)
#' and prepares it for use with a trained model object. It performs data cleaning, type conversion,
#' one-hot encoding, and ensures compatibility with the model's requirements.
#'
#' @param data A data frame containing the external validation data. Either this or `object_stats` must be provided.
#' @param object_stats A Stat object from which clean data can be extracted. Either this or `data` must be provided.
#' @param object_model A Train_Model object that will be updated with the external validation data.
#' @param group_col The name of the column containing group/class information (optional).
#' @param ... Additional arguments passed to diagnostic and conversion functions.
#'
#' @return The updated `object_model` with the external validation data added to its `split.data` slot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with direct data input
#' model <- Extract_external_validation(
#'   data = validation_data,
#'   object_model = trained_model,
#'   group_col = "diagnosis"
#' )
#'
#' # Example with Stat object input
#' model <- Extract_external_validation(
#'   object_stats = stat_object,
#'   object_model = trained_model
#' )
#' }
#'
Extract_external_validation <- function(
    data = NULL,
    object_stats = NULL,
    object_model = NULL,
    group_col = NULL,
    ...
) {
  if (!is.null(data) && !is.null(object_stats)) {
    stop("Only one of 'data' and 'object_stats' should be provided.")
  }
  if (is.null(data) && is.null(object_stats)) {
    stop("At least one of 'data' and 'object_stats' must be provided.")
  }
  if (is.null(object_model)) {
    stop("'object_model' must be provided.")
  }
  
  clean_symbol_values <- function(data) {
    for(col in colnames(data)) {
      if(is.character(data[[col]])) {
        if(any(grepl("^[<>]", data[[col]]))) {
          data[[col]] <- gsub("[<>]", "", data[[col]])
          data[[col]] <- as.numeric(data[[col]])
          warning(paste("Removed >/< symbols from column", col, "and converted to numeric"))
        }
      }
    }
    return(data)
  }
  
  prepare_data <- function(data, data_name) {
    if (nrow(data) == 0) {
      stop(paste(data_name, "is empty."))
    }
    
    data <- as.data.frame(data)
    data <- clean_symbol_values(data)
    
    if (anyDuplicated(rownames(data))) {
      warning(paste("Duplicate row names found in", data_name, "; they have been made unique."))
      rownames(data) <- make.unique(rownames(data))
    }
    
    if (anyDuplicated(colnames(data))) {
      warning(paste("Duplicate column names found in", data_name, "; they have been made unique."))
      colnames(data) <- make.unique(colnames(data))
    }
    
    if (is.null(colnames(data))) {
      stop(paste(data_name, "is missing column names."))
    }
    
    cat(paste("Data prepared for", data_name, "\n"))
    return(data)
  }
  
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("The 'data' parameter must be a data frame.")
    }
    data.df <- prepare_data(data, "external validation data")
  } else {
    if (!inherits(object_stats, "Stat")) {
      stop("The 'object_stats' parameter must be an instance of class 'Stat'.")
    }
    data.df <- ExtractCleanData(object = object_stats)
    if (is.null(data.df)) {
      stop("Failed to extract clean data from the provided 'Stat' object.")
    }
    data.df <- prepare_data(data.df, "clean data from Stat object")
  }
  
  original_colnames <- colnames(data.df)
  
  variable_types <- diagnose_variable_type(data.df, 
                                           group_col = group_col, 
                                           max_unique_values = 5)
  
  converted_data <- convert_variables(data.df, variable_types)
  
  onehot_data <- one_hot_encode(converted_data, 
                                group_col = group_col,
                                max_unique_values = 5)
  
  final_data <- convert_factors_to_binary(onehot_data)
  
  common_cols <- intersect(colnames(final_data), original_colnames)
  colnames(final_data)[colnames(final_data) %in% common_cols] <- 
    original_colnames[original_colnames %in% common_cols]
  
  object_model@split.data[["external_validation"]] <- final_data
  
  cat("The independent external validation set has been added.\n")
  cat("Updating 'Train_Model' object...\n")
  cat("The 'Train_Model' object has been updated with the following slots:\n")
  cat("- 'split.data' slot updated with external validation data.\n")
  
  return(object_model)
}
