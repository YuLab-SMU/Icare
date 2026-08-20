#' Train_Model S4 Class
#'
#' An S4 class to store cleaned data, grouping information, and various results from the
#' modeling pipeline, including feature selection, model training.
#'
#' @import methods
#' @slot data.df A `data.frame` containing the original input data.
#' @slot clean.df A `data.frame` containing cleaned and preprocessed data.
#' @slot group_col A column name or identifier indicating group labels.
#' @slot process.info A `list` storing information about class balancing.
#' @slot split.data A `list` of training/testing data splits.
#' @slot split.scale.data A `list` of scaled training/testing data.
#' @slot feature.selection A `list` of selected features for modeling.
#' @slot feature.result A `list` containing feature evaluation results.
#' @slot filtered.set A `list` of data or variables filtered by some criteria.
#' @slot train.models A `list` of trained models.
#' @slot all.results A `list` storing evaluation results from all models.
#' @slot best.model.result A `list` containing the best model and its performance.
#'
#' @export
Train_Model <- setClass(
  Class = 'Train_Model',
  slots = c(
    data.df = 'data.frame',
    clean.df = 'data.frame',
    group_col = "ANY",
    process.info = 'list',
    split.data = 'list',
    split.scale.data = 'list',
    feature.selection = 'list',
    feature.result = 'list',
    filtered.set = 'list',
    train.models = 'list',
    all.results = 'list',
    best.model.result = 'list'
  ),
  prototype = list(
    data.df = data.frame(),
    group_col = NULL,
    clean.df = data.frame(),
    process.info = list(),
    split.data = list(),
    split.scale.data = list(),
    feature.result = list(),
    filtered.set = list(),
    all.results = list(),
    feature.selection = list(),
    train.models = list(),
    best.model.result = list()
  ),
  validity = function(object) {
    # Check data frame slots
    if (!is.data.frame(object@data.df)) {
      return("Slot 'data.df' must be a data.frame")
    }
    if (!is.data.frame(object@clean.df)) {
      return("Slot 'clean.df' must be a data.frame")
    }
    
    # Check list slots
    if (!is.list(object@process.info)) {
      return("Slot 'process.info' must be a list")
    }
    if (!is.list(object@split.data)) {
      return("Slot 'split.data' must be a list")
    }
    if (!is.list(object@split.scale.data)) {
      return("Slot 'split.scale.data' must be a list")
    }
    if (!is.list(object@feature.selection)) {
      return("Slot 'feature.selection' must be a list")
    }
    if (!is.list(object@feature.result)) {
      return("Slot 'feature.result' must be a list")
    }
    if (!is.list(object@filtered.set)) {
      return("Slot 'filtered.set' must be a list")
    }
    if (!is.list(object@train.models)) {
      return("Slot 'train.models' must be a list")
    }
    if (!is.list(object@all.results)) {
      return("Slot 'all.results' must be a list")
    }
    if (!is.list(object@best.model.result)) {
      return("Slot 'best.model.result' must be a list")
    }
    
    # Check row name consistency if data exists
    if (nrow(object@data.df) > 0 && nrow(object@clean.df) > 0) {
      if (!setequal(rownames(object@data.df), rownames(object@clean.df))) {
        return("Row names of 'data.df' and 'clean.df' should match")
      }
    }
    
    return(TRUE)
  }
)

#' Extract Clean Data from Stat Object
#'
#' This function extracts the `clean.data` slot from an object of class `Stat`.
#'
#' @param object An object of class `Stat`.
#'
#' @return A `data.frame` containing the cleaned data from the `Stat` object,
#' or `NULL` if extraction fails.
#'
#' @export
ExtractCleanData <- function(object) {
  if (!inherits(object, "Stat")) {
    stop("The input object must be of class 'Stat'.")
  }

  data <- tryCatch(
    slot(object, "clean.data"),
    error = function(e) {
      warning("Error extracting 'clean.data' from the object: ", e$message)
      return(NULL)
    }
  )
  return(data)
}

#' Create a Train_Model Object
#'
#' This function creates a `Train_Model` S4 object based on either a raw data frame or a preprocessed `Stat` object. It is designed to provide a structured object for subsequent modeling analyses.
#'
#' @param data A data.frame containing raw input data. It must not contain missing values. If this parameter is provided, the `object` parameter should be NULL.
#' @param object A `Stat` class object containing preprocessed data. If this parameter is provided, the function will extract `clean.data` and `group_col` from it.
#' @param group_col The column name for the grouping variable. If the input is a `Stat` object, this value will be extracted from the object.
#' @param na.action Character string specifying how to handle NA values. 
#'   Options are "error" (stop if NA found, default for modeling), "omit" (remove rows with NA),
#'   or "allow" (keep NA values - may cause errors in modeling).
#' @param ... Additional parameters for future extensions (currently unused).
#'
#' @returns A `Train_Model` S4 object, which contains the data and grouping information for further modeling.
#'
#' @export
#' @examples
#' # Binary classification with mtcars
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' print(model)
CreateModelObject <- function(
    data = NULL,
    object = NULL,
    group_col = NULL,
    na.action = c("error", "omit", "allow"),
    ...
) {
  
  na.action <- match.arg(na.action)
  
  if (!is.null(data) && !is.null(object)) {
    stop("Only one of 'data' and 'object' should be provided.")
  }
  
  if (is.null(data) && is.null(object)) {
    stop("At least one of 'data' and 'object' must be provided.")
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
    
    # Handle NA values based on na.action
    if (any(is.na(data))) {
      switch(na.action,
             "allow" = {
               message("Note: NA values found in ", data_name, " but kept (na.action='allow')")
             },
             "omit" = {
               warning("NA values found in ", data_name, ". Removing ", sum(is.na(data)), " NA values.")
               data <- na.omit(data)
             },
             "error" = {
               stop("NA values found in ", data_name, ". Set na.action='omit' to remove or na.action='allow' to keep.")
             }
      )
    }
    
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
    data.df <- prepare_data(data, "input data")
  } else {
    if (!inherits(object, "Stat")) {
      stop("The 'object' parameter must be an instance of class 'Stat'.")
    }
    data.df <- ExtractCleanData(object = object)
    if (is.null(data.df)) {
      stop("Failed to extract clean data from the provided 'Stat' object.")
    }
    data.df <- prepare_data(data.df, "clean data from Stat object")
    if (is.null(group_col)) {
      group_col <- object@group_col
    }
  }
  
  # Check group_col BEFORE modifying column names
  original_group_col <- group_col
  if (!is.null(group_col) && !group_col %in% colnames(data.df)) {
    stop(paste("Specified group_col", group_col, "not found in the data columns."))
  }
  
  data.df <- modify_column_names(data.df)
  
  # Update group_col if it was modified by modify_column_names
  if (!is.null(original_group_col)) {
    modified_group_col <- make.names(tolower(gsub(" ", "_", gsub("[^[:alnum:]_]", "", original_group_col))))
    if (modified_group_col != original_group_col && modified_group_col %in% colnames(data.df)) {
      message("Note: group_col '", original_group_col, "' was modified to '", modified_group_col, "'")
      group_col <- modified_group_col
    }
  }
  
  # Create the S4 object
  Train_Model_instance <- new(
    'Train_Model',
    data.df = data.df,
    group_col = group_col
  )
  
  Train_Model_instance@clean.df <- data.df
  
  cat("Model object created successfully.\n")
  return(Train_Model_instance)
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
#' @export
#' @examples
#' \dontrun{
#' data(model_obj_test)
#' # Assuming 'train_obj_test' is a valid Train_Model object
#' extracted <- ExtractModel(model_obj_test)
#' print(extracted$best_model)
#' }
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
    result$best_model <- object@best.model.result[["model"]]
    
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
