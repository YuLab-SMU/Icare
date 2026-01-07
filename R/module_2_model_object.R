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
#'
#' @examples
#' object_model <- new("Train_Model")
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
  )
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
#'
#' @examples
#' clean_df <- ExtractCleanData(my_stat_object)
#'
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
#' @param ... Additional parameters for future extensions (currently unused).
#'
#' @returns A `Train_Model` S4 object, which contains the data and grouping information for further modeling.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a model object from a data frame
#' df <- data.frame(a = c(1, 2), b = c(3, 4))
#' object_model <- CreateModelObject(data = df, group_col = "a")
#'
#' # Create a model object from a Stat object (assuming a Stat object exists)
#' object_model <- CreateModelObject(object = object_stat)
#' }
CreateModelObject <- function(
    data = NULL,
    object = NULL,
    group_col = NULL,
    ...
) {

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

  }
  else {
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

  if (!is.null(group_col) && !group_col %in% colnames(data.df)) {
    stop(paste("Specified group_col", group_col, "not found in the data columns."))
  }
  data.df<-modify_column_names(data.df)
  data.df <- as.data.frame(data.df)
  Train_Model_instance <- new(
    'Train_Model',
    data.df = data.df,
    group_col = group_col
  )

  cat("Model object created successfully.\n")
  return(Train_Model_instance)
}
