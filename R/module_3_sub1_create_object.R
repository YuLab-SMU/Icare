#' Extract Info Data
#'
#' @param object An object containing info.data slot.
#' @export
ExtractInfoData <- function(object) {
  
  data <- tryCatch(
    slot(object, "info.data"),
    error = function(e) {
      warning("Error extracting 'info.data' from the object: ", e$message)
      return(NULL)
    }
  )
  return(data)
}

#' Subtyping Class
#'
#' @import methods
#' @slot clean.data Data frame.
#' @slot info.data Data frame.
#' @slot scale.data Data frame.
#' @slot Optimal.cluster Any.
#' @slot cluster.results Any.
#' @slot visualization.results Any.
#' @slot clustered.data Data frame.
#' @slot evaluation_results List.
#' @export
Subtyping <- setClass(
  Class = 'Subtyping',
  slots = c(
    clean.data = 'data.frame',
    info.data = 'data.frame',
    scale.data = 'data.frame',
    Optimal.cluster = 'ANY',
    cluster.results = 'ANY',
    visualization.results = 'ANY',
    clustered.data = 'data.frame',
    evaluation_results = 'list'
  ),
  prototype = list(
    clean.data = data.frame(),
    info.data = data.frame(),
    scale.data = data.frame(),
    Optimal.cluster = NULL,
    cluster.results = NULL,
    visualization.results = NULL,
    clustered.data = data.frame(),
    evaluation_results = list()
  )
)

#' Create Subtyping Object
#'
#' @param clean.data Clean data.
#' @param info.data Info data.
#' @param scale.data Scale data.
#' @param Optimal.cluster Optimal cluster.
#' @param cluster.results Cluster results.
#' @param visualization.results Visualization results.
#' @param clustered.data Clustered data.
#' @param evaluation_results Evaluation results.
#' @param object Input object (Stat, Subtyping, PrognosiX, Model_data).
#' @export
CreateSubtypingObject <- function(
    clean.data = NULL,
    info.data = data.frame(),
    scale.data = data.frame(),
    Optimal.cluster = NULL,
    cluster.results = NULL,
    visualization.results = NULL,
    clustered.data = data.frame(),
    evaluation_results = list(),
    object = NULL
) {
  if (is.null(clean.data) && is.null(object)) {
    stop("At least one of 'clean.data' or 'object' must be provided.")
  }
  
  if (!is.null(object)) {
    if (!inherits(object, "Stat") && !inherits(object, "Subtyping") && 
        !inherits(object, "PrognosiX") && !inherits(object, "Model_data") && !inherits(object, "Train_Model")) {
      stop("The 'object' parameter must be an instance of class 'Stat', 'Subtyping', 'Model_data'/'Train_Model' or 'PrognosiX'.")
    }
    
    if (inherits(object, "Stat")) {
      clean.data <- ExtractCleanData(object)
      info.data <- ExtractInfoData(object)
      scale.data <- object@scale.data
      if (is.null(clean.data)) {
        stop("Failed to extract clean data from the provided 'Stat' object.")
      }
    }
    
    else if (inherits(object, "Subtyping")) {
      clean.data <- object@clean.data
      info.data <- object@info.data
      scale.data <- object@scale.data
      Optimal.cluster <- object@Optimal.cluster
      cluster.results <- object@cluster.results
      visualization.results <- object@visualization.results
      clustered.data <- object@clustered.data
      evaluation_results<- object@evaluation_results
    }
    
    else if (inherits(object, "PrognosiX")) {
      clean.data <- object@clean.data
      info.data <- object@info.data
      scale.data <- object@scale.data
    }
    else if (inherits(object, "Model_data") || inherits(object, "Train_Model")) {
      clean.data <- object@clean.df
      # scale.data might be different in Train_Model
      if(.hasSlot(object, "scale.data")) scale.data <- object@scale.data
    }
  }
  
  matched_data <- info.data[rownames(clean.data), ]  
  info.data <- matched_data
  
  if (is.null(clean.data)) {
    stop("The 'clean.data' must be provided.")
  }
  
  ensure_numeric_data <- function(data, data_name) {
    if (is.null(data) || nrow(data) == 0) {
      return(data.frame())
    }
    
    row_names <- rownames(data)
    
    data <- as.data.frame(lapply(data, function(x) {
      if (is.factor(x)) {
        levels <- levels(x)
        if (length(levels) == 2) {
          return(as.numeric(x) - 1)  
        } else {
          warning(paste("Factor variable in", data_name, "has more than 2 levels. It will be converted to numeric as-is."))
          return(as.numeric(x))
        }
      } else if (is.character(x)) {
        x <- as.factor(x)
        levels <- levels(x)
        if (length(levels) == 2) {
          return(as.numeric(x) - 1) 
        } else {
          warning(paste("Character variable in", data_name, "has more than 2 levels. It will be converted to numeric as-is."))
          return(as.numeric(x))
        }
      } else if (is.numeric(x)) {
        return(x)
      } else {
        warning(paste("Column in", data_name, "cannot be converted to numeric. It will be kept as-is."))
        return(x)
      }
    }))
    
    rownames(data) <- row_names
    
    na_columns <- colnames(data)[apply(data, 2, function(x) all(is.na(x)))]
    if (length(na_columns) > 0) {
      warning(paste("The following columns in", data_name, "were entirely NA after conversion and will be removed:", paste(na_columns, collapse = ", ")))
      data <- data[, !colnames(data) %in% na_columns]
    }
    
    if (!all(sapply(data, is.numeric))) {
      warning(paste("Some columns in", data_name, "are not numeric. Please check the data."))
    }
    
    return(data)
  }
  
  prepare_data <- function(data, data_name) {
    if (is.null(data) || nrow(data) == 0) {
      return(data.frame())
    }
    
    if (any(is.na(data))) {
      warning(paste("NA values found in", data_name, ". Rows with NA values will be removed."))
      data <- na.omit(data)
      cat(paste("Removed rows with NA values from", data_name, ". New dimensions:", nrow(data), "rows and", ncol(data), "columns.\n"))
    }
    
    if (anyDuplicated(rownames(data))) {
      warning(paste("Duplicate row names found in", data_name, "; they have been made unique."))
      rownames(data) <- make.unique(rownames(data))
    }
    
    if (anyDuplicated(colnames(data))) {
      warning(paste("Duplicate column names found in", data_name, "; they have been made unique."))
      colnames(data) <- make.unique(colnames(data))
    }
    
    if (is.null(colnames(data))){
      stop(paste(data_name, "is missing column names."))
    }
    
    data <- modify_column_names(data)
    
    cat(paste("Data prepared for", data_name, "with", nrow(data), "rows and", ncol(data), "columns.\n"))
    
    return(data)
  }
  
  clean.data <- ensure_numeric_data(clean.data, "clean.data")
  scale.data <- ensure_numeric_data(scale.data, "scale.data")
  clustered.data <- ensure_numeric_data(clustered.data, "clustered.data")
  
  if (any(is.na(clean.data))) {
    warning(paste("NA values found in clean.data. Rows with NA values will be removed."))
    clean.data <- na.omit(clean.data)
    cat(paste("Removed rows with NA values from clean.data. New dimensions:", nrow(clean.data), "rows and", ncol(clean.data), "columns.\n"))
  }
  
  clean.data <- prepare_data(clean.data, "clean.data")
  
  if (nrow(info.data) == 0) {
    info.data <- data.frame(row.names = row.names(clean.data))
    cat("info.data was created from clean.data with", nrow(info.data), "rows.\n")
  } else if (nrow(info.data) > 0 && nrow(clean.data) > 0) {
    if (!identical(rownames(clean.data), rownames(info.data))) {
      warning("Row names of 'clean.data' and 'info.data' are not identical. They will be unified.")
      rownames(info.data) <- rownames(clean.data)
      cat("Row names of info.data have been unified with clean.data.\n")
    }
  }
  
  Subtyping <- new(
    Class = 'Subtyping',
    clean.data = clean.data,
    info.data = info.data,
    scale.data = scale.data,
    Optimal.cluster = Optimal.cluster,
    cluster.results = cluster.results,
    visualization.results = visualization.results,
    clustered.data = clustered.data,
    evaluation_results = evaluation_results
  )
  
  cat("Subtyping object created successfully.\n")
  
  return(Subtyping)
}
