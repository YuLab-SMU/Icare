#' Modify column names of a data frame
#'
#' This function modifies the column names of a data frame by performing the following:
#' 1. Converts all column names to lowercase.
#' 2. Replaces spaces with underscores.
#' 3. Removes any non-alphanumeric characters except for underscores.
#' 4. Ensures all column names are unique by appending numbers to duplicated names.
#'
#' @param data A data frame whose column names need to be modified.
#' @param tolower Logical. If \code{TRUE}, column names are converted to
#'   lowercase before other modifications. Default is \code{FALSE}.
#' @returns A data frame with modified column names.
#' @export
#'
#' @examples
#' data <- data.frame("First Name" = c("John", "Jane"), "Age" = c(25, 30))
#' modified_data <- modify_column_names(data)
#' print(modified_data)

modify_column_names <- function(data,tolower=F) {
  if (is.null(names(data))) stop("Column names cannot be NULL!")
  if(tolower){ new_names <- tolower(names(data)) }   
  new_names=names(data)
  new_names <- gsub(" ", "_", new_names)
  new_names <- gsub("[^[:alnum:]_]", "", new_names)
  if (anyDuplicated(new_names)) {
    warning("Duplicate column names found; made unique.")
    new_names <- make.unique(new_names)
  }
  names(data) <- new_names
  return(data)
}


#' Class to store statistical analysis results
#'
#' This class is designed to store and manage various components of statistical analysis.
#' It includes raw and cleaned data, processed results, descriptive statistics, correlation results,
#' variable types, and metadata related to the analysis.
#'
#' @import methods
#' @slot raw.data A data.frame containing the raw input data.
#' @slot clean.data A data.frame containing the cleaned data after preprocessing.
#' @slot info.data A data.frame containing additional information related to the data.
#' @slot scale.data A data.frame containing the scaled version of the data.
#' @slot meta.featurename A character vector containing metadata feature names.
#' @slot group_col The grouping column in the data, which can be of any type.
#' @slot baseline.table A table that stores baseline information, can be of any type.
#' @slot process.info A list containing information about the processing steps.
#' @slot compute.descriptive A list containing results of descriptive statistics computations.
#' @slot corr.result A list containing the results of correlation analysis.
#' @slot var.result A list containing the results of variance analysis.
#' @slot variable.types A list specifying the types of the variables in the dataset.
#'
#' @returns An object of class 'Stat' containing the slots mentioned above.
#' @export
Stat <- setClass(
  Class = "Stat",
  slots = c(
    raw.data = "data.frame",
    clean.data = "data.frame",
    info.data = "data.frame",
    scale.data = "data.frame",
    meta.featurename = "character",
    group_col = "ANY",
    baseline.table = "ANY",
    process.info = "list",
    compute.descriptive = "list",
    corr.result = "list",
    var.result = "list",
    variable.types = "list"
  ),
  prototype = list(
    raw.data = data.frame(),
    clean.data = data.frame(),
    info.data = data.frame(),
    scale.data = data.frame(),
    meta.featurename = character(0),  
    group_col = NULL,
    baseline.table = NULL,
    process.info = list(),
    compute.descriptive = list(),
    corr.result = list(),
    var.result = list(),
    variable.types = list()
  ),
  validity = function(object) {
    # Check data frame slots
    if (!is.data.frame(object@raw.data)) {
      return("Slot 'raw.data' must be a data.frame")
    }
    if (!is.data.frame(object@clean.data)) {
      return("Slot 'clean.data' must be a data.frame")
    }
    if (!is.data.frame(object@info.data)) {
      return("Slot 'info.data' must be a data.frame")
    }
    if (!is.data.frame(object@scale.data)) {
      return("Slot 'scale.data' must be a data.frame")
    }
    
    # Check character slot
    if (!is.character(object@meta.featurename)) {
      return("Slot 'meta.featurename' must be a character vector")
    }
    
    # Check list slots
    if (!is.list(object@process.info)) {
      return("Slot 'process.info' must be a list")
    }
    if (!is.list(object@compute.descriptive)) {
      return("Slot 'compute.descriptive' must be a list")
    }
    if (!is.list(object@corr.result)) {
      return("Slot 'corr.result' must be a list")
    }
    if (!is.list(object@var.result)) {
      return("Slot 'var.result' must be a list")
    }
    if (!is.list(object@variable.types)) {
      return("Slot 'variable.types' must be a list")
    }
    
    # Check row name consistency if data exists
    if (nrow(object@clean.data) > 0 && nrow(object@info.data) > 0) {
      if (!setequal(rownames(object@clean.data), rownames(object@info.data))) {
        return("Row names of 'clean.data' and 'info.data' must match")
      }
    }
    
    return(TRUE)
  }
)

#' Create a Stat object for statistical analysis
#'
#' This function creates a `Stat` object, which is used to store and manage various components
#' of statistical analysis, such as raw data, cleaned data, additional metadata, and processing
#' information. It performs basic checks and preparation on the input data before creating the object.
#' 
#' @import methods
#' @param raw.data A data.frame containing the raw data for analysis. Defaults to an empty data frame.
#' @param clean.data A data.frame containing the cleaned data for analysis. Defaults to an empty data frame.
#' @param info.data A data.frame containing additional metadata related to the data. Defaults to an empty data frame.
#' @param group_col A character string specifying the name of the column used for grouping. Default is `"group"`.
#' @param ... Additional arguments passed to methods (not used in this function).
#' @param na.action Character string specifying how to handle NA values. 
#'   Options are "allow" (keep NA values, default), "omit" (remove rows with NA),
#'   or "error" (stop if NA values are found).
#' @returns An object of class `Stat`, which contains the processed data and metadata as slots.
#' @export
#'
#' @examples
#' raw <- data.frame(id = 1:10, value = rnorm(10), group = rep(c("A","B"), each=5))
#' inf<-data.frame(age = 1:10, value = rnorm(10), group = rep(c("A","B"), each=5))
#' stat_new <- CreateStatObject(raw.data = raw,  group_col = "group")
#' print(stat_new)
#' stat_new <- CreateStatObject(raw.data = raw, clean.data = raw, group_col = "group")
#' print(stat_new)
#' stat_new <- CreateStatObject(raw.data = raw, group_col = "group",info.data = inf)
#' print(stat_new)
CreateStatObject <- function(
    raw.data = data.frame(),
    clean.data = data.frame(),
    info.data = data.frame(),
    group_col = "group",
    na.action = c("allow", "omit", "error"),
    ...
) {

  na.action <- match.arg(na.action)
  raw.data <- as.data.frame(raw.data)
  clean.data <- as.data.frame(clean.data)
  info.data <- as.data.frame(info.data)
  
  # Store original group_col before any column name modifications
  original_group_col <- group_col
  
  # Handle NA values based on na.action
  handle_na <- function(data, data_name) {
    if (is.null(data) || nrow(data) == 0) return(data)
    if (!any(is.na(data))) return(data)
    
    switch(na.action,
           "allow" = {
             message("Note: NA values found in ", data_name, " but kept (na.action='allow')")
             return(data)
           },
           "omit" = {
             warning("NA values found in ", data_name, ". Removing ", sum(is.na(data)), " NA values.")
             return(na.omit(data))
           },
           "error" = {
             stop("NA values found in ", data_name, ". Set na.action='omit' to remove or na.action='allow' to keep.")
           }
    )
  }
  
  raw.data <- handle_na(raw.data, "raw.data")
  clean.data <- handle_na(clean.data, "clean.data")

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

  if (nrow(raw.data) == 0 && nrow(clean.data) == 0) {
    stop("At least one of 'raw.data' or 'clean.data' must be provided and not empty.")
  }


  if(nrow(raw.data) > 0) {
    raw.data <- clean_symbol_values(raw.data)
  }
  if(nrow(clean.data) > 0) {
    clean.data <- clean_symbol_values(clean.data)
  }

  if (nrow(info.data) > 0) {
    if (nrow(raw.data) > 0 && !setequal(rownames(info.data), rownames(raw.data))) {
      stop("Row names in 'info.data' do not match 'raw.data' (content mismatch).")
    }
    
    if (nrow(clean.data) > 0 && !setequal(rownames(info.data), rownames(clean.data))) {
      stop("Row names in 'info.data' do not match 'clean.data' (content mismatch).")
    }
  }
  
  # Check group_col BEFORE modifying column names
  if (!is.null(original_group_col)) {
    # Check in raw.data or clean.data (whichever has data)
    ref_data <- if(nrow(clean.data) > 0) clean.data else if(nrow(raw.data) > 0) raw.data else NULL
    if (!is.null(ref_data) && !(original_group_col %in% colnames(ref_data))) {
      stop(sprintf("Specified group_col '%s' not found in data columns: %s",
                   original_group_col, paste(colnames(ref_data), collapse = ", ")))
    }
  }

  prepare_data <- function(data, data_name) {
    if (nrow(data) == 0) {
      return(data)
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

    data <- modify_column_names(data)

    cat(paste("Data prepared for", data_name, "\n"))

    return(data)
  }

  raw.data <- prepare_data(raw.data, "raw.data")
  clean.data <- prepare_data(clean.data, "clean.data")
  
  # Update group_col if it was modified by modify_column_names
  if (!is.null(original_group_col)) {
    modified_group_col <- make.names(tolower(gsub(" ", "_", gsub("[^[:alnum:]_]", "", original_group_col))))
    if (modified_group_col != original_group_col) {
      if (modified_group_col %in% colnames(clean.data)) {
        message("Note: group_col '", original_group_col, "' was modified to '", modified_group_col, "'")
        group_col <- modified_group_col
      } else if (nrow(clean.data) > 0 && !(original_group_col %in% colnames(clean.data))) {
        warning("group_col '", original_group_col, "' not found in data columns after name modification")
      }
    }
  }

  if (nrow(raw.data) > 0) {
    meta.featurename <- as.character(colnames(raw.data))
  } else if (nrow(clean.data) > 0) {
    meta.featurename <- as.character(colnames(clean.data))
  } else {
    meta.featurename <- character()
  }

  if (nrow(info.data) == 0) {
    if (nrow(clean.data) > 0) {
      info.data <- data.frame(row.names = row.names(clean.data))
    } else if (nrow(raw.data) > 0) {
      info.data <- data.frame(row.names = row.names(raw.data))
    }
  } else {
    if (nrow(clean.data) > 0) {
      rownames(info.data) <- rownames(clean.data)
    } else if (nrow(raw.data) > 0) {
      rownames(info.data) <- rownames(raw.data)
    }
  }

  cat("Final info.data prepared.\n")

  Stat <- new(
    Class = 'Stat',
    raw.data = raw.data,
    info.data = info.data,
    clean.data = clean.data,
    meta.featurename = meta.featurename,
    process.info = list(),
    group_col = group_col
  )

  cat("Stat object created.\n")

  return(Stat)
}
