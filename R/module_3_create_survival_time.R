#' Convert Two Class to Binary
#'
#' @param data Data frame.
#' @export
convert_two_class_to_binary <- function(data) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }
  
  original_row_names <- rownames(data)
  
  data_converted <- data.frame(lapply(data, function(column) {
    
    unique_values <- unique(column)
    
    if (length(unique_values) == 2 && all(unique_values %in% c(1, 2))) {
      column <- ifelse(column == 1, 0, 1)
    }
    return(column)
  }))
  
  rownames(data_converted) <- original_row_names
  
  return(data_converted)
}

#' Convert to Numeric
#'
#' @param data Data frame.
#' @export
convert_to_numeric <- function(data) {
  data[] <- lapply(data, function(x) as.numeric(x))
  return(data)
}

#' Sub Survival Data
#'
#' @param object Subtyping object.
#' @param time_col Time col.
#' @param status_col Status col.
#' @param keep_all Keep all.
#' @export
Sub_survival_data <- function(
    object,
    time_col = "time",
    status_col = "status",
    keep_all = FALSE
) {
  if (!inherits(object, "Subtyping")) {
    stop("Input must be a 'Subtyping' class object")
  }
  
  clustered.data <- object@clustered.data
  info.data <- object@info.data
  clean.data <- object@clean.data
  
  if (!(time_col %in% colnames(info.data))) {
    stop("Time column '", time_col, "' not found in info.data")
  }
  
  if (!(status_col %in% colnames(info.data))) {
    stop("Status column '", status_col, "' not found in info.data")
  }
  
  if (!identical(rownames(clustered.data), rownames(info.data))) {
    common_rows <- intersect(rownames(clustered.data), rownames(info.data))
    if (length(common_rows) == 0) {
      stop("No matching row names between clustered.data and info.data")
    }
    clustered.data <- clustered.data[common_rows, , drop = FALSE]
    info.data <- info.data[common_rows, , drop = FALSE]
    cat("Matched", length(common_rows), "rows between clustered.data and info.data\n")
  }
  
  info.data[[status_col]] <- as.factor(info.data[[status_col]])
  cat("Creating survival dataset...\n")
  survival.data <- clustered.data
  survival.data[[time_col]] <- info.data[[time_col]]
  survival.data[[status_col]] <- info.data[[status_col]]
  
  
  if (keep_all) {
    na_rows <- !complete.cases(survival.data[, c(time_col, status_col)])
    if (any(na_rows)) {
      warning("Found ", sum(na_rows), " rows with missing survival data")
    }
  } else {
    original_rows <- nrow(survival.data)
    survival.data <- survival.data[complete.cases(survival.data[, c(time_col, status_col)]), ]
    removed_rows <- original_rows - nrow(survival.data)
    if (removed_rows > 0) {
      cat("Removed", removed_rows, "rows with missing survival data\n")
    }
  }
  cat("Standardizing column names...\n")
  colnames(survival.data)[colnames(survival.data) == time_col] <- "time"
  colnames(survival.data)[colnames(survival.data) == status_col] <- "status"
  survival.data <- convert_to_numeric(survival.data)
  survival.data <- convert_two_class_to_binary(survival.data)
  if (inherits(object, "Subtyping")) {
    object@survival.data <- survival.data
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'survival.data' slot updated.\n")
    return(object)
  } 
}
