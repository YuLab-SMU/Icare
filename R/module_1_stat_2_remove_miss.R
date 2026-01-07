#' Remove Variables and Samples with High Missing Values
#'
#' @param data A data frame.
#' @param miss_threshold The threshold percentage of missing values.
#' @param save_dir Directory to save cleaned data.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
remove_high_missing <- function(data, 
                                miss_threshold = 25,
                                save_dir = here::here("StatObject","Data"),
                                save_data = TRUE,
                                csv_filename = "clean_data.csv") {
  if (!is.data.frame(data)) stop("Input must be a data frame.")
  if (miss_threshold < 0 || miss_threshold > 100) {
    stop("miss_threshold must be between 0 and 100.")
  }
  
  
  data[data == '<NA>' | data == 'NA' | data == '' | data == 'NULL'] <- NA
  
  var_missing_percentage <- colMeans(is.na(data)) * 100
  sample_missing_percentage <- rowMeans(is.na(data)) * 100
  
  high_missing_vars <- names(var_missing_percentage[var_missing_percentage >= miss_threshold])
  high_missing_samples <- which(sample_missing_percentage >= miss_threshold)
  
  
  
  if (length(high_missing_vars) > 0) {
    data <- data[, !names(data) %in% high_missing_vars, drop = FALSE]
  }
  
  if (length(high_missing_samples) > 0) {
    data <- data[-high_missing_samples, , drop = FALSE]
  }
  
  if (ncol(data) == 0) stop("No variables remain after removing high-missing variables.")
  if (nrow(data) == 0) stop("No samples remain after removing high-missing samples.")
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    full_path <- file.path(save_dir, csv_filename)
    write.csv(data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")
  }
  
  
  cat("\n=== Data Cleaning Summary ===\n")
  cat("Missing threshold:", miss_threshold, "%\n")
  cat("Removed variables:", length(high_missing_vars), "\n")
  if (length(high_missing_vars) > 0) {
    cat("Removed variable names:", paste(high_missing_vars, collapse = ", "), "\n")
  }
  cat("Removed samples:", length(high_missing_samples), "\n")
  cat("Final data:", nrow(data), "samples,", ncol(data), "variables\n")
  
  return(list(
    miss_threshold=miss_threshold,
    cleaned_data = data,
    high_missing_vars = high_missing_vars
  ))
}


#' Remove Missing Values from Stat Object or Data Frame
#'
#' @param object A Stat object or data frame.
#' @param miss_threshold Missing threshold.
#' @param save_dir Directory to save.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
stat_miss_remove <- function(object, 
                             miss_threshold = 25,
                             save_dir = here::here("StatObject","Data"),
                             save_data = TRUE,
                             csv_filename = "clean_data.csv") {
  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    }    
    group_col <- slot(object, "group_col")
    group_col <- if (length(group_col) == 0) NULL else group_col
  } else if (is.data.frame(object)) {
    data <- object
    group_col <- NULL
  } else {
    stop("Input must be an object of class 'Stat' or a data frame")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  clean_drop_data <- remove_high_missing(data, 
                                         miss_threshold = miss_threshold,
                                         save_dir =save_dir,
                                         save_data = save_data,
                                         csv_filename = csv_filename)
  
  
  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@clean.data <- clean_drop_data$cleaned_data
    object@process.info[["missing_removal"]]<-clean_drop_data
    cat("- 'clean.data' slot updated.\n")
    cat("- 'process.info' slot updated.\n")
    return(object)
  } else {
    return(clean_drop_data)
  }
}
