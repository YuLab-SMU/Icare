#' Detect and Mark Outliers
#'
#' @param data A data frame.
#' @param method Detection method ("zscore" or "iqr").
#' @param threshold Threshold (e.g., 3 for zscore, 1.5 for IQR).
#' @param group_col Group column.
#' @param max_unique_values Max unique values for categorical variables.
#' @param save_dir Directory to save results.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
#' @examples
#' data(iris)
#' detect_and_mark_outliers(iris[, 1:4], method = "zscore", threshold = 3)
detect_and_mark_outliers <- function(data,
                                     method = "zscore",
                                     threshold = 3,
                                     group_col = "group",
                                     max_unique_values = 5,
                                     save_dir = NULL,
                                     save_data = FALSE,
                                     csv_filename = "outlier_info.csv") {
  
  if (!is.data.frame(data)) stop("Input must be a data frame.")
  
  variable_types <- diagnose_variable_type(data, group_col, max_unique_values)
  numeric_vars <- variable_types$numeric_vars
  
  outlier_info <- list()
  outlier_indices <- list()
  
  for (col in numeric_vars) {
    x <- data[[col]]
    if (method == "zscore") {
      z_scores <- abs(scale(x))
      outliers <- which(z_scores > threshold)
    } else if (method == "iqr") {
      Q1 <- quantile(x, 0.25, na.rm = TRUE)
      Q3 <- quantile(x, 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - threshold * IQR
      upper_bound <- Q3 + threshold * IQR
      outliers <- which(x < lower_bound | x > upper_bound)
    } else {
      stop("Invalid method. Choose 'zscore' or 'iqr'.")
    }
    
    if (length(outliers) > 0) {
      outlier_info[[col]] <- list(
        count = length(outliers),
        indices = outliers,
        values = x[outliers]
      )
      outlier_indices[[col]] <- outliers
    }
  }
  
  all_outliers <- sort(unique(unlist(outlier_indices)))
  
  if (save_data) {
    if (is.null(save_dir)) {
      stop("'save_dir' cannot be NULL when 'save_data' is TRUE. Please provide a valid directory path.")
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    
    outlier_summary <- data.frame(
      Variable = names(outlier_info),
      Count = sapply(outlier_info, function(x) x$count)
    )
    full_path <- file.path(save_dir, csv_filename)
    write.csv(outlier_summary, file = full_path, row.names = FALSE)
    cat("Outlier summary saved to:", full_path, "\n")
  }
  
  stats_per_var <- list()
  for (col in numeric_vars) {
    x <- data[[col]]
    if (method == "zscore") {
      stats_per_var[[col]] <- list(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
    } else if (method == "iqr") {
      Q1 <- quantile(x, 0.25, na.rm = TRUE)
      Q3 <- quantile(x, 0.75, na.rm = TRUE)
      stats_per_var[[col]] <- list(Q1 = Q1, Q3 = Q3, IQR = Q3 - Q1)
    }
  }
  
  return(list(
    outlier_info = outlier_info,
    total_outliers = length(all_outliers),
    outlier_indices = all_outliers,
    method = method,
    threshold = threshold,
    stats = stats_per_var  
  ))
}

#' Detect and Mark Outliers in Stat Object
#'
#' @param object Stat object.
#' @param method Detection method.
#' @param threshold Threshold.
#' @param group_col Group column.
#' @param max_unique_values Max unique values.
#' @param save_dir Directory to save results.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
#' @examples
#' stat <- CreateStatObject(clean.data = mtcars, group_col = "cyl")
#' stat <- stat_detect_and_mark_outliers(stat, method = "zscore", threshold = 3)
#' print(stat@process.info$outlier_detection)
stat_detect_and_mark_outliers <- function(object,
                                          method = "zscore",
                                          threshold = 3,
                                          group_col = "group",
                                          max_unique_values = 5,
                                          save_dir = NULL,
                                          save_data = FALSE,
                                          csv_filename = "outlier_info.csv") {
  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    }   
    group_col <- slot(object, "group_col")
    if (length(group_col) == 0) {
      group_col <- NULL
    }
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  outlier_result <- detect_and_mark_outliers(
    data,
    method = method,
    threshold = threshold,
    group_col = group_col,
    max_unique_values = max_unique_values,
    save_dir =save_dir,
    save_data = save_data,
    csv_filename = csv_filename
  )
  
  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@process.info[["outlier_detection"]] <- outlier_result
    cat("- 'process.info' slot updated.\n")
    return(object)
  } else {
    return(outlier_result)
  }
}

#' Extract Outlier Data
#'
#' @param object Stat object.
#' @export
#' @examples
#' stat <- CreateStatObject(raw.data = mtcars, group_col = "cyl")
#' stat <- stat_detect_and_mark_outliers(stat)
#' out <- extract_outlier_data(stat)
#' print(out)
extract_outlier_data <- function(object) {
  if (inherits(object, "Stat")) {
    outlier_info <- object@process.info[["outlier_detection"]]
    if (is.null(outlier_info)) {
      cat("No outlier detection information found in 'process.info'.\n")
      return(NULL)
    }
    return(outlier_info)
  } else {
    stop("Input must be an object of class 'Stat'.")
  }
}

#' Handle Outliers
#'
#' @param data A data frame.
#' @param outlier_info Outlier information list.
#' @param method Handling method ("remove", "winsorize", "impute").
#' @param impute_value Value to impute (e.g., median).
#' @param save_dir Directory to save cleaned data.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
#' @examples
#' data(iris)
#' out_info <- detect_and_mark_outliers(iris[, 1:4], method = "zscore", threshold = 3)
#' cleaned <- handle_outliers(iris[, 1:4], out_info, method = "remove", save_data = FALSE)
#' head(cleaned)
handle_outliers <- function(data,
                            outlier_info,
                            method = "remove",
                            impute_value = "median",
                            save_dir = NULL,
                            save_data = FALSE,
                            csv_filename = "clean_data.csv") {
  
  if (!is.data.frame(data)) stop("Input must be a data frame.")
  if (is.null(outlier_info)) stop("Outlier information is missing.")
  
  method <- match.arg(method, choices = c("remove", "winsorize", "impute"))
  
  cleaned_data <- data
  
  if (method == "remove") {
    outlier_indices <- outlier_info$outlier_indices
    if (length(outlier_indices) > 0) {
      cleaned_data <- data[-outlier_indices, ]
      cat("Removed", length(outlier_indices), "samples with outliers.\n")
    }
  } else {
    for (col in names(outlier_info$outlier_info)) {
      indices <- outlier_info$outlier_info[[col]]$indices
      
      if (method == "winsorize") {
        lower_bound <- quantile(data[[col]], 0.05, na.rm = TRUE)
        upper_bound <- quantile(data[[col]], 0.95, na.rm = TRUE)
        cleaned_data[[col]][indices] <- ifelse(data[[col]][indices] < lower_bound, lower_bound, upper_bound)
        cat("Winsorized outliers in column:", col, "\n")
      } else if (method == "impute") {
        if (impute_value == "median") {
          val <- median(data[[col]], na.rm = TRUE)
        } else if (impute_value == "mean") {
          val <- mean(data[[col]], na.rm = TRUE)
        } else {
          stop("Invalid impute_value. Choose 'median' or 'mean'.")
        }
        cleaned_data[[col]][indices] <- val
        cat("Imputed outliers in column:", col, "with", impute_value, "\n")
      }
    }
  }
  
  if (save_data) {
    if (is.null(save_dir)) {
      stop("'save_dir' cannot be NULL when 'save_data' is TRUE. Please provide a valid directory path.")
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(cleaned_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")
  }
  
  return(cleaned_data)
}

#' Handle Outliers in Stat Object
#'
#' @param object Stat object.
#' @param method Handling method.
#' @param impute_value Impute value.
#' @param save_dir Save directory.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
#' @examples
#' stat <- CreateStatObject(clean.data = mtcars, group_col = "cyl")
#' stat <- stat_detect_and_mark_outliers(stat)
#' stat <- stat_handle_outliers(stat, method = "remove", save_data = FALSE)
#' head(stat@clean.data)
stat_handle_outliers <- function(object,
                                 method = "remove",
                                 impute_value = "median",
                                 save_dir = NULL,
                                 save_data = FALSE,
                                 csv_filename = "clean_data.csv") {
  
  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    }   
    outlier_info <- object@process.info[["outlier_detection"]]
  } else {
    stop("Input must be an object of class 'Stat'.")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  if (is.null(outlier_info)) {
    stop("No outlier detection information found. Run stat_detect_and_mark_outliers first.")
  }
  
  cleaned_data <- handle_outliers(
    data,
    outlier_info,
    method = method,
    impute_value = impute_value,
    save_dir =save_dir,
    save_data = save_data,
    csv_filename = csv_filename
  )
  
  cat("Updating 'Stat' object...\n")
  object@clean.data <- cleaned_data
  object@process.info[["outlier_handling"]] <- list(
    method = method,
    impute_value = impute_value,
    detection_stats = outlier_info$stats,   
    threshold = outlier_info$threshold
  )
  cat("- 'clean.data' slot updated.\n")
  cat("- 'process.info' slot updated.\n")
  
  return(object)
}
