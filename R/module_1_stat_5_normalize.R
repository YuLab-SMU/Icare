#' Log Transformation
#'
#' @param x Numeric vector.
#' @return Log-transformed vector.
#' @export
log_transform <- function(x) {
  if (any(x <= 0, na.rm = TRUE)) {
    warning("Data contains non-positive values. Adding constant before log transformation.")
    x <- x + abs(min(x, na.rm = TRUE)) + 1
  }
  return(log(x))
}

#' Min-Max Scaling
#'
#' @param x Numeric vector.
#' @return Scaled vector.
#' @export
min_max_scale <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

#' Z-Score Standardization
#'
#' @param x Numeric vector.
#' @return Standardized vector.
#' @export
z_score_standardize <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

#' Center Data
#'
#' @param x Numeric vector.
#' @return Centered vector.
#' @export
center_data <- function(x) {
  return(x - mean(x, na.rm = TRUE))
}

#' Scale Data
#'
#' @param x Numeric vector.
#' @return Scaled vector.
#' @export
scale_data <- function(x) {
  return(x / sd(x, na.rm = TRUE))
}

#' Max Abs Scaling
#'
#' @param x Numeric vector.
#' @return Scaled vector.
#' @export
max_abs_scale <- function(x) {
  return(x / max(abs(x), na.rm = TRUE))
}

#' Box-Cox Transformation
#'
#' @param x Numeric vector.
#' @return Transformed vector.
#' @export
boxcox_transform <- function(x) {
  if (any(x <= 0, na.rm = TRUE)) {
    warning("Data contains non-positive values. Adding constant before Box-Cox transformation.")
    x <- x + abs(min(x, na.rm = TRUE)) + 1
  }
  bc <- caret::BoxCoxTrans(x)
  return(predict(bc, x))
}

#' Yeo-Johnson Transformation
#'
#' @param x Numeric vector.
#' @return Transformed vector.
#' @export
yeojohnson_transform <- function(x) {
  yj <- caret::preProcess(data.frame(x), method = "YeoJohnson")
  return(predict(yj, data.frame(x))[, 1])
}

#' Preprocess Data
#'
#' @param data Data frame.
#' @param method Preprocessing method.
#' @param group_col Group column.
#' @param max_unique_values Max unique values.
#' @return Processed data frame.
#' @export
preprocess_data <- function(data, method = "log", group_col = "group", max_unique_values = 5) {
  variable_types <- diagnose_variable_type(data, group_col, max_unique_values)
  numeric_vars <- variable_types$numeric_vars
  
  processed_data <- data
  
  for (col in numeric_vars) {
    x <- data[[col]]
    if (method == "log") {
      processed_data[[col]] <- log_transform(x)
    } else if (method == "min_max") {
      processed_data[[col]] <- min_max_scale(x)
    } else if (method == "z_score") {
      processed_data[[col]] <- z_score_standardize(x)
    } else if (method == "center") {
      processed_data[[col]] <- center_data(x)
    } else if (method == "scale") {
      processed_data[[col]] <- scale_data(x)
    } else if (method == "max_abs") {
      processed_data[[col]] <- max_abs_scale(x)
    } else if (method == "box_cox") {
      processed_data[[col]] <- boxcox_transform(x)
    } else if (method == "yeo_johnson") {
      processed_data[[col]] <- yeojohnson_transform(x)
    } else {
      stop("Invalid method.")
    }
  }
  
  return(processed_data)
}

#' Normalize Data (Auto Selection)
#'
#' @param data Data frame.
#' @param method Normalization method ("auto" or specific).
#' @param group_col Group column.
#' @param max_unique_values Max unique values.
#' @param save_dir Save directory.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @return Normalized data frame.
#' @export
normalize_data <- function(data,
                           method = "auto",
                           group_col = "group",
                           max_unique_values = 5,
                           save_dir = here::here("StatObject","Data"),
                           save_data = TRUE,
                           csv_filename = "scale_data.csv") {
  
  if (!is.data.frame(data)) stop("Input must be a data frame.")
  
  variable_types <- diagnose_variable_type(data, group_col, max_unique_values)
  numeric_vars <- variable_types$numeric_vars
  
  if (length(numeric_vars) == 0) {
    warning("No numeric variables found for normalization.")
    return(data)
  }
  
  normalized_data <- data
  normalization_info <- list()
  
  for (col in numeric_vars) {
    x <- data[[col]]
    
    if (method == "auto") {
      shapiro_test <- tryCatch(shapiro.test(x)$p.value, error = function(e) 0)
      skewness_val <- e1071::skewness(x, na.rm = TRUE)
      
      if (shapiro_test > 0.05) {
        selected_method <- "z_score"
      } else if (abs(skewness_val) > 1) {
        if (all(x > 0, na.rm = TRUE)) {
          selected_method <- "box_cox"
        } else {
          selected_method <- "yeo_johnson"
        }
      } else {
        selected_method <- "min_max"
      }
    } else {
      selected_method <- method
    }
    
    normalized_data[[col]] <- switch(selected_method,
                                     "log" = log_transform(x),
                                     "min_max" = min_max_scale(x),
                                     "z_score" = z_score_standardize(x),
                                     "center" = center_data(x),
                                     "scale" = scale_data(x),
                                     "max_abs" = max_abs_scale(x),
                                     "box_cox" = boxcox_transform(x),
                                     "yeo_johnson" = yeojohnson_transform(x),
                                     stop("Invalid method.")
    )
    
    normalization_info[[col]] <- selected_method
  }
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(normalized_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")
  }
  
  attr(normalized_data, "normalization_info") <- normalization_info
  return(normalized_data)
}

#' Normalize Data in Stat Object
#'
#' @param object Stat object.
#' @param method Normalization method.
#' @param group_col Group column.
#' @param max_unique_values Max unique values.
#' @param save_dir Save directory.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
stat_normalize_process <- function(object,
                                   method = "auto",
                                   group_col = "group",
                                   max_unique_values = 5,
                                   save_dir = here::here("StatObject","Data"),
                                   save_data = TRUE,
                                   csv_filename = "scale_data.csv") {
  
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
  
  normalized_data <- normalize_data(
    data,
    method = method,
    group_col = group_col,
    max_unique_values = max_unique_values,
    save_dir =save_dir,
    save_data = save_data,
    csv_filename = csv_filename
  )
  
  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@scale.data <- normalized_data
    object@process.info[["normalization"]] <- attr(normalized_data, "normalization_info")
    cat("- 'scale.data' slot updated.\n")
    cat("- 'process.info' slot updated.\n")
    return(object)
  } else {
    return(normalized_data)
  }
}

#' Extract Scaled Data
#'
#' @param object Stat object.
#' @export
ExtractScaleData <- function(object) {
  if (inherits(object, "Stat")) {
    return(object@scale.data)
  } else {
    stop("Input must be an object of class 'Stat'.")
  }
}

#' Extract Clean Data
#'
#' @param object Stat object.
#' @export
ExtractCleanData <- function(object) {
  if (inherits(object, "Stat")) {
    return(object@clean.data)
  } else {
    stop("Input must be an object of class 'Stat'.")
  }
}
