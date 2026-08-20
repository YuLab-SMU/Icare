#' Log Transformation
#'
#' @param x Numeric vector.
#' @return Log-transformed vector.
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
min_max_scale <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

#' Z-Score Standardization
#'
#' @param x Numeric vector.
#' @return Standardized vector.
z_score_standardize <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

#' Center Data
#'
#' @param x Numeric vector.
#' @return Centered vector.
center_data <- function(x) {
  return(x - mean(x, na.rm = TRUE))
}

#' Scale Data
#'
#' @param x Numeric vector.
#' @return Scaled vector.
scale_data <- function(x) {
  return(x / sd(x, na.rm = TRUE))
}

#' Max Abs Scaling
#'
#' @param x Numeric vector.
#' @return Scaled vector.
max_abs_scale <- function(x) {
  return(x / max(abs(x), na.rm = TRUE))
}

#' Box-Cox Transformation
#'
#' @param x Numeric vector.
#' @return Transformed vector.
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
yeojohnson_transform <- function(x) {
  yj <- caret::preProcess(data.frame(x), method = "YeoJohnson")
  return(predict(yj, data.frame(x))[, 1])
}

#' Fit/Apply Helpers for Reproducible Normalization
#'
#' Each `.fit_*(x)` function learns the parameters of a normalization method
#' from a numeric vector `x` (typically a training-set column) and returns
#' them as a list. Each matching `.apply_*(x, params)` function reproduces
#' the *exact same* transform on any new vector `x` using those stored
#' parameters, instead of recomputing statistics (mean/sd/min/max/lambda...)
#' from `x` itself.
#'
#' This is what lets `process_new_data()` apply the identical normalization
#' recipe learned on the training data to new/deployment data -- including a
#' single new observation, for which recomputing mean/sd or min/max from
#' scratch would be undefined (NA/NaN). The existing `log_transform()`,
#' `min_max_scale()`, etc. functions above are left untouched; these
#' fit/apply pairs are used internally by `normalize_data()` in addition to
#' them, purely to also capture and store the parameters.
#'
#' @name normalize-fit-apply
#' @keywords internal
NULL

.fit_log <- function(x) {
  shift <- 0
  if (any(x <= 0, na.rm = TRUE)) shift <- abs(min(x, na.rm = TRUE)) + 1
  list(method = "log", shift = shift)
}
.apply_log <- function(x, params) log(x + params$shift)

.fit_min_max <- function(x) {
  list(method = "min_max", min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
}
.apply_min_max <- function(x, params) {
  rng <- params$max - params$min
  if (is.na(rng) || rng == 0) return(rep(0, length(x)))  # degenerate: constant training column
  (x - params$min) / rng
}

.fit_z_score <- function(x) {
  list(method = "z_score", mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
}
.apply_z_score <- function(x, params) {
  if (is.na(params$sd) || params$sd == 0) return(rep(0, length(x)))
  (x - params$mean) / params$sd
}

.fit_center <- function(x) list(method = "center", mean = mean(x, na.rm = TRUE))
.apply_center <- function(x, params) x - params$mean

.fit_scale <- function(x) list(method = "scale", sd = sd(x, na.rm = TRUE))
.apply_scale <- function(x, params) {
  if (is.na(params$sd) || params$sd == 0) return(x)
  x / params$sd
}

.fit_max_abs <- function(x) list(method = "max_abs", max_abs = max(abs(x), na.rm = TRUE))
.apply_max_abs <- function(x, params) {
  if (is.na(params$max_abs) || params$max_abs == 0) return(x)
  x / params$max_abs
}

.fit_box_cox <- function(x) {
  shift <- 0
  if (any(x <= 0, na.rm = TRUE)) shift <- abs(min(x, na.rm = TRUE)) + 1
  bc <- caret::BoxCoxTrans(x + shift)
  list(method = "box_cox", shift = shift, transform = bc)
}
.apply_box_cox <- function(x, params) predict(params$transform, x + params$shift)

.fit_yeo_johnson <- function(x) {
  yj <- caret::preProcess(data.frame(x = x), method = "YeoJohnson")
  list(method = "yeo_johnson", transform = yj)
}
.apply_yeo_johnson <- function(x, params) predict(params$transform, data.frame(x = x))[, 1]

#' Preprocess Data
#'
#' @param data Data frame.
#' @param method Preprocessing method.
#' @param group_col Group column.
#' @param max_unique_values Max unique values.
#' @return Processed data frame.
#' @export
#' @examples
#' data(iris)
#' iris_num <- iris[, 1:4]
#' proc <- preprocess_data(iris_num, method = "z_score")
#' head(proc)
preprocess_data <- function(data, method = "log", group_col = "group", max_unique_values = 5) {
  variable_types <- diagnose_variable_type(
    data = data,
    group_col = group_col,
    max_unique_values = max_unique_values,
    treat_low_card_numeric_as_categorical = TRUE   
  )
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

#' Normalize Data
#'
#' This function applies normalization to numeric variables in a data frame.
#' Supports both automatic method selection and user-specified methods.
#' Normalization parameters are stored as attributes for later reproduction.
#'
#' @param data A data frame containing the data to normalize.
#' @param method Normalization method. Either \code{"auto"} for automatic
#'   selection based on data characteristics, or one of the specific methods:
#'   \code{"log"}, \code{"min_max"}, \code{"z_score"}, \code{"center"},
#'   \code{"scale"}, \code{"max_abs"}, \code{"box_cox"}, \code{"yeo_johnson"}.
#'   Default is \code{"auto"}.
#' @param group_col Character string specifying the grouping column name.
#'   This column is excluded from normalization. Default is \code{"group"}.
#' @param max_unique_values Integer. Maximum number of unique values for a
#'   variable to be considered categorical. Default is \code{5}.
#' @param save_dir Character string specifying the directory to save the
#'   normalized data. If \code{NULL}, uses the default output directory.
#'   Default is \code{NULL}.
#' @param save_data Logical. If \code{TRUE}, writes the normalized data to
#'   a CSV file. Default is \code{TRUE}.
#' @param csv_filename Character string. Name of the output CSV file.
#'   Default is \code{"scale_data.csv"}.
#'
#' @return A normalized data frame with two attributes:
#'   \itemize{
#'     \item \code{normalization_info}: Named list of methods used per column.
#'     \item \code{normalization_params}: Named list of fitted parameters
#'       for each column, enabling reproduction on new data.
#'   }
#'
#' @details
#' When \code{method = "auto"}, the function selects a normalization method
#' based on the following criteria:
#' \itemize{
#'   \item Constant or near-constant columns are skipped (method = "none").
#'   \item For large samples (>5000), Anderson-Darling test is used for
#'     normality assessment; otherwise Shapiro-Wilk test is used.
#'   \item If the data is approximately normal, \code{"z_score"} is selected.
#'   \item If skewness > 1 or < -1, \code{"box_cox"} (for positive data) or
#'     \code{"yeo_johnson"} (for data with non-positive values) is selected.
#'   \item Otherwise, \code{"min_max"} scaling is applied.
#' }
#'
#' The function stores the fitted parameters (mean, sd, min, max, lambda, shift)
#' as attributes of the returned data frame. These parameters can be retrieved
#' and applied to new data using \code{process_new_data()} or
#' \code{Sub_apply_norm_params()}.
#'
#' @export
#' @examples
#' data(iris)
#' iris_num <- iris[, 1:4]
#' norm <- normalize_data(iris_num, method = "min_max")
#' summary(norm)
normalize_data <- function(data,
                           method = "auto",
                           group_col = "group",
                           max_unique_values = 5,
                           save_dir = NULL,
                           save_data = TRUE,
                           csv_filename = "scale_data.csv") {
  
  # Set default save directory if not provided
  if (is.null(save_dir)) save_dir <- get_output_dir("StatObject", "Data")
  
  # Validate input
  if (!is.data.frame(data)) stop("Input must be a data frame.")
  
  # Identify variable types (numeric vs categorical)
  variable_types <- diagnose_variable_type(
    data = data,
    group_col = group_col,
    max_unique_values = max_unique_values,
    treat_low_card_numeric_as_categorical = TRUE
  )
  numeric_vars <- variable_types$numeric_vars
  
  # Return early if no numeric variables found
  if (length(numeric_vars) == 0) {
    warning("No numeric variables found for normalization.")
    return(data)
  }
  
  # Initialize containers for normalization metadata
  normalized_data <- data
  normalization_info   <- list()   # column -> method name
  normalization_params <- list()   # column -> fitted parameters
  
  # Loop through each numeric variable
  for (col in numeric_vars) {
    x <- data[[col]]
    
    # ---- Determine normalization method ----
    if (method == "auto") {
      # Automatic method selection based on data characteristics
      
      # Skip constant or near-constant columns
      if (length(unique(na.omit(x))) < 2 || sd(x, na.rm = TRUE) == 0) {
        selected_method <- "none"   
        warning("Column '", col, "' is constant or near-constant; normalization skipped.")
      } else {
        # Normality test: Anderson-Darling for large samples, Shapiro-Wilk otherwise
        n_nonmiss <- sum(!is.na(x))
        if (n_nonmiss > 5000) {
          if (requireNamespace("nortest", quietly = TRUE)) {
            p_val <- tryCatch(nortest::ad.test(x)$p.value, error = function(e) 0)
          } else {
            warning("Package 'nortest' not available; using skewness only for large samples.")
            p_val <- 0 
          }
        } else {
          p_val <- tryCatch(shapiro.test(x)$p.value, error = function(e) 0)
        }
        
        # Calculate skewness
        skewness_val <- e1071::skewness(x, na.rm = TRUE)
        
        # Select method based on normality and skewness
        if (p_val > 0.05) {
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
      }
    } else {
      # User-specified method
      selected_method <- method
    }
    
    # ---- Apply normalization or skip ----
    if (selected_method == "none") {
      normalized_data[[col]] <- x
      next
    }
    
    # ---- Fit parameters and apply transformation ----
    # The .fit_* functions learn parameters from the data; .apply_* functions
    # apply the transformation using those parameters. This enables exact
    # reproduction on new data (e.g., for prediction/deployment).
    params <- switch(selected_method,
                     "log"         = .fit_log(x),
                     "min_max"     = .fit_min_max(x),
                     "z_score"     = .fit_z_score(x),
                     "center"      = .fit_center(x),
                     "scale"       = .fit_scale(x),
                     "max_abs"     = .fit_max_abs(x),
                     "box_cox"     = .fit_box_cox(x),
                     "yeo_johnson" = .fit_yeo_johnson(x),
                     stop("Invalid method.")
    )
    
    normalized_data[[col]] <- switch(selected_method,
                                     "log"         = .apply_log(x, params),
                                     "min_max"     = .apply_min_max(x, params),
                                     "z_score"     = .apply_z_score(x, params),
                                     "center"      = .apply_center(x, params),
                                     "scale"       = .apply_scale(x, params),
                                     "max_abs"     = .apply_max_abs(x, params),
                                     "box_cox"     = .apply_box_cox(x, params),
                                     "yeo_johnson" = .apply_yeo_johnson(x, params)
    )
    
    # Store metadata for reproducibility
    normalization_info[[col]]   <- selected_method
    normalization_params[[col]] <- params
  }
  
  # Save normalized data to CSV if requested
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    full_path <- file.path(save_dir, csv_filename)
    write.csv(normalized_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")
  }
  
  # Attach normalization metadata as attributes for downstream use
  attr(normalized_data, "normalization_info")   <- normalization_info
  attr(normalized_data, "normalization_params") <- normalization_params
  
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
#' @examples
#' stat <- CreateStatObject(clean.data = mtcars, group_col = "cyl")
#' stat <- stat_normalize_process(stat, method = "z_score")
#' head(stat@scale.data)
stat_normalize_process <- function(object,
                                   method = "auto",
                                   group_col = "group",
                                   max_unique_values = 5,
                                   save_dir = NULL,
                                   save_data = TRUE,
                                   csv_filename = "scale_data.csv") {
  if (is.null(save_dir)) save_dir <- get_output_dir("StatObject", "Data")
  
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
    save_dir = save_dir,
    save_data = save_data,
    csv_filename = csv_filename
  )
  
  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@scale.data <- normalized_data
    object@process.info[["normalization"]]        <- attr(normalized_data, "normalization_info")
    object@process.info[["normalization_params"]] <- attr(normalized_data, "normalization_params")
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
#' @examples
#' stat <- CreateStatObject(clean.data = mtcars, group_col = "cyl")
#' stat <- stat_normalize_process(stat, method = "z_score")
#' scaled <- ExtractScaleData(stat)
#' head(scaled)
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
#' @examples
#' stat <- CreateStatObject(clean.data = mtcars, group_col = "cyl")
#' clean <- ExtractCleanData(stat)
#' head(clean)
ExtractCleanData <- function(object) {
  if (inherits(object, "Stat")) {
    return(object@clean.data)
  } else {
    stop("Input must be an object of class 'Stat'.")
  }
}
