#' Logarithmic Transformation of Data
#'
#' Applies a logarithmic transformation to the input data with a small offset to handle zeros.
#' This transformation is commonly used to stabilize variance and make the data more normally distributed.
#' @import stats
#' @param x A numeric vector or matrix containing non-negative values to be transformed.
#'          Negative values will cause an error.
#'
#' @returns A list containing three components:
#'   \item{scaled_data}{The log-transformed data (log(x + 1e-8))}
#'   \item{normalize_method}{A character string indicating the transformation method ("log_transform")}
#'   \item{normalize_info}{A list containing transformation details including the offset value used (1e-8)}
#'
#' @export
#'
#' @examples
#' # Basic usage with positive values
#' x <- c(1, 10, 100, 1000)
#' result <- log_transform(x)
#' print(result$scaled_data)
#'
#' # Handling zeros
#' y <- c(0, 1, 10, 100)
#' result <- log_transform(y)
#'
#' \dontrun{
#' # Will throw an error due to negative values
#' z <- c(-1, 0, 1, 2)
#' log_transform(z)
#' }
log_transform <- function(x) {
  if (any(x < 0, na.rm = TRUE)) {
    stop("Input contains negative values. Log transformation is not possible.")
  }
  scaled_data <- log(x + 1e-8)
  return(list(scaled_data = scaled_data,
              normalize_method = "log_transform",
              normalize_info = list(offset = 1e-8)))
}

#' Min-Max Normalization
#'
#' Scales numeric data to [0,1] range. Handles constant columns by adding small noise.
#'
#' @import stats
#' @param x Numeric vector, matrix or data frame to be scaled.
#' @return List with: scaled data (data.frame), method name, and scaling parameters.
#'
#' @details Uses formula: (x - min(x))/(max(x) - min(x)). Adds noise (SD=1e-8) to
#' constant columns to avoid division by zero.
#'
#' @export
#' @examples
#' # Basic usage
#' min_max_scale(data.frame(a=1:5, b=10:14))
#'
#' # With constant column (triggers warning)
#' min_max_scale(data.frame(a=c(1,1,1), b=c(2,4,6)))
min_max_scale <- function(x) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  min_vals <- apply(x, 2, min, na.rm = TRUE)
  max_vals <- apply(x, 2, max, na.rm = TRUE)
  
  constant_cols <- which(max_vals == min_vals)
  if (length(constant_cols) > 0) {
    warning(paste("Columns", paste(names(constant_cols), collapse = ", "),
                  "are constant. Adding small noise."))
    for (col in names(constant_cols)) {
      x[[col]] <- x[[col]] + rnorm(length(x[[col]]), sd = 1e-8)
    }
    min_vals <- apply(x, 2, min, na.rm = TRUE)
    max_vals <- apply(x, 2, max, na.rm = TRUE)
  }
  
  scaled_data <- as.data.frame(sweep(sweep(x, 2, min_vals, "-"), 2, (max_vals - min_vals), "/"))
  return(list(scaled_data = scaled_data,
              normalize_method = "min_max_scale",
              normalize_info = list(min_vals = min_vals,
                                    max_vals = max_vals)))
}





#' Z-Score Standardization
#'
#' Standardizes numeric data to have mean=0 and standard deviation=1. Handles constant
#' columns by setting their standard deviation to 1.
#'
#' @import stats
#' @param x Numeric vector, matrix or data frame to be standardized.
#' @return List containing:
#' \item{scaled_data}{Standardized data (data.frame)}
#' \item{normalize_method}{Method name ("z_score_standardize")}
#' \item{normalize_info}{List with mean and standard deviation values for each column}
#'
#' @details Uses formula: (x - mean(x))/sd(x). For constant columns (sd=0),
#' sets sd to 1 to avoid division by zero.
#'
#' @export
#' @examples
#' # Basic usage
#' z_score_standardize(data.frame(a=1:5, b=c(10,20,30,40,50)))
#'
#' # With constant column
#' z_score_standardize(data.frame(a=c(1,1,1), b=c(2,4,6)))
z_score_standardize <- function(x) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  mean_vals <- apply(x, 2, mean, na.rm = TRUE)
  sd_vals <- apply(x, 2, sd, na.rm = TRUE)
  sd_vals[sd_vals == 0] <- 1
  
  standardized_data <- as.data.frame(sweep(sweep(x, 2, mean_vals, "-"), 2, sd_vals, "/"))
  return(list(scaled_data = standardized_data,
              normalize_method = "z_score_standardize",
              normalize_info = list(mean_vals = mean_vals,
                                    sd_vals = sd_vals)))
}

#' Data Centering
#'
#' Centers numeric data by subtracting the mean from each value (mean=0).
#'
#' @param x Numeric vector, matrix or data frame to be centered. Non-numeric columns
#' will be converted to numeric if possible.
#' @import stats
#' @return List containing:
#' \item{scaled_data}{Centered data (data.frame)}
#' \item{normalize_method}{Method name ("center_data")}
#' \item{normalize_info}{List with mean values used for centering}
#'
#' @details Performs centering using the formula: x - mean(x). This transformation
#' shifts the data distribution to have mean zero while preserving the original scale.
#'
#' @export
#' @examples
#' # Basic usage
#' center_data(data.frame(a = 1:5, b = c(10, 20, 30, 40, 50)))
#'
#' # With matrix input
#' center_data(matrix(1:10, ncol = 2))
center_data <- function(x) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  mean_vals <- apply(x, 2, mean, na.rm = TRUE)
  centered_data <- sweep(x, 2, mean_vals, "-")
  return(list(scaled_data = centered_data,
              normalize_method = "center_data",
              normalize_info = list(center = mean_vals)))
}

#' Data Scaling
#'
#' Scales numeric data to unit variance (standard deviation=1). Handles constant
#' columns by setting their standard deviation to 1.
#' @import stats
#' @param x Numeric vector, matrix or data frame to be scaled. Non-numeric columns
#' will be converted to numeric if possible.
#'
#' @return List containing:
#' \item{scaled_data}{Scaled data (data.frame)}
#' \item{normalize_method}{Method name ("scale_data")}
#' \item{normalize_info}{List with standard deviation values used for scaling}
#'
#' @details Performs scaling using the formula: x / sd(x). For constant columns (sd=0),
#' sets sd to 1 to avoid division by zero. This transformation preserves the mean while
#' adjusting the scale of the data.
#'
#' @export
#' @examples
#' # Basic usage
#' scale_data(data.frame(a = 1:5, b = c(1, 2, 3, 4, 5)))
#'
#' # With constant column
#' scale_data(data.frame(a = c(1,1,1), b = c(2,4,6)))
scale_data <- function(x) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  scale_vals <- apply(x, 2, sd, na.rm = TRUE)
  scale_vals[scale_vals == 0] <- 1
  scaled_data <- sweep(x, 2, scale_vals, "/")
  return(list(scaled_data = scaled_data,
              normalize_method = "scale_data",
              normalize_info = list(scale = scale_vals)))
}

#' Maximum Absolute Value Scaling
#'
#' Scales numeric data by dividing each feature by its maximum absolute value, resulting in
#' values between -1 and 1. Handles constant columns by setting their scale factor to 1.
#'
#'
#' @param x Numeric vector, matrix or data frame to be scaled. Non-numeric columns will be
#' automatically converted to numeric if possible.
#'
#' @import stats
#' @return A list containing:
#' \item{scaled_data}{The scaled data (data.frame) with values between -1 and 1}
#' \item{normalize_method}{Character string indicating the scaling method ("max_abs_scale")}
#' \item{normalize_info}{List containing the maximum absolute values used for scaling}
#'
#' @details
#' This scaling method uses the formula: x / max(abs(x)). It preserves the sign of the original
#' data while scaling all values to the range [-1, 1]. For constant zero columns (max_abs = 0),
#' the function sets the scale factor to 1 to avoid division by zero.
#'
#' @export
#' @examples
#' # Basic usage with positive and negative values
#' max_abs_scale(data.frame(a = -5:5, b = seq(-10, 10, by = 2)))
#'
#' # With constant column (will be handled automatically)
#' max_abs_scale(data.frame(a = c(0, 0, 0), b = c(0.5, -0.5, 0.25)))
#'
#' # With matrix input
#' max_abs_scale(matrix(c(-3, -2, -1, 0, 1, 2), ncol = 2))
max_abs_scale <- function(x) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  max_abs <- apply(x, 2, function(col) max(abs(col), na.rm = TRUE))
  max_abs[max_abs == 0] <- 1
  scaled_data <- sweep(x, 2, max_abs, "/")
  return(list(scaled_data = scaled_data,
              normalize_method = "max_abs_scale",
              normalize_info = list(max_abs = max_abs)))
}

#' Box-Cox Transformation for Positive Data
#'
#' Applies a Box-Cox transformation to a numeric vector to stabilize variance and make the data more normally distributed.
#' If the Box-Cox transformation fails or the optimal lambda is close to 0, a log transformation is applied instead.
#'
#' @importFrom MASS boxcox
#' @import stats
#' @param x A numeric vector. Values must be positive and non-missing for transformation. Missing or non-positive values will be ignored in lambda estimation.
#'
#' @returns A list containing:
#' \describe{
#'   \item{scaled_data}{A numeric vector of transformed data, with the same length as \code{x}.}
#'   \item{normalize_method}{A character string indicating the method used: \code{"boxcox_transform"}.}
#'   \item{normalize_info}{A list containing the lambda used in the Box-Cox transformation.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- abs(rnorm(100))
#' result <- boxcox_transform(x)
boxcox_transform <- function(x) {
  
  x_clean <- x[!is.na(x) & x > 0]
  if (length(x_clean) < 3) {
    stop("Insufficient positive values for Box-Cox transform")
  }
  
  
  bc <- tryCatch(
    boxcox(lm(x_clean ~ 1)),
    error = function(e) {
      warning("Box-Cox failed, using log transform instead")
      return(list(x = 0, y = 0))
    }
  )
  
  lambda <- if (max(bc$y) - min(bc$y) < 1e-6) 0 else bc$x[which.max(bc$y)]
  if (abs(lambda) < 1e-6) {
    scaled <- log(x + 1e-8)
  } else {
    scaled <- (x^lambda - 1) / lambda
  }
  
  return(list(
    scaled_data = scaled,
    normalize_method = "boxcox_transform",
    normalize_info = list(lambda = lambda)
  ))
}

#' Yeo-Johnson Transformation for Numeric Data
#'
#' Applies the Yeo-Johnson transformation to a numeric vector, which can handle both positive and negative values.
#' If the transformation fails, an identity transform is used as fallback.
#'
#' @importFrom car powerTransform yjPower
#' @import stats
#' @param x A numeric vector. Missing values will be ignored during lambda estimation, but retained in the output.
#'
#' @returns A list containing:
#' \describe{
#'   \item{scaled_data}{A numeric vector of transformed data, matching the length of \code{x}.}
#'   \item{normalize_method}{A character string indicating the method used: \code{"yeojohnson_transform"}.}
#'   \item{normalize_info}{A list containing the lambda parameter used in the Yeo-Johnson transformation.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' result <- yeojohnson_transform(x)
yeojohnson_transform <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) < 3) {
    stop("Insufficient values for Yeo-Johnson transform")
  }
  
  pt <- tryCatch(
    powerTransform(x_clean, family = "yjPower"),
    error = function(e) {
      warning("Yeo-Johnson failed, using identity transform instead")
      return(list(lambda = 1))
    }
  )
  
  lambda <- pt$lambda
  scaled <- yjPower(x, lambda)
  
  return(list(
    scaled_data = scaled,
    normalize_method = "yeojohnson_transform",
    normalize_info = list(lambda = lambda)
  ))
}

#' Preprocess Numeric Data for Transformation
#'
#' This function preprocesses a numeric vector by handling missing values and ensuring compatibility with transformations such as log or Box-Cox.
#' Missing values are imputed using the median of the available data. If the specified method requires strictly positive values (e.g., log or Box-Cox),
#' a small offset is added to make all values positive.
#'
#' @param x A numeric vector to be preprocessed. May contain missing or non-positive values.
#' @param method A character string indicating the intended transformation method. Supported values: \code{"log_transform"}, \code{"boxcox_transform"}.
#'
#' @returns A numeric vector with missing values imputed and, if applicable, all values adjusted to be strictly positive.
#'
#' @export
#'
#' @examples
#' x <- c(-1, 0, 1, NA, 5)
#' x_preprocessed <- preprocess_data(x, method = "log_transform")
#' x_preprocessed

preprocess_data <- function(x, method) {
  
  if (any(is.na(x))) {
    na_count <- sum(is.na(x))
    x[is.na(x)] <- median(x, na.rm = TRUE)
    warning(paste(na_count, "NA values replaced with median in preprocessing"))
  }
  
  if (method %in% c("log_transform", "boxcox_transform")) {
    min_val <- min(x, na.rm = TRUE)
    if (min_val <= 0) {
      offset <- abs(min_val) + 1e-8
      x <- x + offset
      warning(paste("Added offset", offset, "to make values positive for", method))
    }
  }
  
  return(x)
}


#' Normalize Numeric Columns in a Data Frame
#'
#' This function performs normalization (or transformation) on numeric columns of a data frame.
#' It supports various normalization methods and can automatically select the most suitable method
#' for each column based on distributional properties (e.g., via Shapiro-Wilk test for normality).
#'
#' @importFrom car powerTransform yjPower
#' @importFrom MASS boxcox
#' @import stats
#' @param data A data frame containing the variables to be normalized. Both numeric and non-numeric columns are allowed.
#' @param normalize_method A character string indicating the normalization method to use. Options include:
#'   \itemize{
#'     \item \code{"auto"} — Automatically select normalization method based on normality test.
#'     \item \code{"log_transform"} — Logarithmic transformation.
#'     \item \code{"min_max_scale"} — Min-max scaling.
#'     \item \code{"z_score_standardize"} — Z-score standardization.
#'     \item \code{"max_abs_scale"} — Maximum absolute scaling.
#'     \item \code{"center_data"} — Centering the data (mean = 0).
#'     \item \code{"scale_data"} — Scaling (variance = 1).
#'     \item \code{"boxcox_transform"} — Box-Cox transformation (requires positive values).
#'     \item \code{"yeojohnson_transform"} — Yeo-Johnson transformation (handles zero and negative values).
#'   }
#' @param group_col Optional. A character string specifying the name of the grouping column.
#'   If not found in the data, will be set to \code{NULL}.
#' @param max_unique_values An integer. Variables with fewer than or equal to this number of unique values
#'   will be treated as categorical and excluded from normalization. Default is 5.
#' @param alpha Significance level used in the Shapiro-Wilk normality test when \code{normalize_method = "auto"}. Default is 0.05.
#'
#' @returns A list containing:
#' \describe{
#'   \item{scaled_data}{A data frame of the same dimensions as input, with numeric columns normalized.}
#'   \item{normalize_method}{The specified method used for normalization.}
#'   \item{method_info}{A list describing the normalization applied to each numeric column.}
#'   \item{alpha_threshold}{The alpha used for normality testing (if applicable).}
#'   \item{group_col}{The name of the grouping column, if specified.}
#'   \item{group_info}{Details about the groups if \code{group_col} is provided.}
#'   \item{timestamp}{Timestamp of when normalization was performed.}
#'   \item{summary}{A string summarizing the normalization status across columns.}
#' }
#'
#' @export
#'
#' @examples
#' data <- data.frame(
#'   group = rep(c("A", "B"), each = 5),
#'   x = c(1:5, 2:6),
#'   y = c(10, 20, NA, 40, 50, 60, 70, 80, 90, 100),
#'   z = rep(1, 10)
#' )
#' result <- normalize_data(data, normalize_method = "auto", group_col = "group")
#' head(result$scaled_data)

normalize_data <- function(data,
                           normalize_method = "auto",
                           group_col = "group",
                           max_unique_values = 5,
                           alpha = 0.05,
                           save_dir = here::here("StatObject","Data"),
                           save_data = TRUE,
                           csv_filename = "normalize_data.csv"
                           
) {
  
  stopifnot(is.data.frame(data))
  if (length(group_col) == 0 || is.null(group_col) || !group_col %in% colnames(data)) {
    group_col <- NULL
  }
  
  valid_methods <- c("auto", "log_transform", "min_max_scale", "z_score_standardize",
                     "max_abs_scale", "center_data", "scale_data",
                     "boxcox_transform", "yeojohnson_transform")
  stopifnot(normalize_method %in% valid_methods)
  
  variable_types <- diagnose_variable_type(data, group_col = group_col,
                                           max_unique_values = max_unique_values)
  num_cols <- variable_types$numeric_vars
  count_cols <- setdiff(names(data), num_cols)
  
  if (length(num_cols) == 0) {
    stop("No numeric columns to normalize.")
  }
  
  cat("Processing", length(num_cols), "numeric columns:", paste(num_cols, collapse = ", "), "\n")
  
  scaled_data <- data.frame(matrix(nrow = nrow(data), ncol = length(num_cols)))
  colnames(scaled_data) <- num_cols
  method_info <- list()
  row_names <- rownames(data)
  
  for (col in num_cols) {
    x <- data[[col]]
    x_clean <- na.omit(x)
    
    if (length(x_clean) < 3) {
      warning(sprintf("Column '%s' has too few non-NA values (%d), skipping normalization.",
                      col, length(x_clean)))
      scaled_data[[col]] <- x
      method_info[[col]] <- list(method = "none", reason = "insufficient data")
      next
    }
    
    if (normalize_method == "auto") {
      shapiro_test <- tryCatch(
        shapiro.test(x_clean),
        error = function(e) list(p.value = NA)
      )
      p_value <- shapiro_test$p.value
      is_normal <- !is.na(p_value) && p_value > alpha
      
      if (is_normal) {
        current_method <- "z_score_standardize"
        reason <- paste("normal distribution (Shapiro-Wilk p =", round(p_value, 3), ")")
      } else if (all(x_clean > 0) && !all(x_clean == 0)) {
        current_method <- "boxcox_transform"
        reason <- paste("positive skewed distribution (Shapiro-Wilk p =", round(p_value, 3), ")")
      } else {
        current_method <- "yeojohnson_transform"
        reason <- paste("non-normal distribution with zero/negative values (Shapiro-Wilk p =",
                        round(p_value, 3), ")")
      }
      cat(sprintf("- Column '%s': %s (%s)\n", col, current_method, reason))
    } else {
      current_method <- normalize_method
    }
    
    x_preprocessed <- tryCatch(
      preprocess_data(x, current_method),
      error = function(e) {
        warning(sprintf("Preprocessing failed for column '%s': %s", col, e$message))
        x
      }
    )
    
    result <- tryCatch(
      {
        switch(
          current_method,
          "log_transform" = log_transform(x_preprocessed),
          "min_max_scale" = min_max_scale(x_preprocessed),
          "z_score_standardize" = z_score_standardize(x_preprocessed),
          "max_abs_scale" = max_abs_scale(x_preprocessed),
          "center_data" = center_data(x_preprocessed),
          "scale_data" = scale_data(x_preprocessed),
          "boxcox_transform" = boxcox_transform(x_preprocessed),
          "yeojohnson_transform" = yeojohnson_transform(x_preprocessed),
          stop("Unknown normalization method")
        )
      },
      error = function(e) {
        warning(sprintf("Primary method %s failed for column '%s', trying z-score: %s",
                        current_method, col, e$message))
        tryCatch(
          z_score_standardize(x_preprocessed),
          error = function(e) {
            warning(sprintf("Fallback method also failed for column '%s': %s", col, e$message))
            list(scaled_data = x,
                 normalize_method = "failed",
                 normalize_info = list(error = e$message))
          }
        )
      }
    )
    
    scaled_col <- if (is.data.frame(result$scaled_data)) {
      result$scaled_data[[1]]
    } else {
      result$scaled_data
    }
    scaled_data[[col]] <- scaled_col
    
    method_info[[col]] <- c(
      list(
        method = result$normalize_method,
        status = ifelse(result$normalize_method == "failed", "failed", "success"),
        original_method = current_method
      ),
      result$normalize_info,
      if (normalize_method == "auto") list(
        shapiro_p = p_value,
        is_normal = is_normal,
        auto_selection_reason = reason
      )
    )
  }
  
  final_data <- data
  final_data[, num_cols] <- scaled_data
  
  rownames(final_data) <- row_names
  
  normalized_result <- list(
    scaled_data = final_data,
    normalize_method = normalize_method,
    method_info = method_info,
    alpha_threshold = if (normalize_method == "auto") alpha else NULL,
    group_col = group_col
  )
  
  if (!is.null(group_col)) {
    normalized_result$group_info <- list(
      group_col = group_col,
      groups = unique(data[[group_col]])
    )
  }
  
  success_count <- sum(sapply(method_info, function(x) x$status == "success"))
  normalized_result$summary <- sprintf(
    "Processed %d columns (%d successful, %d skipped/failed)",
    length(num_cols),
    success_count,
    length(num_cols) - success_count
  )
  nmdat <- normalized_result$scaled_data
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(nmdat, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  
  cat(normalized_result$summary, "\n")
  return(normalized_result)
}


#' Perform Normalization on Statistical Object or Data Frame
#'
#' This function normalizes the numeric variables in the provided object (either a `Stat` object
#' or a `data.frame`) using the specified normalization method. It automatically detects the
#' appropriate normalization strategy if `normalize_method = "auto"`.
#'
#' @importFrom car powerTransform yjPower
#' @importFrom MASS boxcox
#' @import stats
#' @import methods
#' @param object An object of class `Stat` (with slots `clean.data`, `group_col`, `scale.data`,
#'        and `normalization.info`) or a `data.frame`.
#' @param normalize_method The normalization method to apply. One of: `"auto"`,
#'        `"log_transform"`, `"min_max_scale"`, `"z_score_standardize"`, `"max_abs_scale"`,
#'        `"center_data"`, `"scale_data"`, `"boxcox_transform"`, `"yeojohnson_transform"`.
#' @param group_col The column name for grouping (optional, used for diagnostics).
#' @param max_unique_values Maximum unique values for categorical variable detection. Default is 5.
#' @param alpha Significance level for normality test in `"auto"` mode. Default is 0.05.
#'
#' @returns
#' If `object` is a `Stat` object, returns the updated object with `scale.data` and `normalization.info` slots updated.
#' If `object` is a `data.frame`, returns a list containing the normalized data and metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For data.frame input
#' result <- stat_normalize_process(my_data, normalize_method = "z_score_standardize")
#'
#' # For Stat object
#' my_stat <- stat_normalize_process(my_stat_object)
#' }

stat_normalize_process <- function(object,
                                   normalize_method = "auto",
                                   group_col = "group",
                                   max_unique_values = 5,
                                   alpha = 0.05,
                                   save_dir = here::here("StatObject","Data"),
                                   save_data = TRUE,
                                   csv_filename = "normalize_data.csv") {
  cat("Input object class:", class(object), "\n")
  
  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    group_col <- slot(object, "group_col")
    
    if (length(group_col) == 0 || is.null(group_col) || !group_col %in% colnames(data)) {
      cat("Group column is not valid, setting to NULL.\n")
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
  
  cat("Starting normalization process...\n")
  
  nm_result <- normalize_data(data,
                              normalize_method = normalize_method,
                              group_col = group_col,
                              max_unique_values = max_unique_values,
                              alpha = alpha,
                              save_dir =save_dir,
                              save_data = save_data,
                              csv_filename = csv_filename)
  
  nmdat <- nm_result$scaled_data
  
  stopifnot(is.data.frame(nmdat))
  stopifnot(all(names(data) %in% names(nmdat)))
  stopifnot(nrow(nmdat) == nrow(data))
  if (inherits(object, "Stat")) {
    if (!is.null(slotNames(object)) && "scale.data" %in% slotNames(object)) {
      object@scale.data <- nmdat
      object@process.info[["normalization.info"]] <- nm_result
      cat("Updating 'Stat' object...\n")
      cat("The 'Stat' object has been updated with the following slots:\n")
      cat("- 'scale.data' slot updated\n")
      cat("- 'process.info' slot updated\n")
    } else {
      stop("The 'Stat' object does not have required slots.")
    }
    return(object)
  }
  
  cat("Normalization complete, returning normalized data frame.\n")
  return(nm_result)
}
