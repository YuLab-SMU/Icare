#' Apply Previous Normalization Parameters to New Data
#'
#' This function applies the normalization parameters learned from a previous dataset (training set)
#' to a new dataset. The normalization method can be chosen automatically based on the previous dataset's
#' settings or manually specified. Supported methods include log transformation, min-max scaling, z-score
#' standardization, and more. The function returns the normalized data and a summary of the operations performed.
#'
#' @importFrom car yjPower
#' @import stats
#' @param new_data A data frame containing the new dataset that needs to be normalized.
#' @param prev_normalization_info A list containing the normalization method and parameters used during the
#' training phase. This should include `method_info` (which stores the method parameters for each variable)
#' and `normalize_method` (which stores the normalization method applied).
#' @param normalize_method A string specifying the normalization method to be used. If "auto", the method
#' will be selected based on the training dataset's normalization parameters. If any other valid method is
#' specified, that method will be applied directly to the new data. Defaults to "auto".
#' @param group_col A string specifying the column name that represents the grouping variable in the new data.
#' Used for diagnosing variable types and handling categorical variables. Defaults to "group".
#' @param alpha A numeric value used for thresholding purposes if the method is "auto". Defaults to 0.05.
#'
#' @returns A list containing the following components:
#' - `scaled_data`: The normalized data frame, with the same structure as `new_data`, but with the numeric
#'   columns transformed according to the selected normalization method.
#' - `normalize_method`: The method applied for normalization.
#' - `method_info`: A list of method details for each column, including method, parameters, and status (success/fail).
#' - `alpha_threshold`: If "auto" was used for normalization, this field contains the alpha threshold applied.
#' - `prev_normalization_info`: The original normalization information passed in for reference.
#' - `timestamp`: The time when the normalization was applied.
#' - `group_info`: Information about the group column (if specified), including the column name and unique groups.
#' - `summary`: A summary string detailing the number of columns processed, successful operations, and failures.
#'
#' @export
#'
#' @examples
#' # Example usage of the function
#' prev_norm_info <- list(
#'   normalize_method = "z_score_standardize",
#'   method_info = list(
#'     column1 = list(method = "z_score_standardize", mean_vals = 10, sd_vals = 5)
#'   )
#' )
#' new_data <- data.frame(
#'   column1 = c(12, 15, 18, 10),
#'   column2 = c(20, 21, 22, 19)
#' )
#' result <- apply_previous_normalization(new_data, prev_norm_info)
#'
#' # Check the normalized data and summary
#' result$scaled_data
#' result$summary
apply_previous_normalization <- function(new_data,
                                         prev_normalization_info,
                                         normalize_method = "auto",
                                         group_col = "group",
                                         alpha = 0.05) {

  stopifnot(is.data.frame(new_data))
  if (!is.list(prev_normalization_info) ||
      is.null(prev_normalization_info$method_info) ||
      is.null(prev_normalization_info$normalize_method)) {
    stop("Invalid previous normalization info. Must contain method_info and normalize_method.")
  }


  valid_methods <- c("auto", "log_transform", "min_max_scale", "z_score_standardize",
                     "max_abs_scale", "center_data", "scale_data",
                     "boxcox_transform", "yeojohnson_transform")
  stopifnot(normalize_method %in% valid_methods)

  variable_types <- diagnose_variable_type(new_data, group_col = group_col, max_unique_values = 5)
  num_cols <- variable_types$numeric_vars
  count_cols <- setdiff(names(new_data), num_cols)

  if (length(num_cols) == 0) {
    stop("No numeric columns found in the new data.")
  }

  scaled_data <- new_data[, num_cols, drop = FALSE]
  method_info <- list()

  train_method_info <- prev_normalization_info$method_info
  train_method <- prev_normalization_info$normalize_method

  for (col in num_cols) {
    x <- new_data[[col]]

    if (normalize_method == "auto") {
      if (col %in% names(train_method_info)) {
        current_method <- train_method_info[[col]]$method
        reason <- if (!is.null(train_method_info[[col]]$auto_selection_reason)) {
          train_method_info[[col]]$auto_selection_reason
        } else {
          paste("using same method as training:", current_method)
        }
      } else {
        current_method <- "z_score_standardize"
        reason <- "new column, defaulting to z-score"
      }
      cat(sprintf("- Column '%s': %s (%s)\n", col, current_method, reason))
    } else {
      current_method <- normalize_method
      reason <- paste("manually specified method:", current_method)
    }

    params <- if (col %in% names(train_method_info)) {
      train_method_info[[col]]
    } else {
      warning(sprintf("No normalization parameters found for column '%s' in training data", col))
      list()
    }

    result <- tryCatch(
      {
        switch(
          current_method,
          "log_transform" = {
            offset <- if (!is.null(params$offset)) params$offset else 1e-8
            list(
              scaled_data = log(x + offset),
              normalize_method = "log_transform",
              normalize_info = list(offset = offset)
            )
          },
          "min_max_scale" = {
            if (!col %in% names(train_method_info) ||
                is.null(train_method_info[[col]]$min_vals) ||
                is.null(train_method_info[[col]]$max_vals)) {
              stop("Missing min/max values from training data")
            }
            min_val <- train_method_info[[col]]$min_vals
            max_val <- train_method_info[[col]]$max_vals
            if (max_val == min_val) max_val <- min_val + 1e-8
            list(
              scaled_data = (x - min_val) / (max_val - min_val),
              normalize_method = "min_max_scale",
              normalize_info = list(min_vals = min_val, max_vals = max_val)
            )
          },
          "z_score_standardize" = {
            if (!col %in% names(train_method_info) ||
                is.null(train_method_info[[col]]$mean_vals) ||
                is.null(train_method_info[[col]]$sd_vals)) {
              stop("Missing mean/sd values from training data")
            }
            mean_val <- train_method_info[[col]]$mean_vals
            sd_val <- train_method_info[[col]]$sd_vals
            if (sd_val == 0) sd_val <- 1
            list(
              scaled_data = (x - mean_val) / sd_val,
              normalize_method = "z_score_standardize",
              normalize_info = list(mean_vals = mean_val, sd_vals = sd_val)
            )
          },
          "max_abs_scale" = {
            if (!col %in% names(train_method_info) || is.null(train_method_info[[col]]$max_abs)) {
              stop("Missing max_abs values from training data")
            }
            max_abs <- train_method_info[[col]]$max_abs
            if (max_abs == 0) max_abs <- 1
            list(
              scaled_data = x / max_abs,
              normalize_method = "max_abs_scale",
              normalize_info = list(max_abs = max_abs)
            )
          },
          "center_data" = {
            if (!col %in% names(train_method_info) || is.null(train_method_info[[col]]$center)) {
              stop("Missing center values from training data")
            }
            center <- train_method_info[[col]]$center
            list(
              scaled_data = x - center,
              normalize_method = "center_data",
              normalize_info = list(center = center)
            )
          },
          "scale_data" = {
            if (!col %in% names(train_method_info) || is.null(train_method_info[[col]]$scale)) {
              stop("Missing scale values from training data")
            }
            scale_val <- train_method_info[[col]]$scale
            if (scale_val == 0) scale_val <- 1
            list(
              scaled_data = x / scale_val,
              normalize_method = "scale_data",
              normalize_info = list(scale = scale_val)
            )
          },
          "boxcox_transform" = {
            if (!col %in% names(train_method_info) || is.null(train_method_info[[col]]$lambda)) {
              stop("Missing lambda value from training data")
            }
            lambda <- train_method_info[[col]]$lambda
            offset <- if (!is.null(train_method_info[[col]]$offset)) train_method_info[[col]]$offset else 0
            x_adj <- x + offset
            if (abs(lambda) < 1e-6) {
              scaled <- log(x_adj)
            } else {
              scaled <- (x_adj^lambda - 1) / lambda
            }
            list(
              scaled_data = scaled,
              normalize_method = "boxcox_transform",
              normalize_info = list(lambda = lambda, offset = offset)
            )
          },
          "yeojohnson_transform" = {
            if (!col %in% names(train_method_info) || is.null(train_method_info[[col]]$lambda)) {
              stop("Missing lambda value from training data")
            }
            lambda <- train_method_info[[col]]$lambda
            scaled <- yjPower(x, lambda)
            list(
              scaled_data = scaled,
              normalize_method = "yeojohnson_transform",
              normalize_info = list(lambda = lambda)
            )
          },
          stop("Unknown normalization method")
        )
      },
      error = function(e) {
        warning(sprintf("Failed to normalize column '%s' with %s: %s",
                        col, current_method, e$message))
        list(scaled_data = x,
             normalize_method = "failed",
             normalize_info = list(error = e$message))
      }
    )

    scaled_data[[col]] <- result$scaled_data
    method_info[[col]] <- c(
      list(
        method = result$normalize_method,
        status = ifelse(result$normalize_method == "failed", "failed", "success"),
        original_method = current_method
      ),
      result$normalize_info,
      list(reason = reason)
    )
  }

  if (length(count_cols) > 0) {
    scaled_data <- cbind(scaled_data, new_data[, count_cols, drop = FALSE])
  }

  normalized_result <- list(
    scaled_data = scaled_data,
    normalize_method = normalize_method,
    method_info = method_info,
    alpha_threshold = if (normalize_method == "auto") alpha else NULL,
    prev_normalization_info = prev_normalization_info,
    timestamp = Sys.time()
  )

  if (!is.null(group_col) && group_col %in% colnames(new_data)) {
    normalized_result$group_info <- list(
      group_col = group_col,
      groups = unique(new_data[[group_col]])
    )
  }

  success_count <- sum(sapply(method_info, function(x) x$status == "success"))
  normalized_result$summary <- sprintf(
    "Processed %d columns (%d successful, %d skipped/failed) using training set parameters",
    length(num_cols),
    success_count,
    length(num_cols) - success_count
  )

  cat(normalized_result$summary, "\n")
  return(normalized_result)
}






#' Normalize Dataset Columns
#'
#' This function normalizes datasets (training, testing, validation) using specified methods
#' and optionally saves the normalized data to CSV files.
#'
#' @param object Either a Train_Model object or a list containing datasets
#' @param normalize_method Normalization method ("auto", "scale", "range", etc.)
#' @param group_col Grouping column name (default: "group")
#' @param max_unique_values Maximum unique values for categorical detection (default: 5)
#' @param alpha Significance level for normality tests (default: 0.05)
#' @param save_data Whether to save normalized data as CSV files (default: TRUE)
#' @param train_filename Filename for normalized training data (default: "normalized_train.csv")
#' @param test_filename Filename for normalized test data (default: "normalized_test.csv")
#' @param val_filename Filename for normalized validation data (default: "normalized_val.csv")
#' @param ext_val_filename Filename for external validation data (default: "normalized_ext_val.csv")
#' @param save_dir Directory to save CSV files (default: here::here("ModelData","Data"))
#' @return Normalized data in object or list format
#' @export
NormalizeData <- function(object,
                          normalize_method = "auto",
                          group_col = "group",
                          max_unique_values = 5,
                          alpha = 0.05,
                          save_data = TRUE,
                          train_filename = "normalized_train.csv",
                          test_filename = "normalized_test.csv",
                          val_filename = "normalized_val.csv",
                          ext_val_filename = "normalized_ext_val.csv",
                          save_dir = here::here("ModelData","Data")) {
  
  cat("\n=== Starting Data Normalization Process ===\n")
  cat("Normalization method:", normalize_method, "\n")
  
  if (inherits(object, 'Train_Model')) {
    group_col <- object@group_col %||% group_col
    training <- slot(object, "split.data")[["training"]]
    testing <- slot(object, "split.data")[["testing"]]
    validation <- if ("validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["validation"]] else NULL
    external_validation <- if ("external_validation" %in% names(slot(object, "split.data"))) 
      slot(object, "split.data")[["external_validation"]] else NULL
    
    if (length(group_col) == 0 || is.null(group_col) || !group_col %in% colnames(training)) {
      cat("-> Group column is not valid, setting to NULL.\n")
      group_col <- NULL
    }
  } else if (is.list(object)) {
    training <- object$training
    testing <- object$testing
    validation <- if ("validation" %in% names(object)) object$validation else NULL
    external_validation <- if ("external_validation" %in% names(object)) object$external_validation else NULL
  } else {
    stop("Input must be an object of class 'Train_Model' or a list")
  }
  
  if (is.null(training) || nrow(training) == 0) {
    stop("No valid training data found in the input")
  }
  if (is.null(testing) || nrow(testing) == 0) {
    stop("No valid test data found in the input")
  }
  
  cat("\n[Step 1] Normalizing TRAINING data...\n")
  normalize_info <- normalize_data(
    training,
    normalize_method = normalize_method,
    group_col = group_col,
    alpha = alpha
  )
  training_norm <- normalize_info$scaled_data
  cat("-> Training data normalized. Variables processed:", ncol(training_norm), "\n")
  
  cat("\n[Step 2] Applying normalization to TEST data...\n")
  testing_norm <- apply_previous_normalization(
    new_data = testing,
    prev_normalization_info = normalize_info,
    normalize_method = normalize_method,
    group_col = group_col,
    alpha = alpha
  )
  cat("-> Test data normalized using training parameters.\n")
  
  validation_norm <- NULL
  if (!is.null(validation)) {
    cat("\n[Step 3] Applying normalization to VALIDATION data...\n")
    validation_norm <- apply_previous_normalization(
      new_data = validation,
      prev_normalization_info = normalize_info,
      normalize_method = normalize_method,
      group_col = group_col,
      alpha = alpha
    )
    cat("-> Validation data normalized using training parameters.\n")
  }
  
  external_norm <- NULL
  if (!is.null(external_validation)) {
    cat("\n[Step 4] Applying normalization to EXTERNAL VALIDATION data...\n")
    external_norm <- apply_previous_normalization(
      new_data = external_validation,
      prev_normalization_info = normalize_info,
      normalize_method = normalize_method,
      group_col = group_col,
      alpha = alpha
    )
    cat("-> External validation data normalized using training parameters.\n")
  }
  
  # Save normalized data to CSV if requested
  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("Created directory:", save_dir, "\n")
      }
      
      write.csv(training_norm, file.path(save_dir, train_filename), row.names = FALSE)
      write.csv(testing_norm$scaled_data, file.path(save_dir, test_filename), row.names = FALSE)
      
      if (!is.null(validation_norm)) {
        write.csv(validation_norm$scaled_data, file.path(save_dir, val_filename), row.names = FALSE)
      }
      
      if (!is.null(external_norm)) {
        write.csv(external_norm$scaled_data, file.path(save_dir, ext_val_filename), row.names = FALSE)
      }
      
      cat("\nSaved normalized datasets to:\n")
      cat("-", file.path(save_dir, train_filename), "\n")
      cat("-", file.path(save_dir, test_filename), "\n")
      if (!is.null(validation_norm)) cat("-", file.path(save_dir, val_filename), "\n")
      if (!is.null(external_norm)) cat("-", file.path(save_dir, ext_val_filename), "\n")
      
    }, error = function(e) {
      warning("Failed to save CSV files: ", e$message)
    })
  }
  
  if (inherits(object, "Train_Model")) {
    object@split.scale.data <- list(
      training = training_norm,
      testing = testing_norm$scaled_data,
      validation = if (!is.null(validation_norm)) validation_norm$scaled_data else NULL,
      external_validation = if (!is.null(external_norm)) external_norm$scaled_data else NULL,
      normalize_method = normalize_method,
      normalization_info = normalize_info,
      testing_normalization_info = testing_norm,
      validation_normalization_info = validation_norm,
      external_normalization_info = external_norm
    )
    
    cat("\nUpdated 'Train_Model' object contains:\n")
    cat("- Normalized datasets in split.scale.data slot:\n")
    cat("  * Training\n  * Testing\n")
    if (!is.null(validation_norm)) cat("  * Validation\n")
    if (!is.null(external_norm)) cat("  * External validation\n")
    return(object)
  } else {
    return(list(
      training = training_norm,
      testing = testing_norm$scaled_data,
      validation = if (!is.null(validation_norm)) validation_norm$scaled_data else NULL,
      external_validation = if (!is.null(external_norm)) external_norm$scaled_data else NULL,
      training_normalization_info = normalize_info,
      testing_normalization_info = testing_norm,
      validation_normalization_info = validation_norm,
      external_normalization_info = external_norm
    ))
  }
}