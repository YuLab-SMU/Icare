#' Process New Data Using Existing Stat Object
#'
#' This function applies the preprocessing steps (missing value handling, outlier handling, normalization)
#' from an existing trained `Stat` object to a new dataset.
#'
#' @param stat_object A trained `Stat` object containing processing information.
#' @param new_data A new data frame to be processed.
#' @param group_col Group column name in the new data (optional).
#' @param max_unique_values Max unique values for variable diagnosis.
#' @param save_dir Directory to save results.
#' @param save_data Logical.
#' @return A processed data frame.
#' @export
#' @examples
#' # Train a Stat object with processing steps
#' stat <- CreateStatObject(clean.data = mtcars, group_col = "cyl")
#' stat <- stat_normalize_process(stat, method = "z_score")
#' 
#' # Create new data (first 3 rows)
#' new_data <- mtcars[1:3, ]
#' processed <- process_new_data(stat, new_data, save_data = FALSE)
#' head(processed)
process_new_data <- function(stat_object,
                             new_data,
                             group_col = "group",
                             max_unique_values = 5,
                             save_dir = NULL,
                             save_data = TRUE) {
  
  if (!inherits(stat_object, "Stat")) stop("stat_object must be of class 'Stat'.")
  if (!is.data.frame(new_data)) stop("new_data must be a data frame.")
  
  process_info <- stat_object@process.info
  if (length(process_info) == 0) warning("No processing info found in stat_object.")
  
  # 1. Variable Type Diagnosis (on new data)
  variable_types <- diagnose_variable_type(new_data, group_col, max_unique_values)
  
  processed_data <- new_data
  
  # 2. Missing Value Handling
  if (!is.null(process_info$missing_removal)) {
    cat("Applying missing value removal...\n")
    miss_threshold <- process_info$missing_removal$miss_threshold
    high_missing_vars <- process_info$missing_removal$high_missing_vars
    
    if (length(high_missing_vars) > 0) {
      processed_data <- processed_data[, !names(processed_data) %in% high_missing_vars, drop = FALSE]
    }
    # Note: We typically don't remove samples from new data based on training threshold, 
    # but we could check sample missingness. For now, we only remove variables dropped in training.
  }
  
  if (!is.null(process_info$missing_info)) {
    cat("Applying missing value imputation...\n")
    impute_info <- process_info$missing_info$imputation_info
    impute_method <- impute_info$impute_method
    imputation_values <- impute_info$imputation_values
    
    for (col in names(imputation_values)) {
      if (col %in% names(processed_data)) {
        fill_val <- imputation_values[[col]]$used_value
        if (any(is.na(processed_data[[col]]))) {
          processed_data[[col]][is.na(processed_data[[col]])] <- fill_val
        }
      }
    }
  }
  
  # 3. Outlier Handling
  if (!is.null(process_info$outlier_handling)) {
    cat("Applying outlier handling to new data using training thresholds...\n")
    method <- process_info$outlier_handling$method
    threshold <- process_info$outlier_handling$threshold
    stats <- process_info$outlier_handling$detection_stats
    
    if (!is.null(stats)) {
      for (col in names(stats)) {
        if (col %in% colnames(processed_data) && is.numeric(processed_data[[col]])) {
          x <- processed_data[[col]]
          if (method == "zscore") {
            z <- (x - stats[[col]]$mean) / stats[[col]]$sd
            outliers <- which(abs(z) > threshold)
          } else if (method == "iqr") {
            lower <- stats[[col]]$Q1 - threshold * stats[[col]]$IQR
            upper <- stats[[col]]$Q3 + threshold * stats[[col]]$IQR
            outliers <- which(x < lower | x > upper)
          } else {
            next
          }
          if (length(outliers) > 0) {
            cat("  New outliers detected in column", col, ":", length(outliers), "values\n")
          }
        }
      }
    } else {
      cat("No outlier statistics found in training object. Skipping outlier handling on new data.\n")
    }
  }
  
  # 4. Normalization
  if (!is.null(process_info$normalization)) {
    cat("Applying normalization...\n")
    norm_methods <- process_info$normalization
    # NEW: per-column fitted parameters (mean/sd/min/max/lambda/shift) learned
    # on the training data by normalize_data(). Older Stat objects saved
    # before this field existed will have this as NULL -- handled below.
    norm_params  <- process_info$normalization_params
    
    for (col in names(norm_methods)) {
      if (col %in% names(processed_data) && is.numeric(processed_data[[col]])) {
        method <- norm_methods[[col]]
        x <- processed_data[[col]]
        params <- if (!is.null(norm_params)) norm_params[[col]] else NULL
        
        if (!is.null(params)) {
          # Correct path: reproduce the exact training-time transform using
          # the stored parameters, so new data ends up on the same scale the
          # model was trained on. This also works for a single new
          # observation (n = 1), where recomputing mean/sd or min/max from
          # scratch would be undefined.
          processed_data[[col]] <- switch(method,
                                          "log"         = .apply_log(x, params),
                                          "min_max"     = .apply_min_max(x, params),
                                          "z_score"     = .apply_z_score(x, params),
                                          "center"      = .apply_center(x, params),
                                          "scale"       = .apply_scale(x, params),
                                          "max_abs"     = .apply_max_abs(x, params),
                                          "box_cox"     = .apply_box_cox(x, params),
                                          "yeo_johnson" = .apply_yeo_johnson(x, params),
                                          x # Default
          )
        } else {
          # Backward-compatible fallback for Stat objects saved before
          # normalization_params was introduced: recompute independently on
          # new_data, as before. This reproduces the original inconsistency,
          # so it is loudly flagged -- re-running stat_normalize_process() on
          # the training Stat object (and re-saving it) removes this warning.
          warning("No stored normalization parameters found for column '", col,
                  "' (older Stat object, saved before this fix). Falling back ",
                  "to recomputing '", method, "' independently on new_data -- ",
                  "results may not be on the same scale as at training time. ",
                  "Re-run stat_normalize_process() on the training Stat object ",
                  "and re-save it to fix this permanently.")
          processed_data[[col]] <- switch(method,
                                          "log" = log_transform(x),
                                          "min_max" = min_max_scale(x),
                                          "z_score" = z_score_standardize(x),
                                          "center" = center_data(x),
                                          "scale" = scale_data(x),
                                          "max_abs" = max_abs_scale(x),
                                          "box_cox" = boxcox_transform(x),
                                          "yeo_johnson" = yeojohnson_transform(x),
                                          x # Default
          )
        }
      }
    }
  }
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    full_path <- file.path(save_dir, "new_data_processed.csv")
    write.csv(processed_data, file = full_path, row.names = FALSE)
    cat("Processed new data saved to:", full_path, "\n")
  }
  
  return(processed_data)
}
