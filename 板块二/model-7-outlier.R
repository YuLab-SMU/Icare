#' Detect and Mark Outliers in Model Data
#'
#' This function performs outlier detection across training, testing, validation, and external validation datasets.
#' It first establishes normal ranges from the training data, then applies these ranges to identify outliers in other datasets.
#' Results are stored in the Train_Model object's process.info slot.
#'
#' @param object A Train_Model object containing the data to analyze
#' @param save_outliers Logical indicating whether to save outlier information (default: TRUE)
#' @param group_col Column name containing group/class information (default: "group")
#' @param palette_name Color palette name for visualization (default: "Royal1")
#' @param save_plots Logical indicating whether to save outlier plots (default: TRUE)
#' @param save_dir Directory path to save output files (default: here("ModelData", "pre_outlier_info"))
#' @param plot_display_num Number of plots to display interactively (default: 1)
#' @param sub_var Optional subset of variables to analyze (default: NULL)
#' @param plot_width Width of saved plots in inches (default: 5)
#' @param plot_height Height of saved plots in inches (default: 5)
#' @param base_size Base font size for plots (default: 14)
#' @param custom_ranges Optional list of custom ranges for variables (default: NULL)
#' @param max_unique_values Maximum unique values for categorical variables (default: 5)
#'
#' @return The updated Train_Model object with outlier information stored in process.info slot
#'
#' @export
#'
#' @import here 
#' @import methods 
#' @import stats 
#'
#' @examples
#' 
#' # Customized analysis
#' object_model <- ModelDetectOutliers(
#'   object_model,
#'   group_col = "group",
#'   palette_name = "Royal1",
#'   plot_display_num = 3,
#'   custom_ranges = list("AGE" = c(0, 100)
#' )
ModelDetectOutliers <- function(object,
                                save_outliers = TRUE,
                                group_col = "group",
                                palette_name = "Royal1",
                                save_plots = TRUE,
                                save_dir = here("ModelData", "pre_outlier_info"),
                                plot_display_num = 1,
                                sub_var = NULL,
                                plot_width = 5,
                                plot_height = 5,
                                base_size = 14,
                                custom_ranges = NULL,
                                max_unique_values = 5) {
  
  if (!inherits(object, "Train_Model")) {
    stop("Input must be a Train_Model object")
  }
  
  train_data <- slot(object, "split.data")[["training"]]
  test_data <- slot(object, "split.data")[["testing"]]
  validation <- if ("validation" %in% names(slot(object, "split.data"))) 
    slot(object, "split.data")[["validation"]] else NULL
  external_validation <- if ("external_validation" %in% names(slot(object, "split.data"))) 
    slot(object, "split.data")[["external_validation"]] else NULL
  
  group_col <- object@group_col %||% group_col
  
  if (is.null(train_data) || !is.data.frame(train_data) || nrow(train_data) == 0) {
    stop("No valid training data found in the Train_Model object")
  }
  if (is.null(test_data) || !is.data.frame(test_data) || nrow(test_data) == 0) {
    stop("No valid test data found in the Train_Model object")
  }
  
  cat("=== Starting outlier detection process ===\n")
  
  cat("\n[Step 1] Detecting outliers in TRAINING set...\n")
  save_dir1 <- file.path(save_dir, "train_outlier")
  outlier_result_train <- detect_and_mark_outliers(
    data = train_data,
    group_col = group_col,
    palette_name = palette_name,
    save_plots = save_plots,
    save_dir = save_dir1,
    plot_display_num = plot_display_num,
    sub_var = sub_var,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    custom_ranges = custom_ranges,
    max_unique_values = max_unique_values
  )
  
  normal_ranges_list <- if (!is.null(outlier_result_train$normal_ranges)) {
    with(
      outlier_result_train[["normal_ranges"]],
      setNames(
        lapply(seq_len(nrow(outlier_result_train[["normal_ranges"]])), function(i) {
          c(lower_bound[i], upper_bound[i])
        }),
        variable
      )
    )
  } else {
    NULL
  }
  
  cat("\n[Step 2] Applying TRAINING set's outlier ranges to TEST set...\n")
  save_dir2 <- file.path(save_dir, "test_outlier")
  outlier_result_test <- detect_and_mark_outliers(
    data = test_data,
    group_col = group_col,
    palette_name = palette_name,
    save_plots = save_plots,
    save_dir = save_dir2,
    plot_display_num = plot_display_num,
    sub_var = sub_var,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    custom_ranges = normal_ranges_list,
    max_unique_values = max_unique_values
  )
  
  outlier_result_validation <- NULL
  if (!is.null(validation)) {
    cat("\n[Step 3] Applying TRAINING set's outlier ranges to VALIDATION set...\n")
    save_dir3 <- file.path(save_dir, "validation_outlier")
    outlier_result_validation <- detect_and_mark_outliers(
      data = validation,
      group_col = group_col,
      palette_name = palette_name,
      save_plots = save_plots,
      save_dir = save_dir3,
      plot_display_num = plot_display_num,
      sub_var = sub_var,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size,
      custom_ranges = normal_ranges_list,
      max_unique_values = max_unique_values
    )
  }
  
  outlier_result_external <- NULL
  if (!is.null(external_validation)) {
    cat("\n[Step 4] Applying TRAINING set's outlier ranges to EXTERNAL VALIDATION set...\n")
    save_dir4 <- file.path(save_dir, "external_validation_outlier")
    outlier_result_external <- detect_and_mark_outliers(
      data = external_validation,
      group_col = group_col,
      palette_name = palette_name,
      save_plots = save_plots,
      save_dir = save_dir4,
      plot_display_num = plot_display_num,
      sub_var = sub_var,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size,
      custom_ranges = normal_ranges_list,
      max_unique_values = max_unique_values
    )
  }
  
  if (inherits(object, "Train_Model")) {
    object@process.info[["outlier_result"]] <- list(
      outlier_result_train = outlier_result_train,
      outlier_result_test = outlier_result_test
    )
    
    if (!is.null(outlier_result_validation)) {
      object@process.info[["outlier_result"]]$outlier_result_validation <- outlier_result_validation
    }
    if (!is.null(outlier_result_external)) {
      object@process.info[["outlier_result"]]$outlier_result_external <- outlier_result_external
    }
    
    cat("\nThe 'Train_Model' object has been updated with:\n")
    cat("- 'process.info' slot has been updated with outlier annotation information.\n")
    
  }
  
  return(object)
}


#' Apply Specified Outlier Handling Method to Marked Data
#'
#' Processes data that has been marked with outliers (columns ending with "_outlier") 
#' according to the specified handling method. Supports removal, replacement, capping, 
#' or keeping outliers.
#'
#' @param data_marked A data frame containing the original data plus outlier indicator 
#'        columns (ending with "_outlier").
#' @param handling_result A named list containing:
#'        \describe{
#'          \item{method}{The outlier handling method ("remove", "replace", "keep", or "capping")}
#'          \item{replaced_values}{For "replace" method: named list of replacement values}
#'          \item{capping_bounds}{For "capping" method: named list with lower/upper bounds}
#'        }
#'
#' @return A list containing:
#'         \describe{
#'           \item{cleaned_data}{Processed data frame with outliers handled}
#'           \item{method}{The method actually used}
#'           \item{replaced_values}{For "replace": the replacement values used (if applicable)}
#'           \item{capping_bounds}{For "capping": the bounds applied (if applicable)}
#'           \item{parameters}{Additional parameters from handling_result (if any)}
#'         }
#'
#' @export
#'
#' @examples
#' # Example with replacement method
#' marked_data <- data.frame(
#'   value = c(1, 5, 100),
#'   value_outlier = c(FALSE, FALSE, TRUE)
#' )
#' 
#' handling <- list(
#'   method = "replace",
#'   replaced_values = list(value = median(marked_data$value))
#' )
#' 
#' apply_outlier_handling(marked_data, handling)
#'
#' # Example with capping method
#' handling <- list(
#'   method = "capping",
#'   capping_bounds = list(
#'     value = list(lower_bound = 1, upper_bound = 10)
#'   )
#' )
#' 
#' apply_outlier_handling(marked_data, handling)
apply_outlier_handling <- function(data_marked, 
                                   handling_result) {
  
  
  if (!is.list(handling_result) || is.null(handling_result$method)) {
    stop("handling_result must be a list with a 'method' element")
  }
  
  method <- handling_result$method
  
  if (method == "remove") {
    cleaned_data <- handle_outliers_remove(data_marked)
    return(list(cleaned_data = cleaned_data, method = method))
    
  } else if (method == "replace") {
    if (is.null(handling_result$replaced_values)) {
      stop("For 'replace' method, handling_result must contain 'replaced_values'")
    }
    
    replaced_values <- handling_result$replaced_values
    outlier_cols <- grep("_outlier$", colnames(data_marked), value = TRUE)
    
    for (var in outlier_cols) {
      original_var <- sub("_outlier$", "", var)
      
      if (!original_var %in% names(replaced_values)) {
        warning(paste("No replacement value found for variable:", original_var))
        next
      }
      
      data_marked[[original_var]][data_marked[[var]]] <- replaced_values[[original_var]]
    }
    
    cleaned_data <- data_marked %>% dplyr::select(-ends_with("_outlier"))
    
    return(list(
      cleaned_data = cleaned_data,
      method = method,
      replaced_values = replaced_values
    ))
    
  } else if (method == "keep") {
    return(list(cleaned_data = data_marked, method = method))
    
  } else if (method == "capping") {
    if (is.null(handling_result$capping_bounds)) {
      stop("For 'capping' method, handling_result must contain 'capping_bounds'")
    }
    
    capping_bounds <- handling_result$capping_bounds
    outlier_cols <- grep("_outlier$", colnames(data_marked), value = TRUE)
    
    for (var in outlier_cols) {
      original_var <- sub("_outlier$", "", var)
      
      if (!original_var %in% names(capping_bounds)) {
        warning(paste("No capping bounds found for variable:", original_var))
        next
      }
      
      bounds <- capping_bounds[[original_var]]
      lower_bound <- bounds$lower_bound
      upper_bound <- bounds$upper_bound
      
      data_marked[[original_var]] <- pmax(
        pmin(data_marked[[original_var]], upper_bound), 
        lower_bound
      )
    }
    
    cleaned_data <- data_marked %>% dplyr::select(-ends_with("_outlier"))
    
    return(list(
      cleaned_data = cleaned_data,
      method = method,
      capping_bounds = capping_bounds,
      parameters = handling_result$parameters
    ))
    
  } else {
    stop("Unknown method in handling_result")
  }
}


#' Handle Outliers in Model Data
#'
#' Applies specified outlier handling method to training, test, validation, and external validation datasets
#' in a Train_Model object. The handling is first applied to training data, then the same parameters are
#' automatically applied to other datasets to ensure consistency.
#'
#' @param object A Train_Model object containing the datasets to process and outlier detection results
#' @param save_cleaned Logical indicating whether to save cleaned datasets (default: TRUE)
#' @param handle_method Outlier handling method, one of: "replace", "remove", "keep", "capping" (default: "replace")
#' @param lower_quantile Lower quantile threshold for outlier detection (default: 0.05)
#' @param upper_quantile Upper quantile threshold for outlier detection (default: 0.95)
#' @param group_col Column name containing group information (default: "group")
#'
#' @return Returns the modified Train_Model object with:
#'          - Cleaned datasets in the split.data slot
#'          - Outlier handling details stored in process.info slot
#'
#' @export
#'
#' @examples
#' # After detecting outliers with ModelDetectOutliers()
#' object_model <- ModelHandleOutliers(
#'   object = object_model,
#'   handle_method = "replace"
#' )
#'
#'
ModelHandleOutliers <- function(object,
                                save_cleaned = TRUE,
                                handle_method = c("replace", "remove", "keep", "capping"),
                                lower_quantile = 0.05,
                                upper_quantile = 0.95,
                                group_col = "group",
                                save_data = TRUE,
                                train_filename = "train_data.csv",
                                test_filename = "test_data.csv",
                                val_filename = "val_data.csv",
                                ext_val_filename = "ext_val_data.csv",
                                save_dir = here::here("ModelData","Data")) {
  
  handle_method <- match.arg(handle_method)
  
  cat("\n=== Starting outlier handling process ===\n")
  cat("Selected handling method:", handle_method, "\n")
  
  if (!inherits(object, "Train_Model")) {
    stop("Input must be a Train_Model object")
  }
  
  train_data <- slot(object, "split.data")[["training"]]
  test_data <- slot(object, "split.data")[["testing"]]
  validation <- if ("validation" %in% names(slot(object, "split.data"))) 
    slot(object, "split.data")[["validation"]] else NULL
  external_validation <- if ("external_validation" %in% names(slot(object, "split.data"))) 
    slot(object, "split.data")[["external_validation"]] else NULL
  
  group_col <- object@group_col %||% group_col
  outlier_result <- object@process.info[["outlier_result"]]
  
  if (is.null(train_data) || !is.data.frame(train_data) || nrow(train_data) == 0) {
    stop("No valid training data found in the Train_Model object")
  }
  if (is.null(test_data) || !is.data.frame(test_data) || nrow(test_data) == 0) {
    stop("No valid test data found in the Train_Model object")
  }
  if (is.null(outlier_result)) {
    stop("No outlier detection results found. Run ModelDetectOutliers() first.")
  }
  
  cat("\n[Step 1] Handling outliers in TRAINING data...\n")
  train_marked <- outlier_result[["outlier_result_train"]][["data_marked"]]
  handle_outliers_train <- handle_outliers(
    data = train_marked,
    handle_method = handle_method,
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile
  )
  cleaned_train <- handle_outliers_train$cleaned_data
  
  cat("\n[Step 2] Applying same treatment to TEST data...\n")
  test_marked <- outlier_result[["outlier_result_test"]][["data_marked"]]
  handle_outliers_test <- apply_outlier_handling(
    data_marked = test_marked,
    handling_result = handle_outliers_train
  )
  cleaned_test <- handle_outliers_test$cleaned_data
  
  handle_outliers_validation <- NULL
  cleaned_validation <- NULL
  validation_marked <- outlier_result[["outlier_result_validation"]][["data_marked"]]
  
  if (!is.null(validation_marked)) {
    cat("\n[Step 3] Applying same treatment to VALIDATION data...\n")
    handle_outliers_validation <- apply_outlier_handling(
      data_marked = validation_marked,
      handling_result = handle_outliers_train
    )
    cleaned_validation <- handle_outliers_validation$cleaned_data
  }
  
  handle_outliers_external <- NULL
  cleaned_external <- NULL
  external_marked <- outlier_result[["outlier_result_external"]][["data_marked"]]
  
  if (!is.null(external_validation)) {
    cat("\n[Step 4] Applying same treatment to EXTERNAL VALIDATION data...\n")
    handle_outliers_external <- apply_outlier_handling(
      data_marked = external_marked,
      handling_result = handle_outliers_train
    )
    cleaned_external <- handle_outliers_external$cleaned_data
  }
  
  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("\nCreated directory:", save_dir, "\n")
      }
      
      save_cleaned_data <- function(data, filename) {
        if (!is.null(data)) {
          full_path <- file.path(save_dir, filename)
          write.csv(data, file = full_path, row.names = FALSE)
          cat("- Saved", filename, "to", save_dir, "\n")
        }
      }
      
      cat("\n[Step 5] Saving cleaned datasets to CSV...\n")
      save_cleaned_data(cleaned_train, train_filename)
      save_cleaned_data(cleaned_test, test_filename)
      save_cleaned_data(cleaned_validation, val_filename)
      save_cleaned_data(cleaned_external, ext_val_filename)
      
    }, error = function(e) {
      warning("Failed to save CSV files: ", e$message)
    })
  }
  
  object@split.data[["training"]] <- cleaned_train
  object@split.data[["testing"]] <- cleaned_test
  
  update_list <- list(
    handle_outliers_train = handle_outliers_train,
    handle_outliers_test = handle_outliers_test,
    handling_method = handle_method,
    csv_files = if (save_data) list(
      training = file.path(save_dir, train_filename),
      testing = file.path(save_dir, test_filename),
      validation = if (!is.null(cleaned_validation)) file.path(save_dir, val_filename) else NULL,
      external_validation = if (!is.null(cleaned_external)) file.path(save_dir, ext_val_filename) else NULL
    ) else NULL
  )
  
  if (!is.null(handle_outliers_validation)) {
    update_list$handle_outliers_validation <- handle_outliers_validation
    object@split.data[["validation"]] <- cleaned_validation
  }
  if (!is.null(handle_outliers_external)) {
    update_list$handle_outliers_external <- handle_outliers_external
    object@split.data[["external_validation"]] <- cleaned_external
  }
  
  object@process.info[["outlier_result"]] <- modifyList(
    object@process.info[["outlier_result"]] %||% list(),
    update_list
  )
  
  cat("\nUpdated Train_Model object contains:\n")
  cat("- Cleaned versions of all datasets in split.data slot\n")
  cat("- Outlier handling results for:\n")
  cat("  * Training data\n  * Test data\n")
  if (!is.null(handle_outliers_validation)) cat("  * Validation data\n")
  if (!is.null(handle_outliers_external)) cat("  * External validation data\n")
  
  return(object)
}