#' Impute Missing Values in Data
#'
#' @param data A data frame.
#' @param group_col Group column.
#' @param impute_method Method ("mice" or "median_mode").
#' @param m Number of imputations.
#' @param max_unique_values Max unique values for categorical.
#' @param return_imputation_info Logical.
#' @param save_dir Save directory.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
impute_missing_values <- function(data,
                                  group_col = "group",
                                  impute_method = "mice",
                                  m = 5,
                                  max_unique_values = 5,
                                  return_imputation_info = T,
                                  save_dir = here::here("StatObject","Data"),
                                  save_data = TRUE,
                                  csv_filename = "clean_data.csv") {
  if (!is.data.frame(data)) stop("Input must be a data frame.")
  impute_method <- match.arg(impute_method, choices = c("mice", "median_mode"))
  
  imputation_info <- list(
    impute_method = impute_method,
    original_na_counts = colSums(is.na(data)),
    variable_types = NULL,
    imputation_params = list(
      m = if(impute_method == "mice") m else NULL,
      max_unique_values = max_unique_values
    ),
    imputation_values = list(),
    imputation_details = list()
  )
  
  variable_types <- diagnose_variable_type(data, group_col, max_unique_values)
  numeric_vars <- variable_types$numeric_vars
  categorical_vars <- variable_types$categorical_vars
  imputation_info$variable_types <- variable_types
  
  calculate_mode <- function(x) {
    ux <- unique(na.omit(x))
    if(length(ux) == 0) return(NA)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  mic_median_impute <- function(data, m) {
    cat("Performing multiple imputation using MICE...\n")
    
    imputation_info$imputation_params$mice_method <<- "pmm"
    imputation_info$imputation_params$maxit <<- 5
    imputation_info$imputation_params$seed <<- 123
    
    imp_data <- mice::mice(data, m = m, method = 'pmm', maxit = 5, seed = 123)
    imp_data_data <- mice::complete(imp_data, action = 3)
    
    imputation_info$imputation_details$mice_imp <<- imp_data
    
    cat("Imputed data (MICE) generated. Now performing median/mode imputation for remaining NAs...\n")
    
    imp_data_data[] <- lapply(names(imp_data_data), function(col) {
      x <- imp_data_data[[col]]
      if (col %in% numeric_vars) {
        median_val <- median(x, na.rm = TRUE)
        imputation_info$imputation_values[[col]] <<- list(
          method = if(any(is.na(x))) "median (post-MICE)" else "MICE only",
          value = median_val,
          used_value = median_val
        )
        ifelse(is.na(x), median_val, x)
      } else if (col %in% categorical_vars) {
        mode_val <- calculate_mode(x)
        imputation_info$imputation_values[[col]] <<- list(
          method = if(any(is.na(x))) "mode (post-MICE)" else "MICE only",
          value = mode_val,
          used_value = mode_val
        )
        ifelse(is.na(x), mode_val, x)
      } else {
        x
      }
    })
    
    for (col in categorical_vars) {
      imp_data_data[[col]] <- as.factor(imp_data_data[[col]])
    }
    
    return(imp_data_data)
  }
  
  if (impute_method == "median_mode") {
    cat("Using simple median and mode imputation method...\n")
    
    for (col in categorical_vars) {
      mode_val <- calculate_mode(data[[col]])
      imputation_info$imputation_values[[col]] <- list(
        method = "mode",
        value = mode_val,
        used_value = mode_val
      )
      data[[col]][is.na(data[[col]])] <- mode_val
    }
    
    for (col in numeric_vars) {
      median_val <- median(data[[col]], na.rm = TRUE)
      imputation_info$imputation_values[[col]] <- list(
        method = "median",
        value = median_val,
        used_value = median_val
      )
      data[[col]][is.na(data[[col]])] <- median_val
    }
  } else {
    data <- mic_median_impute(data, m)
  }
  
  imputation_info$final_na_counts <- colSums(is.na(data))
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  imputation_method = impute_method
  cat("Imputation completed successfully.\n")
  
  if (return_imputation_info) {
    return(list(
      imputed_data = data,
      imputation_method =imputation_method,
      imputation_info = imputation_info
    ))
  } else {
    return(data)
  }
}


#' Process Missing Values for Stat Object or Data Frame
#'
#' @param object Stat object or data frame.
#' @param m MICE imputations.
#' @param impute_method Method.
#' @param group_col Group column.
#' @param miss_threshold Missing threshold.
#' @param max_unique_values Max unique values.
#' @param return_imputation_info Logical.
#' @param save_dir Save directory.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
stat_miss_processed <- function(object,
                                m = 5,
                                impute_method = "mice",
                                group_col = "group",
                                miss_threshold = 20,
                                max_unique_values = 5,
                                return_imputation_info = TRUE,
                                save_dir = here::here("StatObject","Data"),
                                save_data = TRUE,
                                csv_filename = "clean_data.csv") {
  impute_method <- match.arg(impute_method, choices = c("mice", "median_mode"))
  
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
  
  variable_types <- diagnose_variable_type(data, 
                                           group_col = group_col, 
                                           max_unique_values = max_unique_values)
  if (length(variable_types$numeric_vars) == 0 && 
      length(variable_types$categorical_vars) == 0) {
    stop("No valid variables found after variable type diagnosis")
  }
  
  impute_result <- impute_missing_values(
    data = data,
    group_col = group_col,
    impute_method = impute_method,
    m = m,
    return_imputation_info = return_imputation_info,
    save_dir =save_dir,
    save_data = save_data,
    csv_filename = csv_filename)
  
  
  if (return_imputation_info) {
    clean_miss_data <- impute_result$imputed_data
    missing_info <- impute_result$imputation_info
  } else {
    clean_miss_data <- impute_result
    missing_info <- NULL
  }
  
  if (ncol(clean_miss_data) == 0 || nrow(clean_miss_data) == 0) {
    stop("No data remains after missing value imputation")
  }
  

  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@clean.data <- clean_miss_data
    cat("- 'clean.data' slot updated.\n")
    if (return_imputation_info && !is.null(impute_result)) {
      object@process.info[["missing_info"]] <- impute_result
      cat("- 'process.info' slot updated.\n")
    }
    return(object) 
  } else {
    return(impute_result)}
  
} 
