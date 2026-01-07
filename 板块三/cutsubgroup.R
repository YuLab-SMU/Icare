
find_roc_cutoff <- function(data,
                            variables,
                            time_col = "time",
                            status_col = "status",
                            seed = 123) {
  set.seed(seed)
  
  if (!(status_col %in% names(data)) | !(time_col %in% names(data))) {
    stop("Error: time or status column does not exist in the data.")
  }
  
  
  new_data <- data
  
  
  for (var in variables) {
    cat("Processing variable:", var, "\n")
    
    if (is.numeric(data[[var]]) && var != time_col) {
      valid_data <- data[!is.na(data[[status_col]]) & !is.na(data[[var]]), ]
      
      cat("Valid data points for", var, ":", nrow(valid_data), "\n")
      
      if (nrow(valid_data) < 2) {
        warning("Variable", var, "has insufficient valid data points for ROC analysis. Skipping.")
        next
      }
      
      roc_result <- tryCatch(
        roc(valid_data[[status_col]], valid_data[[var]]),
        warning = function(w) {
          warning("Warning in ROC analysis: ", w$message)
          NULL
        },
        error = function(e) {
          warning("Error in ROC analysis: ", e$message)
          NULL
        }
      )
      
      if (is.null(roc_result)) {
        next
      }
      
      cut_off <- tryCatch(
        pROC::coords(roc_result, "best", ret = "threshold"),
        error = function(e) {
          warning("Error in determining threshold: ", e$message)
          NULL
        }
      )
      
      if (is.null(cut_off)) {
        next
      }
      
      if (is.list(cut_off) || length(cut_off) > 1) {
        cut_off <- cut_off[[1]]
      }
      
      cat("Optimal cutoff for", var, ":", cut_off, "\n")
      
      new_column <- ifelse(data[[var]] > cut_off, paste0(">", round(cut_off, 2)), paste0("<=", round(cut_off, 2)))
      unique_levels <- unique(new_column)
      
      new_data[[var]] <- factor(new_column, levels = unique_levels)
      
      cat("Finished processing variable:", var, "\n\n")
    } else {
      cat("Skipping non-numeric variable:", var, "\n")
    }
  }
  
  return(new_data)
}

library(maxstat)
library(boot)
library(survival)

find_maxstat_cutoff<- function(data,
                               variables,
                               time_col = "time",
                               status_col = "status",
                               n_iterations = 10000,
                               min_group_size = 0.3,
                               seed = 123) {
  
  set.seed(seed)
  
  if (!(time_col %in% colnames(data)) | !(status_col %in% colnames(data))) {
    stop("Error: time or status column does not exist in the data.")
  }
  
  time <- data[[time_col]]
  status <- data[[status_col]]
  
  new_data <- data
  
  for (var in variables) {
    cat("Processing variable:", var, "\n")
    
    if (is.numeric(data[[var]]) && var != time_col) {
      valid_data <- data[!is.na(data[[status_col]]) & !is.na(data[[var]]), ]
      
      cat("Valid data points for", var, ":", nrow(valid_data), "\n")
      
      if (nrow(valid_data) < 2) {
        warning("Variable", var, "has insufficient valid data points for MaxStat analysis. Skipping.")
        next
      }
      
      optimal_cutoffs <- numeric(n_iterations)
      
      for (i in 1:n_iterations) {
        sample_indices <- sample(1:nrow(valid_data), size = floor(0.7 * nrow(valid_data)), replace = TRUE)
        bootstrap_data <- valid_data[sample_indices, ]
        
        maxstat_result <- tryCatch(
          maxstat.test(Surv(bootstrap_data[[time_col]], bootstrap_data[[status_col]]) ~ bootstrap_data[[var]],
                       data = bootstrap_data),
          error = function(e) {
            warning("Error in MaxStat analysis: ", e$message)
            NULL
          }
        )
        
        if (!is.null(maxstat_result)) {
          optimal_cutoffs[i] <- maxstat_result$statistic
        }
      }
      
      cutoff_table <- table(optimal_cutoffs)
      final_cutoff <- as.numeric(names(sort(cutoff_table, decreasing = TRUE)[1]))
      
      if (final_cutoff == 0) {
        warning(paste("Optimal cutoff for", var, "is zero. Consider adjusting maxstat parameters."))
        final_cutoff <- mean(valid_data[[var]], na.rm = TRUE)
      }
      
      cat("Optimal cutoff for", var, ":", final_cutoff, "\n")
      
      new_column <- ifelse(data[[var]] > final_cutoff, paste0(">", round(final_cutoff, 2)), paste0("<=", round(final_cutoff, 2)))
      new_data[[var]] <- factor(new_column, levels = c(paste0("<=", round(final_cutoff, 2)), paste0(">", round(final_cutoff, 2))))
      
      cat("Finished processing variable:", var, "\n\n")
    } else {
      cat("Skipping non-numeric variable:", var, "\n")
    }
  }
  
  return(new_data)
}

perform_subgroup_analysis <- function(object,
                                      status_col="status",
                                      id_col = "id",
                                      time_col = "time",
                                      max_unique_values  = 5,
                                      methods = "roc",
                                      seed =123,
                                      save_dir = here::here("SurObj","Data"),
                                      save_data = TRUE,
                                      csv_filename = "subgroup_data.csv" ) {
  
  cat("Performing subgroup analysis...\n")
  
  if (inherits(object, 'SurObj')) {
    data <- slot(object, "survival.data")
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'SurObj' or a data frame.")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("The survival.data in the SurObj object is empty.")
  }
  
  cat("Diagnosing variable types...\n")
  diagnostics <- diagnose_variable_type(data, 
                                        max_unique_values=max_unique_values)
  
  cat("Processing variables by method: ", methods, "\n")
  
  if ("maxstat" %in% methods) {
    new_data <- find_maxstat_cutoff(data, 
                                    variables=diagnostics$numeric_vars,
                                    time_col =time_col,
                                    status_col = status_col)  } 
  else if ("roc" %in% methods) {
    new_data <- find_roc_cutoff(data, 
                                variables=diagnostics$numeric_vars, 
                                status_col = status_col, 
                                
                                time_col=time_col)  } else {
                                  stop("Invalid method specified. Choose either 'roc' or 'maxstat'.")
                                }
  
  variable_types <- diagnose_variable_type(new_data, max_unique_values=max_unique_values)
  new_data <- convert_two_class_to_binary(new_data)
  new_data <- convert_variables(new_data, variable_types = variable_types)
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(new_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  
  if (inherits(object, 'SurObj')) {
    cat("Updating 'SurObj' object...\n")
    
    object@sub.data <- new_data
    
    cat("The 'SurObj' object has been updated with the following slots:\n")
    cat("- 'sub.data' slot updated.\n")
    
    return(object)
  }
  
  return(new_data)
}



