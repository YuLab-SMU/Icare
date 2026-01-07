process_new_data <- function(
    object,
    new_data,
    group_col = "group",
    miss_threshold = 25,
    save_plots = FALSE,
    save_dir = here::here("ModelData", "new_data"),
    plot_width = 5,
    plot_height = 5,
    base_size = 14,
    palette_name = "Royal1",
    data_type = "clean",
    save_data = FALSE,
    csv_filename = "new_data.csv",
    plot_display_num = 1,
    max_unique_values = 5,
    suppress_warnings = TRUE) {
  
  if (!inherits(object, "Stat")) {
    stop("Input object must be a Stat object. Received: ", class(object)[1])
  }
  
  if (is.null(new_data) || nrow(new_data) == 0) {
    stop("new_data must be a non-empty data frame. Current dimensions: ", 
         paste(dim(new_data), collapse = " x "))
  }
  reference_vars <- colnames(object@clean.data)
  # 2. 初始化设置 ----
  cat("\n=== Starting New Data Processing ===\n")
  
  processed_data <- new_data
  cat("\n--- Variable Processing ---\n")
  
  processed_data <- stat_onehot_encode(object=processed_data, 
                                       group_col = group_col, 
                                       max_unique_values = max_unique_values,
                                       save_dir =save_dir,
                                       save_data = save_data,
                                       csv_filename = csv_filename)
  
  processed_data<-stat_convert_variables(object =processed_data,
                                         group_col = group_col, 
                                         max_unique_values = max_unique_values,
                                         save_dir =save_dir,
                                         save_data = save_data,
                                         csv_filename = csv_filename)
  processed_data <- processed_data[, reference_vars, drop = FALSE]
  
  cat("\n--- Missing Value Processing ---\n")
  if (!is.null(object@process.info[["missing_removal"]])) {
    missing_removal <- object@process.info[["missing_removal"]]
    missing_info <- object@process.info[["missing_info"]]
    
    processed_data <- as.data.frame(lapply(processed_data, function(x) {
      x[x %in% c("<NA>", "NA", "", "NULL", "NaN")] <- NA
      x
    }), stringsAsFactors = FALSE)
    
    miss_threshold <- missing_removal$miss_threshold %||% miss_threshold
    missing_ratio <- rowMeans(is.na(processed_data)) * 100
    data_rm <- which(missing_ratio >= miss_threshold)
    
    if (length(data_rm) > 0) {
      cat("- Removed", length(data_rm), "samples with missingness >=", 
          miss_threshold)
      processed_data <- processed_data[-data_rm, , drop = FALSE]
    } else {
      cat("- No samples exceeded missingness threshold (max =", 
          max(missing_ratio, na.rm = TRUE), "%)\n")
    }
    }

    if (!is.null(missing_info$imputation_info)) {
      processed_data <- suppressWarnings(
        tryCatch({
        suppressWarnings(
            apply_imputation(new_data = processed_data, 
                             imputation_info=missing_info$imputation_info)
          )
          
        }, error = function(e) {
          cat("Imputation failed with error: ", e$message, "\nReturning data without imputation")
          processed_data
        })
      )
     
      cat("\n- Applied imputation method:", missing_info$imputation_method %||% "unknown", "\n")
    } else {
      cat("No imputation information found in the object. Skipping imputation.")
    }
  
  
  cat("\n--- Outlier Processing ---\n")
  if (!is.null(object@process.info[["outlier_info"]])) {
    outlier_info <- object@process.info[["outlier_info"]]
    normal_ranges_list <- outlier_info$detect_info$normal_ranges
    
    if (!is.null(outlier_info$detect_info)) {
      outlier_result <- tryCatch({
        detect_and_mark_outliers(
          data = processed_data,
          group_col = group_col,
          palette_name = palette_name,
          save_plots = save_plots,
          save_dir = save_dir,
          plot_display_num = plot_display_num,
          plot_width = plot_width,
          plot_height = plot_height,
          base_size = base_size,
          custom_ranges = normal_ranges_list,
          max_unique_values = max_unique_values
        )
      }, error = function(e) {
        cat("Outlier detection failed: ", e$message)
        list(data_marked = processed_data)
      })
      
      
    } else {
      cat("No outlier detection parameters found. Skipping outlier detection.")
    }
    data_marked <- outlier_result$data_marked
    
    if (!is.null(outlier_info$handle_info)) {
        handle_result <-  tryCatch({apply_outlier_handling(
          data_marked=data_marked,
          handling_result = outlier_info$handle_info
        )
      
        
      }, error = function(e) {
        cat("Outlier handling failed: ", e$message, "\n")
        data_marked  
      })
      cat("- Applied outlier handling:", outlier_info$handle_info$method %||% "unknown", "\n")
    } else {
      cat("No outlier handling method specified. Skipping outlier handling.")
    }
  } else {
    cat("Outlier processing information not found in the object. Skipping outlier processing.")
  }
  
  cat("\n--- Data Normalization ---\n")
  if (data_type == "scale" && !is.null(object@process.info[["normalization.info"]])) {
    norm_info <- object@process.info[["normalization.info"]]
    processed_data <- tryCatch({
      norm_data  <- apply_previous_normalization(
        new_data = processed_data,
        prev_normalization_info = norm_info,
        normalize_method = norm_info[["normalize_method"]],
        group_col = group_col)
      norm_data[["scaled_data"]]
    }, error = function(e) {
      cat("Normalization failed: ", e$message, "\nReturning unscaled data")
      processed_data
    })
    cat("- Applied normalization method:", norm_info$normalize_method %||% "unknown", "\n")
  } else if (data_type == "scale") {
    cat("Normalization requested but no normalization info found. Returning unscaled data.")
  }
  
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      cat("- Created output directory:", save_dir, "\n")
    }
    
    data_path <- file.path(save_dir, csv_filename)
    write.csv(processed_data, data_path, row.names = TRUE)
    cat("- Saved processed data to:", normalizePath(data_path), "\n")
  }
  
  cat("\n=== Processing Completed ===\n")

  
  return(processed_data)
}
