detect_and_mark_outliers <- function(data,
                                     group_col = "group",
                                     palette_name = "Royal1",
                                     save_plots = FALSE,
                                     save_dir = here("StatObject"),
                                     plot_display_num = 1,
                                     sub_var = NULL,
                                     plot_width = 5,
                                     plot_height = 5,
                                     base_size = 14,
                                     custom_ranges = NULL,
                                     max_unique_values = 5) {
  stopifnot(is.data.frame(data))
  
  original_rownames <- rownames(data)
  variable_type<-diagnose_variable_type(data)
  numeric_vars<-variable_type$numeric_vars
  numeric_data <- data[,numeric_vars]
  
  detect_outliers_with_range <- function(x, var_name = NULL) {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    
    if (!is.null(custom_ranges) && !is.null(var_name) && var_name %in% names(custom_ranges)) {
      lower_bound <- custom_ranges[[var_name]][1]
      upper_bound <- custom_ranges[[var_name]][2]
      custom_flag <- TRUE
    } else {
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      custom_flag <- FALSE
    }
    
    outliers <- x < lower_bound | x > upper_bound
    
    return(list(
      outliers = outliers,
      range_info = data.frame(
        variable = ifelse(is.null(var_name), NA, var_name),
        IQR = IQR,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        is_custom_range = custom_flag
      )
    ))
  }
  
  range_info_list <- list()
  outlier_vars <- c()
  
  for (var in names(numeric_data)) {
    result <- detect_outliers_with_range(numeric_data[[var]], var)
    if (any(result$outliers, na.rm = TRUE)) {
      outlier_vars <- c(outlier_vars, var)
      range_info_list[[var]] <- result$range_info
    }
  }
  
  if (length(outlier_vars) == 0) {
    cat("No variables with outliers found.")
    return(NULL)
  }
  
  if (!is.null(sub_var)) {
    outlier_vars <- intersect(outlier_vars, sub_var)
    range_info_list <- range_info_list[outlier_vars]
  }
  
  range_info_df <- bind_rows(range_info_list)
  
  cat("Detected outlier variables:", paste(outlier_vars, collapse = ","), "\n")
  
  for (var in outlier_vars) {
    result <- detect_outliers_with_range(data[[var]], var)
    data <- data %>%
      mutate(!!paste0(var, "_outlier") := result$outliers)
  }
  
  if (!is.null(group_col)) {
    data[[group_col]] <- as.factor(data[[group_col]])
  }
  data_with_outliers <- data %>% dplyr::select(all_of(outlier_vars), all_of(group_col))
  
  melted_data <- pivot_longer(data_with_outliers, cols = all_of(outlier_vars), names_to = "variable", values_to = "value")
  colors<-as.vector(wes_palette(palette_name))
  all_plots <- list()
  
  for (i in seq(1, length(outlier_vars), by = plot_display_num)) {
    selected_vars <- outlier_vars[i:min(i + plot_display_num - 1, length(outlier_vars))]
    
    if (!is.null(group_col)) {
      p <- ggplot(melted_data %>% dplyr::filter(variable %in% selected_vars), aes(x = variable, y = value, fill = !!sym(group_col))) +
        geom_boxplot(outlier.shape = 19, outlier.colour = "red", outlier.size = 1) +
        scale_fill_manual(values = wes_palette(palette_name)) +
        labs(title = "Boxplot of Variables with Outliers",
             x = "Variable",
             y = "Value") +
        theme_prism(base_size = base_size) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.position = "top",
          legend.title = element_blank()
        )
    } else {
      p <- ggplot(melted_data %>% dplyr::filter(variable %in% selected_vars), aes(x = variable, y = value)) +
        geom_boxplot(outlier.shape = 19, outlier.colour = "red", outlier.size = 1, fill = colors[1]) +
        labs(title = "Boxplot of Variables with Outliers",
             x = "Variable",
             y = "Value") +
        theme_prism(base_size = base_size) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12)
        )
    }
    
    all_plots[[length(all_plots) + 1]] <- p
  }
  cat("Generated ", length(all_plots), " plots for outlier visualization.\n")
  print(all_plots[[1]])
  
  if (save_plots) {
    outlier_dir <- file.path(save_dir)
    
    if (!is.null(sub_var) && length(sub_var) > 0) {
      sub_var_name <- paste(sub_var, collapse = "_")
      sub_dir <- file.path(outlier_dir, paste("sub_var_", sub_var_name, sep = ""))
      if (!dir.exists(sub_dir)) {
        dir.create(sub_dir, recursive = TRUE)
      }
      outlier_dir <- sub_dir
    } else {
      if (!dir.exists(outlier_dir)) {
        dir.create(outlier_dir, recursive = TRUE)
      }
    }
    
    for (i in 1:length(all_plots)) {
      ggsave(filename = file.path(outlier_dir, paste0("boxplot_outliers_batch_", i, ".pdf")),
             plot = all_plots[[i]],
             width = plot_width,
             height = plot_height, device = "pdf")
    }
    cat("Saved plots to: ", outlier_dir, "\n")
  }
  
  rownames(data) <- original_rownames
  
  return(list(
    data_marked = data,
    plots = all_plots,
    normal_ranges = range_info_df
  ))
  
}


stat_detect_and_mark_outliers <- function(object,
                                          save_outliers = TRUE,
                                          group_col = "group",
                                          palette_name = "Royal1",
                                          handle_method = c("replace", "remove", "keep", "capping"),
                                          lower_quantile = 0.05,
                                          upper_quantile = 0.95,
                                          save_plots = TRUE,
                                          save_dir = here("StatObject", "pre_outlier_info"),
                                          plot_display_num = 1,
                                          sub_var = NULL,
                                          plot_width = 5,
                                          plot_height = 5,
                                          base_size = 14,
                                          custom_ranges = NULL,
                                          max_unique_values = 5) {
  
  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    }
    group_col <- object@group_col
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
  
  cat("Processing data with", nrow(data), "rows.\n")
  
  outlier_result <- detect_and_mark_outliers(
    data = data,
    group_col = group_col,
    palette_name = palette_name,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_display_num = plot_display_num,
    sub_var = sub_var,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    custom_ranges = custom_ranges,
    max_unique_values = max_unique_values
  )
  
  
  
  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@process.info[["outlier_info"]][["detect_info"]]<- outlier_result
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'process.info' slot updated with pre-handling information\n")
    
    return(object)
  }else{
    return(list(
      marked_data = outlier_result$data_marked,
      handled_data = handled_data$cleaned_data,
      handling_method = handled_data$method,
      plots = outlier_result$plots,
      normal_ranges = outlier_result$normal_ranges,
      handling_details = handled_data[!names(handled_data) %in% "cleaned_data"]
    ))
  }
}


extract_outlier_data <- function(object,
                                 slot_name = "outlier.marked") {
  if (!inherits(object, "Stat")) {
    stop("Input must be an object of class 'Stat'.")
  }
  
  if (!slot_name %in% slotNames(object)) {
    stop(paste("Slot", slot_name, "does not exist in the object."))
  }
  
  outlier_info <- slot(object, slot_name)
  
  if (is.null(outlier_info) || !("data_marked" %in% names(outlier_info))) {
    stop("No marked outlier data found.")
  }
  
  data_marked <- outlier_info[["data_marked"]]
  return(data_marked)
}


handle_outliers_remove <- function(data) {
  cleaned_data <- data %>%
    dplyr::filter(dplyr::if_all(tidyselect::ends_with("_outlier"), ~ !.))
  
  return(list(cleaned_data = cleaned_data,
              method = "remove"))
}


handle_outliers_replace <- function(data) {
  
  replaced_values <- list()
  outlier_cols <- grep("_outlier$", colnames(data), value = TRUE)
  
  for (var in outlier_cols) {
    original_var <- sub("_outlier$", "", var)
    
    median_value <- median(data[[original_var]], na.rm = TRUE)
    
    data[[original_var]][data[[var]]] <- median_value
    
    replaced_values[[original_var]] <- median_value
  }
  
  cleaned_data <- data %>% dplyr::select(-ends_with("_outlier"))
  
  return(list(
    cleaned_data = cleaned_data,
    method = "replace",
    replaced_values = replaced_values  
  ))
}

handle_outliers_keep <- function(data) {
  return(list(cleaned_data=data,
              method ="keep"))}



handle_outliers_capping <- function(data, 
                                    lower_quantile = 0.05, 
                                    upper_quantile = 0.95) {
  capping_bounds <- list()
  
  outlier_cols <- grep("_outlier$", colnames(data), value = TRUE)
  
  for (var in outlier_cols) {
    original_var <- sub("_outlier$", "", var)
    
    lower_bound <- quantile(data[[original_var]], lower_quantile, na.rm = TRUE)
    upper_bound <- quantile(data[[original_var]], upper_quantile, na.rm = TRUE)
    
    capping_bounds[[original_var]] <- list(
      lower_bound = lower_bound,
      upper_bound = upper_bound
    )
    
    data[[original_var]] <- pmax(
      pmin(data[[original_var]], upper_bound), 
      lower_bound
    )
  }
  
  cleaned_data <- data %>% dplyr::select(-ends_with("_outlier"))
  
  return(list(
    cleaned_data = cleaned_data,
    method = "capping",
    capping_bounds = capping_bounds,  
    parameters = list(
      lower_quantile = lower_quantile,
      upper_quantile = upper_quantile
    )
  ))
}

handle_outliers <- function(data,
                            handle_method = c("replace", "remove", "keep", "capping"),
                            lower_quantile = 0.05,
                            upper_quantile = 0.95,
                            save_dir = here::here("StatObject","Data"),
                            save_data = TRUE,
                            csv_filename = "clean_data.csv") {
  
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }
  
  if (!any(grepl("_outlier$", colnames(data)))) {
    stop("No '_outlier' columns found in the data.")
  }
  
  handle_method <- match.arg(handle_method)
  
  if (handle_method == "remove") {
    data <- handle_outliers_remove(data)
  } else if (handle_method == "replace") {
    data <- handle_outliers_replace(data)
  } else if (handle_method == "keep") {
    data <- handle_outliers_keep(data)
  } else if (handle_method == "capping") {
    data <- handle_outliers_capping(data,
                                    lower_quantile = lower_quantile,
                                    upper_quantile = upper_quantile)
  } else {
    stop("Unknown method!")
  }
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  return(data)
}


stat_handle_outliers <- function(object,
                                 save_cleaned = TRUE,
                                 handle_method = c("replace", "remove", "keep", "capping"),
                                 lower_quantile = 0.05,
                                 upper_quantile = 0.95,
                                 group_col = "group",
                                 save_dir = here::here("StatObject","Data"),
                                 save_data = TRUE,
                                 csv_filename = "clean_data.csv") {
  
  handle_method <- match.arg(handle_method, choices = c("replace", "remove", "keep", "capping"))
  
  if (inherits(object, "Stat")) {
    cat("'object' is of class 'Stat'.\n Extracting outlier data...\n")
    data <-  object@process.info[["outlier_info"]][["detect_info"]][["data_marked"]]
    group_col <- object@group_col
    if (length(group_col) == 0) {
      group_col <- NULL
    }
  } else if (is.data.frame(object)) {
    cat(" 'object' is a data frame.\n")
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input.")
  }
  
  cat("Handling outliers with method: ", handle_method,"\n")
  original_rownames <- rownames(data)
  handle_outliers_result <- handle_outliers(data,
                                  handle_method = handle_method,
                                  save_dir =save_dir,
                                  save_data = save_data,
                                  csv_filename = csv_filename,
                                  lower_quantile = lower_quantile,
                                  upper_quantile = upper_quantile)
  
  cleaned_data<-handle_outliers_result[["cleaned_data"]]
  rownames(cleaned_data) <- original_rownames
  cat("Outlier handling completed. Cleaned data has ", nrow(cleaned_data), " rows.","\n")
  
  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@process.info[["outlier_info"]][["handle_info"]] <- handle_outliers_result
    object@clean.data <- cleaned_data
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'outlier.marked' slot updated.\n")
    cat("- 'process.info' slot updated.\n")
    
    return(object)
  }
  
  return(cleaned_data)
}


