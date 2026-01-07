convert_variables <- function(data, 
                              variable_types,
                              save_dir = here::here("StatObject","Data"),
                              save_data = F,
                              group_col =NULL,
                              csv_filename = "clean_data.csv") {
  stopifnot(is.data.frame(data))
  all_factor_cols <- names(data)
  
  if (!is.null(group_col) && group_col %in% all_factor_cols) {
    all_factor_cols <- all_factor_cols[all_factor_cols != group_col]
  }
  for (col in all_factor_cols) {
    if (col %in% variable_types$categorical_vars) {
      data[[col]] <- factor(data[[col]])
      cat("Converted", col, "to factor.\n")
    } else {
      data[[col]] <- as.numeric(data[[col]])
      cat("Converted ", col, " to numeric.\n")
    }
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



one_hot_encode <- function(data, 
                           group_col = "group",
                           max_unique_values = 5,
                           save_dir = here::here("StatObject","Data"),
                           save_data = TRUE,
                           csv_filename = "clean_data.csv") {
  
  if (!is.data.frame(data)) stop("Input must be a data frame")
  
  if (length(group_col) == 0 || !is.character(group_col) || !(group_col %in% colnames(data))) {
    cat("Group column is not valid, setting to NULL.\n")
    group_col <- NULL
  }
  
  if (!is.numeric(max_unique_values) || max_unique_values <= 0) {
    stop("max_unique_values must be a positive numeric value")
  }
  
  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)
  vars_to_encode <- variable_types$vars_to_encode
  
  encoded_data <- data
  row_names <- rownames(data)
  
  for (var in vars_to_encode) { 
    unique_values <- unique(data[!is.na(data[, var]), var])
    cat("Encoding variable:", var, "with unique values:", paste(unique_values, collapse = ", "), "\n")
    for (value in unique_values) {
      col_name <- paste(var, value, sep = "_")
      encoded_data[, col_name] <- as.integer(data[, var] == value)
    }
    encoded_data[, var] <- NULL
  }
  
  if (!is.null(group_col) && group_col %in% names(encoded_data)) {
    group <- encoded_data[[group_col]]
    encoded_data[[group_col]] <- NULL
    encoded_data <- cbind(encoded_data, group)
    colnames(encoded_data)[ncol(encoded_data)] <- group_col
  }
  
  rownames(encoded_data) <- row_names
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    full_path <- file.path(save_dir, csv_filename)
    write.csv(encoded_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  return(encoded_data)
}
