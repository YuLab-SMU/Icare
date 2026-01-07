
diagnose_variable_type <- function(data,
                                   group_col = "group",
                                   max_unique_values = 5) {
  numeric_vars <- vector("list")
  categorical_vars <- vector("list")
  vars_to_encode <- vector()
  is_group_col_present <- !is.null(group_col) && group_col %in% names(data)
  for (col in names(data)) {
    if (!is_group_col_present || col != group_col) {
      unique_values <- length(unique(data[[col]]))
      if (unique_values <= max_unique_values) {
        categorical_vars[[col]] <- col
        if (unique_values > 2) {
          vars_to_encode <- c(vars_to_encode, col)
        }
      } else if (is.numeric(data[[col]])) {
        numeric_vars[[col]] <- col
      }
    }
  }
  numeric_vars <- unlist(numeric_vars)
  categorical_vars <- unlist(categorical_vars)
  return(list(numeric_vars = numeric_vars, categorical_vars = categorical_vars, vars_to_encode = vars_to_encode))
}


convert_two_class_to_binary <- function(data) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }
  
  original_row_names <- rownames(data)
  
  data_converted <- data.frame(lapply(data, function(column) {
    
    unique_values <- unique(column)
    
    if (length(unique_values) == 2 && all(unique_values %in% c(1, 2))) {
      column <- ifelse(column == 1, 0, 1)
    }
    return(column)
  }))
  
  rownames(data_converted) <- original_row_names
  
  return(data_converted)
}


convert_variables <- function(data, variable_types) {
  stopifnot(is.data.frame(data))
  for (col in names(data)) {
    if (col %in% variable_types$numeric_vars) {
      data[[col]] <- as.numeric(data[[col]])
      cat("Converted ", col, " to numeric.\n")
    } else {
      data[[col]] <- factor(data[[col]])
      cat("Converted", col, "to factor.\n")
    }
  }
  return(data)
}
convert_to_numeric <- function(data) {
  # Use apply() to apply as.numeric to all columns
  data[] <- lapply(data, function(x) as.numeric(x))
  return(data)
}

PrepareDataForPrognosiX <- function(object ,
                                    group_col = NULL,
                                    max_unique_values = 5) {
  cat("Starting transformation of factors to dummy variables... \n")

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting the data slot... \n")
    survival_data<-slot(object, "survival.data")
    status_col <- slot(object, "status_col")
    sub_data<-slot(object, "sub.data")
  }

  if (nrow(survival_data) == 0 && nrow(sub_data) == 0) {
    cat("Both of the extracted datasets are empty. Skipping transformation.\n")
    return(object)
  }

  survival_data <- convert_to_numeric(survival_data)
  sub_data <- convert_to_numeric(sub_data)

  survival_data <- convert_two_class_to_binary(survival_data)
  sub_data <- convert_two_class_to_binary(sub_data)

  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object with transformed data... \n")
    slot(object, "survival.data") <- survival_data
    slot(object, "sub.data") <-sub_data
    cat("The 'PrognosiX' object has been updated with the following:\n")
    cat("- 'survival.data' slot updated.\n")
    cat("- 'sub.data' slot updated.\n")

    return(object)
  }

  cat("Returning transformed data frame.  \n")
  return(data)
}
