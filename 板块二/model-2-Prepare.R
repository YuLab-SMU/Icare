#' Convert Factors to Binary or Dummy Variables
#'
#' Converts binary factors to 0/1 numeric values and multi-level factors to dummy variables.
#' Skips factors with fewer than 2 levels and provides informative warnings.
#'
#' @param data A data frame containing variables to convert
#' @return A data frame with converted variables
#' @export
convert_factors_to_binary <- function(data) {
  cat("Starting conversion of factors to binary or dummy variables...\n")

  original_row_names <- rownames(data)
  converted_columns <- list()

  for(col_name in colnames(data)) {
    column <- data[[col_name]]

    if(is.factor(column) || is.character(column)) {
      unique_vals <- unique(column[!is.na(column)])
      n_levels <- length(unique_vals)

      if(n_levels < 2) {
        warning("Skipping column '", col_name, "' - has fewer than 2 levels (",
                n_levels, ")")
        converted_columns[[col_name]] <- column
        next
      }

      if(n_levels == 2) {
        levels <- sort(unique_vals)
        converted <- rep(NA, length(column))
        converted[column == levels[1] & !is.na(column)] <- 0
        converted[column == levels[2] & !is.na(column)] <- 1
        converted_columns[[col_name]] <- converted
      } else {
        tryCatch({
          temp_df <- data.frame(column = column)
          mm <- model.matrix(~ column - 1, data = temp_df, na.action = na.pass)
          colnames(mm) <- gsub("^column", paste0(col_name, "_"), colnames(mm))
          converted_columns <- c(converted_columns, as.list(as.data.frame(mm)))
        }, error = function(e) {
          warning("Failed to create dummy variables for '", col_name, "': ", e$message)
          converted_columns[[col_name]] <- column
        })
      }
    } else {
      converted_columns[[col_name]] <- column
    }
  }

  data_converted <- as.data.frame(converted_columns)
  rownames(data_converted) <- original_row_names

  cat("Conversion completed.\n")
  return(data_converted)
}

#' Prepare Data for Modeling
#'
#' This function prepares raw data for machine learning by performing type conversion,
#' one-hot encoding, and binary encoding of categorical variables. It handles both
#' data frames and specialized model objects.
#'
#' @param object An input object of class 'Train_Model' or a data frame
#' @param group_col Name of the grouping column (default: "group")
#' @param max_unique_values Maximum unique values to consider a variable categorical (default: 5)
#' @param save_data Logical indicating whether to save output as CSV (default: FALSE)
#' @param csv_filename Name for output CSV file (default: "prepared_data.csv")
#' @param save_dir Directory path for saving CSV (default: here::here("ModelData","Data"))
#'
#' @return Depending on input:
#'   - For 'Train_Model' objects: Returns updated object with prepared data
#'   - For data frames: Returns prepared data frame
#' @export
#'
#' @examples
#' \dontrun{
#' # Using with a data frame
#' data(iris)
#' prepared_iris <- PrepareData(iris, group_col = "Species")
#'
#' # Using with Train_Model object
#' model <- new("Train_Model", data.df = iris)
#' object_model <- PrepareData(object_model)
#'
#' # Save output to CSV
#' PrepareData(object_model, save_data = TRUE)
#' }
#' @importFrom dplyr intersect
PrepareData <- function(object = NULL,
                        group_col = "group",
                        max_unique_values = 5,
                        save_data = TRUE,
                        csv_filename = "prepared_data.csv",
                        save_dir = here::here("ModelData","Data")) {

  if (inherits(object, 'Train_Model')) {
    cat("Input is a 'Train_Model' object. Extracting the data slot...\n")
    original_data <- slot(object, "data.df")
    return_object_modelect <- TRUE
  } else if (is.data.frame(object)) {
    cat("Input is a data frame.\n")
    original_data <- object
    return_object_modelect <- FALSE
  } else {
    stop("Input must be an object of class 'Train_Model' or a data frame.")
  }

  if (is.null(original_data) || nrow(original_data) == 0) {
    stop("No valid data found in the input.")
  }

  original_colnames <- colnames(original_data)

  variable_types <- diagnose_variable_type(original_data,
                                           group_col = group_col,
                                           max_unique_values = max_unique_values)


  
  # 检查group_col是否存在且有两个类别
  if (!is.null(group_col) && group_col %in% names(original_data)) {
    categories <- unique(na.omit(original_data[[group_col]]))
    
    if (length(categories) == 2) {
      # 定义映射规则
      mapping <- list(
        original = sort(categories),
        binary   = c(0, 1)
      )
      names(mapping$binary) <- mapping$original
      
      cat("\nMapping rules for '", group_col, "':\n")
      cat("1. Original category '", mapping$original[1], "' → 0\n")
      cat("2. Original category '", mapping[["original"]][2], "' → 1\n")
      cat("(Note: Categories are sorted alphabetically before assigning 0/1)\n")
      
      # 实际执行映射
      original_data[[group_col]] <- ifelse(
        original_data[[group_col]] == mapping$original[1], 0, 1
      )
      
      cat("\nActual mapping applied to '", group_col, "':\n")
      cat("- '", mapping$original[1], "' → ", 0, "\n")
      cat("- '", mapping$original[2], "' → ", 1, "\n")
      
    } else {
      cat("\nWarning: '", group_col, "' has", length(categories), 
          "categories, expected 2 for binary conversion.\n")
    }
  }
  
  converted_data <- convert_variables(original_data, 
                                      variable_types,
                                      group_col = group_col)
  
  
  onehot_data <- one_hot_encode(converted_data,
                                group_col = group_col,
                                max_unique_values = max_unique_values)
  
  converted_data <- convert_factors_to_binary(onehot_data)
  converted_data <- data.frame(converted_data)
  final_data<-converted_data
  common_cols <- dplyr::intersect(colnames(final_data), original_colnames)
  colnames(final_data)[colnames(final_data) %in% common_cols] <-
    original_colnames[original_colnames %in% common_cols]

  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("Created directory:", save_dir, "\n")
      }

      full_path <- file.path(save_dir, csv_filename)
      write.csv(final_data, file = full_path, row.names = FALSE)
      cat("\nData successfully saved as CSV file:", full_path, "\n")
    }, error = function(e) {
      warning("Failed to save CSV file: ", e$message)
    })
  }

  if (inherits(object, 'Train_Model')) {
    slot(object, "clean.df") <- final_data
    object@process.info[["mapping"]] <- mapping
    cat("Updating 'Train_Model' object...\n")
    cat("The 'Train_Model' object has been updated with the following slots:\n")
    cat("- 'clean.df' slot updated.\n")
    return(object)
  } else {
    cat("Returning transformed data frame with preserved column names.\n")
    return(final_data)
  }
}
