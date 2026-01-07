#' Plot Missing Data Distribution
#'
#' This function generates density plots to visualize the distribution of missing data in the given dataset.
#' It provides a variable-wise and sample-wise missing data distribution plot, as well as a combined plot.
#' The function also saves the plots as PNG files if specified and returns a summary of missing data statistics.
#
#' @import ggplot2
#' @import here
#' @import wesanderson
#' @import stats
#' @import here 
#' @importFrom ggprism theme_prism
#' @importFrom grDevices pdf
#' @param data A data.frame containing the data to be analyzed. Missing values should be represented as `<NA>` or `NA`.
#' @param palette_name A character string specifying the palette name for the plot colors. Default is `"Royal1"`.
#'                   Supported palettes are from the `wesanderson` package (e.g., "Royal1", "Zissou1", "GrandBudapest1").
#' @param alpha A numeric value between 0 and 1 specifying the transparency level of the density plots. Default is `0.9`.
#' @param save_plots A logical value indicating whether to save the plots as PNG files. Default is `TRUE`.
#' @param save_dir A character string specifying the directory where the plots will be saved. Default is `here("StatObject")`,
#'                which uses the `here` package to generate the path.
#' @param plot_width A numeric value specifying the width of the saved plots. Default is `5`.
#' @param plot_height A numeric value specifying the height of the saved plots. Default is `5`.
#' @param base_size A numeric value specifying the base size for text and elements in the plot. Default is `14`.
#'
#' @returns A list containing the following:
#'   - `var_plot_obj`: The ggplot object for the variable-wise missing data distribution plot.
#'   - `sample_plot_obj`: The ggplot object for the sample-wise missing data distribution plot.
#'   - `combined_plot`: The ggplot object for the combined missing data distribution plot.
#'   - `total_missing_values`: The total number of missing values in the data.
#'   - `total_variables`: The total number of variables (columns) in the data.
#'   - `total_samples`: The total number of samples (rows) in the data.
#'   - `variables_with_missing`: The number of variables with at least one missing value.
#'   - `samples_with_missing`: The number of samples with at least one missing value.
#'   - `max_missing_variable`: The highest percentage of missing values in any variable (column).
#'   - `min_missing_variable`: The lowest percentage of missing values in any variable (column).
#'   - `max_missing_sample`: The highest percentage of missing values in any sample (row).
#'   - `min_missing_sample`: The lowest percentage of missing values in any sample (row).
#'
#' @export
#'
#' @examples
#' # Generate missing data plots for the mtcars dataset
#' missing_info <- plot_missing_data(data = mtcars, save_plots = TRUE)
#'
#' # Generate missing data plots with customized parameters
#' missing_info <- plot_missing_data(data = mtcars, palette_name = "Zissou1", alpha = 0.7, plot_width = 6, plot_height = 6)
plot_missing_data <- function(data,
                              palette_name = 'Royal1',
                              alpha = 0.9,
                              save_plots = TRUE,
                              save_dir = here("StatObject"),
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 14,
                              save_data = TRUE,
                              var_filename = "var_missing_data.csv",
                              sample_filename = "sample_missing_data.csv") {

  colors <- wes_palette(n = 3, name = palette_name, type = "discrete")
  colors <- as.list(colors)

  primary_color <- colors[[1]]
  secondary_color <- colors[[2]]
  tertiary_color <- colors[[3]]

  data[data == '<NA>'] <- NA

  var_missing_percentage <- colMeans(is.na(data)) * 100
  sample_missing_percentage <- rowMeans(is.na(data)) * 100

  if (all(var_missing_percentage == 0) && all(sample_missing_percentage == 0)) {
    cat("No missing values in the data.")
    return(list(var_plot_obj = NULL, sample_plot_obj = NULL, combined_plot = NULL))
  }

  if (mean(var_missing_percentage) < 5 && mean(sample_missing_percentage) < 5) {
    cat("The percentage of missing data is low (below 5%).")
  }

  var_missing_df <- data.frame(Variable = names(var_missing_percentage), Missing_Percentage = var_missing_percentage)
  sample_missing_df <- data.frame(Sample = 1:nrow(data), Missing_Percentage = sample_missing_percentage)

  var_plot_obj <- ggplot(data = var_missing_df, aes(x = Missing_Percentage, fill = "Variable-wise")) +
    geom_density(alpha = alpha, color = NA) +
    labs(x = "Missing Percentage", y = "Density", title = "Variable Missing Data") +
    scale_fill_manual(values = primary_color) +
    theme_prism(base_size = base_size) +
    theme(legend.position = "top")

  sample_plot_obj <- ggplot(data = sample_missing_df, aes(x = Missing_Percentage, fill = "Sample-wise")) +
    geom_density(alpha = alpha, color = NA) +
    labs(x = "Missing Percentage", y = "Density", title = "Sample Missing Data") +
    scale_fill_manual(values = secondary_color) +
    theme_prism(base_size = base_size) +
    theme(legend.position = "none")

  combined_plot <- var_plot_obj +
    geom_density(data = sample_missing_df, aes(x = Missing_Percentage, fill = "Sample-wise"), alpha = alpha, color = NA) +
    labs(y = "Density", title = "Overall Missing Data Distribution") +
    scale_fill_manual(values = c(primary_color, secondary_color)) +
    guides(fill = guide_legend(title = NULL)) +
    theme_prism(base_size = base_size)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    cat("\nAll plots saved successfully to: \n",save_dir,"\n")
    
    ggsave(filename = file.path(save_dir, "variable_missing_data_plot.pdf"), plot = var_plot_obj,
           width = plot_width, height = plot_height, device = "pdf")
    ggsave(filename = file.path(save_dir, "sample_missing_data_plot.pdf"), plot = sample_plot_obj,
           width = plot_width, height = plot_height, device = "pdf")
    ggsave(filename = file.path(save_dir, "combined_missing_data_plot.pdf"), plot = combined_plot,
           width = plot_width, height = plot_height, device = "pdf")
  }
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path1 <- file.path(save_dir, var_filename)
    write.csv(var_missing_df, file = full_path1, row.names = FALSE)
    full_path2 <- file.path(save_dir, sample_filename)
    write.csv(sample_missing_df, file = full_path2, row.names = FALSE)
    cat("Saved variable missing rate data to:", full_path1, "\n")
    cat("Saved sample missing rate data to:", full_path2, "\n")
  }
  
  missing_info <- list(
    var_plot_obj = var_plot_obj,
    sample_plot_obj = sample_plot_obj,
    combined_plot = combined_plot,
    total_missing_values = sum(is.na(data)),
    total_variables = ncol(data),
    total_samples = nrow(data),
    variables_with_missing = sum(var_missing_percentage > 0),
    samples_with_missing = sum(sample_missing_percentage > 0),
    max_missing_variable = max(var_missing_percentage),
    min_missing_variable = min(var_missing_percentage),
    max_missing_sample = max(sample_missing_percentage),
    min_missing_sample = min(sample_missing_percentage)
  )
  print(combined_plot)
  return(missing_info)
}


#' Plot Missing Data for a Stat Object or Data Frame
#'
#' This function generates and saves missing data plots (variable-wise and sample-wise distributions) for
#' an object of class 'Stat' or a data frame. It also updates the 'Stat' object with missing data information
#' and returns the object or the missing data summary, depending on the input type.
#'
#' @import ggplot2
#' @import here
#' @import wesanderson
#' @import stats
#' @import methods
#' @import here 
#' @importFrom ggprism theme_prism
#' @importFrom grDevices pdf
#' @param object An object of class 'Stat' or a data frame. If it is a 'Stat' object, the function uses the
#'               `raw.data` slot for missing data analysis. If it is a data frame, it directly performs the analysis.
#' @param palette_name A character string specifying the palette name for the plot colors. Default is `"Royal1"`.
#' @param alpha A numeric value between 0 and 1 specifying the transparency level of the density plots. Default is `0.9`.
#' @param save_plots A logical value indicating whether to save the plots as PNG files. Default is `TRUE`.
#' @param save_dir A character string specifying the directory where the plots will be saved. Default is `"here('StatObject', 'missing_info')"` (within the `StatObject` folder).
#' @param plot_width A numeric value specifying the width of the saved plots. Default is `5`.
#' @param plot_height A numeric value specifying the height of the saved plots. Default is `5`.
#'
#' @returns If the input is a 'Stat' object, returns the updated 'Stat' object with the missing data information stored
#'          in the `process.info` slot. If the input is a data frame, returns the missing data summary as a list.
#'
#' @export
#'
#' @examples
#' # Generate missing data plots for a Stat object
#' updated_stat <- state_plot_missing_data(stat_object, save_plots = TRUE)
#'
#' # Generate missing data plots for a data frame
#' missing_info <- state_plot_missing_data(data_frame, save_plots = FALSE)
state_plot_missing_data <- function(
    object,
    palette_name = 'Royal1',
    alpha = 0.9,
    save_plots = TRUE,
    save_dir = here("StatObject","missing_info"),
    plot_width = 5,
    plot_height = 5,
    save_data = TRUE,
    var_filename = "var_missing_data.csv",
    sample_filename = "sample_missing_data.csv") {

  if (inherits(object, "Stat")) {
    data <- slot(object, "raw.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input.")
  }

  missing_info <- plot_missing_data(data,
                                    palette_name = palette_name,
                                    alpha = alpha,
                                    save_plots = save_plots,
                                    save_dir = save_dir,
                                    plot_width = plot_width,
                                    plot_height = plot_height,
                                    save_data = save_data,
                                    var_filename = var_filename,
                                    sample_filename = sample_filename)

  print(missing_info)

  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@process.info[["missing_count"]] <- missing_info
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'process.info' slot updated.\n")
    return(object)
  }
  return(missing_info)
}

#' Extract Raw Data from Stat Object
#'
#' This function extracts the 'raw.data' slot from an object of class 'Stat'.
#' If the object is not of class 'Stat' or does not contain a 'raw.data' slot,
#' it will return NULL.
#'
#' @param object An object of class 'Stat' which contains a slot named 'raw.data'.
#'               This should be a valid Stat object.
#'
#' @returns Returns the 'raw.data' slot of the Stat object if it exists.
#'          If the object does not have a 'raw.data' slot or the slot is empty,
#'          it returns NULL.
#'
#' @export
#'
#' @examples
#' # Assuming 'stat_object' is a valid Stat object
#' raw_data <- ExtractRawData(stat_object)
#'
#' # If the object does not have raw.data, it will return NULL
#' missing_data <- ExtractRawData(non_stat_object)
ExtractRawData <- function(object) {
  data <- tryCatch(slot(object, "raw.data"), error = function(e) NULL)
  return(data)
}


#' Diagnose Variable Types in Data
#'
#' This function analyzes a data frame to classify variables into numeric,
#' categorical, and variables that need encoding based on their unique values.
#' It also ensures that the specified grouping column (if provided) is excluded from
#' the analysis.
#'
#' @param data A data frame that contains the data to be analyzed. Each column represents a variable.
#' @param group_col A character string specifying the name of the grouping column (default is "group").
#'                  This column will be excluded from the analysis.
#' @param max_unique_values A numeric value specifying the maximum number of unique values a column can have
#'                           to be considered as categorical. Columns with fewer unique values than this threshold
#'                           will be treated as categorical variables (default is 5).
#'
#' @returns A list containing three elements:
#'   - `numeric_vars`: A character vector of numeric variables.
#'   - `categorical_vars`: A character vector of categorical variables.
#'   - `vars_to_encode`: A character vector of categorical variables that have more than 3 unique values
#'                       and should be encoded.
#'
#' @export
#'
#' @examples
#' # Example 1: Diagnose variables in a data frame
#' result <- diagnose_variable_type(data_frame, group_col = "group", max_unique_values = 5)
#' print(result)
#'
#' # Example 2: Diagnose variables without a grouping column
#' result_no_group <- diagnose_variable_type(data_frame, group_col = NULL)
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
        if (unique_values > 3) {
          vars_to_encode <- c(vars_to_encode, col)
        }
      } else if (is.numeric(data[[col]])) {
        numeric_vars[[col]] <- col
      }
    }
  }
  numeric_vars <- unlist(numeric_vars)
  categorical_vars <- unlist(categorical_vars)
  return(list(numeric_vars = numeric_vars,
              categorical_vars = categorical_vars,
              vars_to_encode = vars_to_encode))
}

#' Diagnose Variable Types for 'Stat' Objects or Data Frames
#'
#' This function analyzes the variable types (numeric, categorical, and those needing encoding)
#' of a data frame or an object of class "Stat". If the input is a "Stat" object, it extracts
#' the raw data and group column from the object. If the input is a data frame, it directly uses
#' the provided data for diagnosis. It updates the "Stat" object with the diagnosed variable types.
#'
#' @param object An object of class "Stat" or a data frame. If the object is of class "Stat",
#'               the raw data and group column will be extracted from the object.
#'               If it is a data frame, the function directly operates on the data.
#' @param group_col A character string specifying the name of the grouping column (default is "group").
#'                  This column will be excluded from the analysis if present.
#' @param max_unique_values A numeric value specifying the maximum number of unique values a column can have
#'                           to be considered as categorical. Columns with fewer unique values than this threshold
#'                           will be treated as categorical variables (default is 5).
#'
#' @returns If the input is a "Stat" object, the updated object with the diagnosed variable types
#'          in the "variable.types" slot. If the input is a data frame, a list containing:
#'   - `numeric_vars`: A character vector of numeric variables.
#'   - `categorical_vars`: A character vector of categorical variables.
#'   - `vars_to_encode`: A character vector of categorical variables that have more than 2 unique values
#'                       and should be encoded.
#'
#' @export
#'
#' @examples
#' # Example 1: Diagnose variables in a "Stat" object
#' stat_obj <- stat_diagnose_variable_type(stat_object, group_col = "group")
#' print(stat_obj)
#'
#' # Example 2: Diagnose variables in a data frame
#' result <- stat_diagnose_variable_type(data_frame)
#' print(result)
stat_diagnose_variable_type <- function(object,
                                        group_col = "group",
                                        max_unique_values = 5) {

  if (inherits(object, "Stat")) {
    group_col = slot(object, "group_col")
    if (length(group_col) == 0) {
      group_col <- NULL
    }
    data <- slot(object, "raw.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }


  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)


  cat("Diagnosed variable types:\n")
  cat("Numeric variables:", length(variable_types$numeric_vars), "\n")
  cat("Categorical variables:", length(variable_types$categorical_vars), "\n")

  if (length(variable_types$numeric_vars) == 0 && length(variable_types$categorical_vars) == 0) {
    stop("No valid variables found after variable type diagnosis")
  }


  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")

    object@variable.types <- variable_types
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'variable.types' slot updated.\n")
    return(object)
  }

  return(variable_types)
}



#' Gaze Analysis for Group Comparison
#'
#' This function performs gaze analysis for group comparison using the specified formula,
#' method, and settings. It provides a detailed result, optionally saving it as a Word document.
#' The gaze analysis can be performed based on group columns in the data, and various settings
#' are available to customize the output.
#'
#' @import here
#' @import officer
#' @import flextable
#' @import stats
#' @import methods
#' @import autoReg
#' @param data A data frame containing the data to analyze.
#' @param formula A formula specifying the model to fit. If NULL, the formula is automatically created
#' based on the provided group columns (default is NULL).
#' @param group_cols A vector of column names representing the grouping variables (default is NULL).
#' @param digits The number of digits to display for the result (default is 1).
#' @param show.p A logical value indicating whether to display p-values in the output (default is TRUE).
#' @param gaze_method An integer between 1 and 5 representing the gaze analysis method to use (default is 3).
#' @param save_word A logical value indicating whether to save the results to a Word document (default is TRUE).
#' @param save_dir The directory to save the Word document (default is the "PrognosiX/gaze_baseline" folder).
#'
#' @returns A data frame or matrix with the gaze analysis results. If `save_word` is TRUE, the results
#' are also saved as a Word document in the specified directory.
#' @export
#'
#' @examples
#' # Example 1: Performing gaze analysis with a formula and custom settings
#' result <- gaze_analysis(data = my_data,
#'                         formula = ~ group + age + gender,
#'                         digits = 2,
#'                         show.p = TRUE,
#'                         gaze_method = 3,
#'                         save_word = TRUE,
#'                         save_dir = "path/to/save")
#'
#' # Example 2: Using the default formula based on group columns
#' result <- gaze_analysis(data = my_data,
#'                         group_cols = c("group"),
#'                         digits = 1,
#'                         show.p = FALSE,
#'                         gaze_method = 1)
gaze_analysis <- function(data,
                          formula = NULL,
                          group_cols = NULL,
                          digits = 1,
                          show.p = TRUE,
                          gaze_method = 3,
                          save_word = TRUE,
                          save_dir = here("PrognosiX", "gaze_baseline")) {

  if (!is.data.frame(data)) stop("The input 'data' must be a data frame.")

  if (is.null(formula)) {
    if (!is.null(group_cols)) {
      if (length(group_cols) == 1 && group_cols %in% colnames(data)) {
        formula <- as.formula(paste(group_cols, "~ ."))
        cat("Using formula with one group:", deparse(formula), "\n")
      } else if (length(group_cols) > 1 && all(group_cols %in% colnames(data))) {
        formula <- as.formula(paste(paste(group_cols, collapse = " + "), "~ ."))
        cat("Using formula with multiple groups:", deparse(formula), "\n")
      } else {
        stop("Group columns not found in data.")
      }
    } else {
      formula <- as.formula("~ .")
      cat("Using default formula: ~ .\n")
    }
  } else if (!inherits(formula, "formula")) {
    stop("The input 'formula' must be a valid formula.")
  }

  if (!is.numeric(digits) || digits < 0 || digits != as.integer(digits))
    stop("The input 'digits' must be a non-negative integer.")

  if (!is.logical(show.p))
    stop("The input 'show.p' must be a logical value (TRUE or FALSE).")

  if (!is.numeric(gaze_method) || gaze_method < 1 || gaze_method > 5 || gaze_method != as.integer(gaze_method))
    stop("The input 'gaze_method' must be an integer between 1 and 5.")

  tryCatch({
    cat("Running gaze analysis with method:", gaze_method, "\n")
    result <- gaze(formula, data, digits = digits, show.p = show.p, method = gaze_method)

    cat("Result type:", class(result), "\n")
    if (is.data.frame(result) || is.matrix(result)) {
      result <- myft(result)
      cat("Gaze analysis completed successfully.\n")

      if (save_word) {
        doc <- read_docx()
        doc <- doc %>%
          body_add_flextable(result) %>%
          body_add_par("Gaze Analysis Results", style = "heading 1")

        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }

        word_filename <- file.path(save_dir, "gaze_analysis.docx")
        print(doc, target = word_filename)
        cat("Word file saved to:", word_filename, "\n")
      }

      return(result)
    } else {
      stop("The result is not a data frame or matrix.")
    }
  }, error = function(e) {
    stop("An error occurred while performing the gaze analysis: ", e$message)
  })
}

#' Statistical Gaze Analysis
#'
#' This function performs gaze analysis on a given dataset or `Stat` object, providing results for
#' group comparisons, statistical significance, and optionally saving the results as a Word document
#' and plots as images. It supports customized formulas and gaze methods for the analysis.
#'
#' @import here
#' @import officer
#' @import flextable
#' @import methods
#' @import stats
#' @param object An object of class 'Stat' or a data frame containing the data to analyze.
#' @param formula A formula specifying the model to fit. If NULL, a default formula is used (default is NULL).
#' @param group_col The column name representing the grouping variable (default is "group").
#' @param digits The number of digits to display for the result (default is 1).
#' @param show.p A logical value indicating whether to display p-values in the output (default is TRUE).
#' @param gaze_method An integer between 1 and 5 representing the gaze analysis method to use (default is 3).
#' @param save_word A logical value indicating whether to save the results to a Word document (default is TRUE).
#' @param save_dir The directory to save the Word document and plot images (default is the "StatObject/gaze_baseline" folder).
#'
#' @returns An updated 'Stat' object with the gaze analysis results if the input is of class 'Stat',
#'         or a data frame containing the results otherwise.
#' @export
#'
#' @examples
#' stat_object <- Stat$new(clean.data = your_data_frame, group_col = "group")
#' updated_stat_object <- stat_gaze_analysis(stat_object,
#'                                          formula = ~ group + age + gender,
#'                                          digits = 2,
#'                                          show.p = TRUE,
#'                                          gaze_method = 3,
#'                                          save_word = TRUE,
#'                                          save_plots = TRUE,
#'                                          save_dir = "path/to/save")
#'
#' # Example 2: Performing gaze analysis on a data frame
#' result <- stat_gaze_analysis(your_data_frame,
#'                              formula = ~ group + age + gender,
#'                              save_word = FALSE,
#'                              save_plots = TRUE)
#'
stat_gaze_analysis <- function(object,
                               formula = NULL,
                               group_col = "group",
                               digits = 1,
                               show.p = TRUE,
                               gaze_method = 3,
                               save_word = TRUE,
                               save_dir = here("StatObject", "gaze_baseline")) {

  cat("Starting stat_gaze_analysis function...\n")

  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    } 
    group_col <- slot(object, "group_col")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (is.null(data) || nrow(data) == 0)
    stop("No valid data found in the input.")

  if (!is.null(group_col) && !group_col %in% colnames(data))
    group_col <- NULL

  cat("Data prepared for gaze analysis. Number of rows:", nrow(data), "\n")

  result <- gaze_analysis(data,
                          formula = formula,
                          group_cols = group_col,
                          digits = digits,
                          show.p = show.p,
                          gaze_method = gaze_method,
                          save_word = save_word,
                          save_dir = save_dir)

  print(result)

  if (save_word) {
    cat("Saving results as Word document...\n")
    doc <- read_docx()
    doc <- doc %>%
      body_add_flextable(result) %>%
      body_add_par("Gaze Analysis Results", style = "heading 1")

    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }

    word_filename <- file.path(save_dir, "gaze_analysis.docx")
    print(doc, target = word_filename)
    cat("Word file saved to:", word_filename, "\n")
  }

  if (inherits(object, "Stat")) {
    object@baseline.table <- result
    cat("Updating 'Stat' object...\n")
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'baseline.table' slot updated.\n")
    return(object)
  }

  return(result)
}
