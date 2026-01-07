#' Modify column names of a data frame
#'
#' This function modifies the column names of a data frame by performing the following:
#' 1. Converts all column names to lowercase.
#' 2. Replaces spaces with underscores.
#' 3. Removes any non-alphanumeric characters except for underscores.
#' 4. Ensures all column names are unique by appending numbers to duplicated names.
#'
#' @param data A data frame whose column names need to be modified.
#'
#' @returns A data frame with modified column names.
#' @export
#'
#' @examples
#' data <- data.frame("First Name" = c("John", "Jane"), "Age" = c(25, 30))
#' modified_data <- modify_column_names(data)
#' print(modified_data)
modify_column_names <- function(data) {
  if (is.null(names(data))) {
    stop("Column names of the data frame cannot be NULL!")
  }
  new_names <- names(data)
  new_names <- gsub(" ", "_", new_names)
  new_names <- gsub("[^[:alnum:]_]", "", new_names)

  if (anyDuplicated(new_names)) {
    warning("Duplicate column names found; they have been made unique.")
    new_names <- make.unique(new_names)
  }

  names(data) <- new_names
  return(data)
}



#' Class to store statistical analysis results
#'
#' This class is designed to store and manage various components of statistical analysis.
#' It includes raw and cleaned data, processed results, descriptive statistics, correlation results,
#' variable types, and metadata related to the analysis.
#'
#' @import methods
#' @slot raw.data A data.frame containing the raw input data.
#' @slot clean.data A data.frame containing the cleaned data after preprocessing.
#' @slot info.data A data.frame containing additional information related to the data.
#' @slot scale.data A data.frame containing the scaled version of the data.
#' @slot meta.featurename A character vector containing metadata feature names.
#' @slot group_col The grouping column in the data, which can be of any type.
#' @slot baseline.table A table that stores baseline information, can be of any type.
#' @slot process.info A list containing information about the processing steps.
#' @slot compute.descriptive A list containing results of descriptive statistics computations.
#' @slot corr.result A list containing the results of correlation analysis.
#' @slot var.result A list containing the results of variance analysis.
#' @slot variable.types A list specifying the types of the variables in the dataset.
#'
#' @returns An object of class 'Stat' containing the slots mentioned above.
#' @export
#' @examples
#' stat_obj <- new("Stat", raw.data = my_raw_data, clean.data = my_clean_data)
Stat <- setClass(
  Class = "Stat",
  slots = c(
    raw.data = "data.frame",
    clean.data = "data.frame",
    info.data = "data.frame",
    scale.data = "data.frame",
    meta.featurename = "character",
    group_col = "ANY",
    baseline.table = "ANY",
    process.info = "list",
    compute.descriptive = "list",
    corr.result = "list",
    var.result = "list",
    variable.types = "list"
  ),
  prototype = list(
    raw.data = data.frame(),
    clean.data = data.frame(),
    info.data = data.frame(),
    scale.data = data.frame(),
    meta.featurename = character(0),  
    group_col = NULL,
    baseline.table = NULL,
    process.info = list(),
    compute.descriptive = list(),
    corr.result = list(),
    var.result = list(),
    variable.types = list()
  )
)

#' Create a Stat object for statistical analysis
#'
#' This function creates a `Stat` object, which is used to store and manage various components
#' of statistical analysis, such as raw data, cleaned data, additional metadata, and processing
#' information. It performs basic checks and preparation on the input data before creating the object
#' @import methods
#' @param raw.data A data.frame containing the raw data for analysis. Defaults to an empty data frame.
#' @param clean.data A data.frame containing the cleaned data for analysis. Defaults to an empty data frame.
#' @param info.data A data.frame containing additional metadata related to the data. Defaults to an empty data frame.
#' @param group_col A character string specifying the name of the column used for grouping. Default is `"group"`.
#' @param ... Additional arguments passed to methods (not used in this function).
#'
#' @returns An object of class `Stat`, which contains the processed data and metadata as slots.
#' @export
#'
#' @examples
#' # Creating a Stat object with example data
#' stat_obj <- CreateStatObject(raw.data = example_raw_data, clean.data = example_clean_data)
#'
#' # Creating a Stat object with metadata
#' stat_obj <- CreateStatObject(raw.data = example_raw_data, clean.data = example_clean_data, info.data = example_info_data)
CreateStatObject <- function(
    raw.data = data.frame(),
    clean.data = data.frame(),
    info.data = data.frame(),
    group_col = "group",
    ...
) {

  raw.data <- as.data.frame(raw.data)
  clean.data <- as.data.frame(clean.data)
  info.data <- as.data.frame(info.data)

  clean_symbol_values <- function(data) {
    for(col in colnames(data)) {
      if(is.character(data[[col]])) {

        if(any(grepl("^[<>]", data[[col]]))) {

          data[[col]] <- gsub("[<>]", "", data[[col]])

          data[[col]] <- as.numeric(data[[col]])
          warning(paste("Removed >/< symbols from column", col, "and converted to numeric"))
        }
      }
    }
    return(data)
  }

  if (nrow(raw.data) == 0 && nrow(clean.data) == 0) {
    stop("At least one of 'raw.data' or 'clean.data' must be provided and not empty.")
  }


  
  if(nrow(raw.data) > 0) {
    raw.data <- clean_symbol_values(raw.data)
  }
  if(nrow(clean.data) > 0) {
    clean.data <- clean_symbol_values(clean.data)
  }

  if (nrow(info.data) > 0) {
    if (nrow(raw.data) > 0 && !setequal(rownames(info.data), rownames(raw.data))) {
      stop("Row names in 'info.data' do not match 'raw.data' (content mismatch).")
    }
    
    if (nrow(clean.data) > 0 && !setequal(rownames(info.data), rownames(clean.data))) {
      stop("Row names in 'info.data' do not match 'clean.data' (content mismatch).")
    }
  }

  prepare_data <- function(data, data_name) {
    if (nrow(data) == 0) {
      return(data)
    }

    if (anyDuplicated(rownames(data))) {
      warning(paste("Duplicate row names found in", data_name, "; they have been made unique."))
      rownames(data) <- make.unique(rownames(data))
    }

    if (anyDuplicated(colnames(data))) {
      warning(paste("Duplicate column names found in", data_name, "; they have been made unique."))
      colnames(data) <- make.unique(colnames(data))
    }

    if (is.null(colnames(data))) {
      stop(paste(data_name, "is missing column names."))
    }

    data <- modify_column_names(data)

    cat(paste("Data prepared for", data_name, "\n"))

    return(data)
  }

  raw.data <- prepare_data(raw.data, "raw.data")
  clean.data <- prepare_data(clean.data, "clean.data")

  if (nrow(raw.data) > 0) {
    meta.featurename <- as.character(colnames(raw.data))
  } else if (nrow(clean.data) > 0) {
    meta.featurename <- as.character(colnames(clean.data))
  } else {
    meta.featurename <- character()
  }

  if (nrow(info.data) == 0) {
    if (nrow(clean.data) > 0) {
      info.data <- data.frame(row.names = row.names(clean.data))
    } else if (nrow(raw.data) > 0) {
      info.data <- data.frame(row.names = row.names(raw.data))
    }
  } else {
    if (nrow(clean.data) > 0) {
      rownames(info.data) <- rownames(clean.data)
    } else if (nrow(raw.data) > 0) {
      rownames(info.data) <- rownames(raw.data)
    }
  }

  cat("Final info.data prepared.\n")

  Stat <- new(
    Class = 'Stat',
    raw.data = raw.data,
    info.data = info.data,
    clean.data = clean.data,
    meta.featurename = meta.featurename,
    process.info = list(),
    group_col = group_col
  )

  cat("Stat object created.\n")

  return(Stat)
}

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


