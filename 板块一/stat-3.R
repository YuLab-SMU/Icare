#' Extract scaled data from an object
#'
#' This function extracts the `scale.data` slot from an object, which typically contains scaled data.
#' If the object does not have the `scale.data` slot, it will return `NULL`.
#'
#' @param object An object that may contain a `scale.data` slot.
#'
#' @returns A data object (e.g., a data frame) containing the scaled data from the `scale.data` slot, or `NULL` if the slot does not exist.
#' @export
#'
#' @examples
#' # Assuming 'stat_object' is an object with a 'scale.data' slot
#' scaled_data <- ExtractScaleData(stat_object)
#' head(scaled_data)
ExtractScaleData <- function(object) {
  data <- tryCatch(slot(object, "scale.data"), error = function(e) NULL)
  return(data)
}

#' Extract clean data from an object
#'
#' This function extracts the `clean.data` slot from an object, which typically contains the cleaned data.
#' If the object does not have the `clean.data` slot, it will return `NULL`.
#'
#' @param object An object that may contain a `clean.data` slot.
#'
#' @returns A data object (e.g., a data frame) containing the cleaned data from the `clean.data` slot, or `NULL` if the slot does not exist.
#' @export
#'
#' @examples
#' # Assuming 'stat_object' is an object with a 'clean.data' slot
#' clean_data <- ExtractCleanData(stat_object)
#' head(clean_data)
ExtractCleanData <- function(object) {
  data <- tryCatch(slot(object, "clean.data"), error = function(e) NULL)
  return(data)
}


#' Identify Highly Correlated Features and Generate Heatmap
#'
#' This function identifies highly correlated features in a data frame and generates a heatmap to visualize
#' the correlation matrix. It also returns a list of the correlated features and the top 5 most correlated feature pairs.
#' The function can optionally save the heatmap plot to a specified directory.
#'
#'
#' @import stats
#' @import wesanderson
#' @importFrom pheatmap pheatmap
#' @import here 
#' @param data A data frame containing the features to be analyzed.
#' @param correlation_threshold A numeric value specifying the threshold for correlation. Features with a correlation
#' greater than this value will be considered correlated (default is 0.95).
#' @param palette_name A string specifying the color palette to be used for the heatmap
#' @param save_plots A logical value indicating whether to save the heatmap plot (default is TRUE).
#' @param save_dir The directory where the heatmap plot will be saved, if `save_plots` is TRUE (default is `here("StatObject", "correlation_info")`).
#' @param heatmap_title A string specifying the title of the heatmap (default is "Correlation Heatmap").
#'
#' @returns A list containing the following components:
#' - `corr_matrix`: The correlation matrix of the features.
#' - `correlated_features`: A character vector of the features that are highly correlated.
#' - `record_collinear`: A data frame with pairs of correlated features and their correlation values.
#' - `top5_pairs`: A data frame containing the top 5 most correlated feature pairs.
#' - `original_data`: The original data frame passed to the function.
#'
#' @export
#'
#' @examples
#' # Assuming 'data' is a data frame of features
#' result <- highly_correlated_features(data, correlation_threshold = 0.95, save_plots = TRUE)
#' head(result$top5_pairs)
#'
highly_correlated_features <- function(data,
                                       correlation_threshold = 0.95,
                                       palette_name = "Zissou1",
                                       save_plots = TRUE,
                                       save_dir = here("StatObject", "correlation_info"),
                                       heatmap_title = "Correlation Heatmap",
                                       plot_width = 10,
                                       plot_height = 8,
                                       save_data = TRUE,
                                       csv_filename = "correlation_heatmap_data.csv"
) {
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }

  data <- as.data.frame(lapply(data, as.numeric))

  if (any(is.na(data))) {
    warning("Missing values found in data. Removing rows with missing values.")
    data <- na.omit(data)
  }

  if (any(sapply(data, is.infinite))) {
    warning("Infinite values found in data. Replacing with NA.")
    data <- data.frame(lapply(data, function(x) {
      if (is.numeric(x)) {
        x[is.infinite(x)] <- NA
      }
      return(x)
    }))
  }

  corr_matrix <- cor(data, use = "pairwise.complete.obs")

  cat("Maximum correlation value:", max(corr_matrix, na.rm = TRUE), "\n")
  cat("Minimum correlation value:", min(corr_matrix, na.rm = TRUE), "\n")

  record_collinear <- data.frame(drop_feature = character(), corr_feature = character(), corr_value = numeric())

  for (i in 1:nrow(corr_matrix)) {
    for (j in 1:ncol(corr_matrix)) {
      if (i != j && corr_matrix[i, j] > correlation_threshold) {
        temp_df <- data.frame(drop_feature = rownames(corr_matrix)[i],
                              corr_feature = colnames(corr_matrix)[j],
                              corr_value = corr_matrix[i, j])
        record_collinear <- rbind(record_collinear, temp_df)
      }
    }
  }

  sorted_record_collinear <- record_collinear[order(abs(record_collinear$corr_value), decreasing = TRUE), ]
  top5_pairs <- head(sorted_record_collinear, 5)

  cat("Top 5 correlated feature pairs:\n")
  print(top5_pairs)

  plot_clustermap_heatmap <- function(corr_matrix, palette_name, title) {
    library(wesanderson)
    color_palette <- wes_palette(palette_name, 100, type = "continuous")

    p <- pheatmap(
      corr_matrix,
      color = color_palette,
      breaks = seq(-1, 1, length.out = 100),
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      main = title
    )

    return(p)
  }

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }

    while (!is.null(dev.list())) dev.off()
    pdf(file.path(save_dir, "correlation_heatmap.pdf"), width = plot_width, height = plot_height)

    p <- plot_clustermap_heatmap(corr_matrix, palette_name, heatmap_title)
    print(p)
    dev.off()

    cat("Heatmap saved to:", file.path(save_dir, "correlation_heatmap.pdf"), "\n")
  } else {
    while (!is.null(dev.list())) dev.off()

    p <- plot_clustermap_heatmap(corr_matrix, palette_name, heatmap_title)
    print(p)
  }
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(corr_matrix, file = full_path, row.names = FALSE)
    cat("Correlation matrix data saved to:", full_path, "\n")  
  }
  
  result <- list(
    corr_matrix = corr_matrix,
    correlated_features = unique(record_collinear$drop_feature),
    record_collinear = record_collinear,
    top5_pairs = top5_pairs,
    original_data = data
  )

  return(result)
}




#' Identify Highly Correlated Features from a 'Stat' Object or Data Frame
#'
#' This function performs correlation analysis on a given data frame or 'Stat' object, identifying
#' features that are highly correlated based on a specified threshold. It returns the result of
#' the correlation analysis and optionally updates the 'Stat' object with the correlation results.
#'
#' @import stats
#' @import wesanderson
#' @importFrom pheatmap pheatmap
#' @import here 
#' @import methods
#' @param object An object of class 'Stat' or a data frame containing the features to be analyzed.
#' @param correlation_threshold A numeric value specifying the threshold for correlation. Features with a correlation
#' greater than this value will be considered correlated (default is 0.95).
#' @param palette_name A string specifying the color palette to be used for the heatmap
#' @param data_type A string indicating which type of data to extract from the 'Stat' object. Options are "clean" or "scale" (default is "scale").
#' @param save_dir The directory where the heatmap plot will be saved, if `save_plots` is TRUE (default is `here("StatObject", "correlation_info")`).
#'
#' @returns If the input is a 'Stat' object, the function returns the updated 'Stat' object with the correlation results added to the 'corr.result' slot.
#' If the input is a data frame, the function returns the results of the correlation analysis, including the correlation matrix,
#' identified correlated features, and other relevant data.
#'
#' @export
#'
#' @examples
#' # Assuming 'stat_object' is an instance of 'Stat' class
#' result <- stat_correlated_features(stat_object, correlation_threshold = 0.9)
#'
#' # If input is a data frame
#' result <- stat_correlated_features(my_data_frame, correlation_threshold = 0.85)
#'

stat_correlated_features <- function(object,
                                     correlation_threshold = 0.95,
                                     palette_name = "Zissou1",
                                     data_type = "scale",
                                     save_plots = TRUE,
                                     save_dir = here("StatObject", "correlation_info"),
                                     save_data = TRUE,
                                     csv_filename = "correlation_heatmap_data.csv") {
  if (inherits(object, "Stat")) {
    cat("Extracting data from 'Stat' object...\n")
    data <- if (data_type == "clean") ExtractCleanData(object) else ExtractScaleData(object)
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input.")
  }

  cat("Running correlation analysis with threshold:", correlation_threshold, "\n")
  re_cor <- highly_correlated_features(data,
                                       correlation_threshold = correlation_threshold,
                                       palette_name = palette_name,
                                       save_dir = save_dir,
                                       save_data = save_data,
                                       csv_filename = csv_filename)

  if (inherits(object, "Stat")) {
    object@corr.result <- re_cor
    cat("Updating 'Stat' object...\n")
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'corr.result' slot updated.\n")
    return(object)
  }

  return(re_cor)
}




#' Plot Top Correlated Features from a Data Frame or Matrix
#'
#' This function calculates the pairwise Pearson correlations of features in the provided
#' data frame or matrix, and plots the top `n` most correlated features using a correlation plot.
#' The plot can be saved as a PDF file. The user can specify the number of top correlated features
#' to be plotted, as well as customize plot settings such as the color palette and plot dimensions.
#'
#' Plot Top Correlated Features
#'
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot
#' @importFrom wesanderson wes_palette
#' @import here 
#' @import methods
#' @import stats
#' @param data A data frame or matrix containing the features to be analyzed.
#' @param top_n An integer specifying the number of top correlated features to be plotted (default is 15).
#' @param group_col A character string specifying the column name for grouping (optional, default is "group").
#' @param palette_name A string specifying the color palette for the plot .
#' @param save_plots A logical value indicating whether to save the plot as a PDF (default is TRUE).
#' @param plot_width A numeric value specifying the width of the saved plot (default is 5).
#' @param plot_height A numeric value specifying the height of the saved plot (default is 5).
#' @param save_dir A string specifying the directory where the plot will be saved (default is `here("StatObject", "correlation_info")`).
#' @param base_size A numeric value specifying the base font size for the plot (default is 14).
#'
#' @returns The function returns `invisible()`, and the plot is saved to the specified directory if `save_plots` is TRUE.
#'
#' @export
#'
#' @examples
#' # Plot the top 10 most correlated features from a data frame
#' plot_top_correlations(my_data_frame, top_n = 10)
#'
#' # Plot and save the top 20 most correlated features with custom plot dimensions
#' plot_top_correlations(my_data_frame, top_n = 20, plot_width = 7, plot_height = 6)
#'
plot_top_correlations <- function(data,
                                  top_n = 15,
                                  group_col = "group",
                                  palette_name = "Zissou1",
                                  save_plots = TRUE,
                                  plot_width = 5,
                                  plot_height = 5,
                                  save_dir = here("StatObject", "correlation_info"),
                                  base_size = 14) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("The input data must be of type data.frame or matrix")
  }

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  if (!is.null(group_col) && group_col %in% colnames(data)) {
    cat("Group column '", group_col, "' found. Excluding it from correlation analysis.\n")
    data <- data[, !colnames(data) %in% group_col, drop = FALSE]
  } else {
    cat("Group column '", group_col, "' not found or is NULL. Proceeding with all columns.\n")
  }

  if (!all(sapply(data, is.numeric))) {
    warning("Non-numeric columns detected. Attempting to convert to numeric.")
    data <- as.data.frame(sapply(data, as.numeric))
  }

  n_features <- ncol(data)
  cat("Number of features in data:", n_features, "\n")

  if (top_n > n_features) {
    warning(paste("top_n is greater than the number of features in the data:", n_features, "Setting top_n to", (n_features - 1)))
    top_n <- n_features - 1
  }

  cor_res <- rcorr(as.matrix(data), type = "pearson")
  cor_matrix <- cor_res$r

  avg_cor <- apply(abs(cor_matrix), 2, mean)
  top_cor_idx <- order(avg_cor, decreasing = TRUE)[1:top_n]

  top_cor_matrix <- cor_matrix[top_cor_idx, top_cor_idx]

  colors <- as.vector(wes_palette(palette_name, 100, type = "continuous"))

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    pdf(file.path(save_dir, "top_correlations_plot.pdf"), width = plot_width, height = plot_height)
  }

  corrplot(top_cor_matrix,
           order = 'AOE',
           col = colors,
           tl.col = "black",
           tl.cex = base_size / 14,
           plotCI = NULL)

  if (save_plots) {
    dev.off()
    cat("Plot saved to:", file.path(save_dir, "top_correlations_plot.pdf"), "\n")
  }

  invisible()
}

#' Correlation Analysis and Plot for Top Correlated Features
#'
#' This function extracts data from a given object or data frame, performs a correlation analysis
#' on the features, and then generates a plot for the top `n` most correlated features. The user can
#' specify a correlation threshold, the number of top correlated features to plot, and customize the plot settings.
#' If the input object is of class "Stat", the function updates the `corr.result` slot of the object with the plot.
#'
#'
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot
#' @importFrom wesanderson wes_palette
#' @import here 
#' @import methods
#' @import stats
#' @param object An object of class 'Stat' or a data frame containing the features to be analyzed.
#' @param correlation_threshold A numeric value specifying the correlation threshold for feature selection (default is 0.9).
#' @param top_n An integer specifying the number of top correlated features to be plotted (default is 15).
#' @param palette_name A string specifying the color palette for the plot.
#' @param data_type A string specifying the type of data to extract from the 'Stat' object. Either "clean" or "scale" (default is "scale").
#' @param save_plots A logical value indicating whether to save the plot as a PDF (default is TRUE).
#' @param plot_width A numeric value specifying the width of the saved plot (default is 5).
#' @param plot_height A numeric value specifying the height of the saved plot (default is 5).
#' @param save_dir A string specifying the directory where the plot will be saved (default is `here("StatObject", "correlation_info")`).
#' @param base_size A numeric value specifying the base font size for the plot (default is 14).
#'
#' @returns If the input is a 'Stat' object, the updated object with the correlation plot added to the `corr.result` slot.
#' Otherwise, the function returns the plot for the top correlated features.
#'
#' @export
#'
#' @examples
#' # Perform correlation analysis and plot top 10 correlated features
#' cor_top_correlations(my_stat_object, top_n = 10)
#'
#' # Perform analysis and save the plot of top 20 correlated features with custom plot dimensions
#' cor_top_correlations(my_stat_object, top_n = 20, plot_width = 7, plot_height = 6)
#'
cor_top_correlations <- function(object,
                                 correlation_threshold = 0.9,
                                 top_n = 15,
                                 palette_name = "Zissou1",
                                 data_type = "scale",
                                 save_plots = TRUE,
                                 plot_width = 5,
                                 plot_height = 5,
                                 save_dir = here("StatObject", "correlation_info"),
                                 base_size = 14) {
  if (inherits(object, "Stat")) {
    data <- if (data_type == "clean") ExtractCleanData(object) else ExtractScaleData(object)
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }

  plot_top <- plot_top_correlations(data,
                                    top_n = top_n,
                                    palette_name = palette_name,
                                    save_plots = save_plots,
                                    plot_width = plot_width,
                                    plot_height = plot_height,
                                    save_dir = save_dir,
                                    base_size = base_size)

  if (inherits(object, "Stat")) {
    object@corr.result[["top_cor_plot"]] <- plot_top
    cat("Updating 'Stat' object...\n")
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'corr.result' slot updated.\n")
    return(object)
  }
  return(plot_top)
}

#' Cross-Analysis of Two Variables with Scatter Plot and Correlation Test
#'
#' This function performs a cross-analysis of two numeric variables in a given data frame. It generates
#' a scatter plot of the two variables, fits a linear regression line, and displays the correlation coefficient
#' along with the p-value. The user can customize the plot appearance and save the plot as an image file.
#'
#' @import ggplot2
#' @import stats
#' @import  wesanderson
#' @importFrom ggprism theme_prism
#' @import here 
#' @importFrom rlang sym
#' @param data A data frame containing the variables to be analyzed.
#' @param vars A character vector of length 2 specifying the names of the two variables for the analysis.
#' @param group_col A character string specifying the name of the grouping column for coloring points in the plot (default is "group").
#' @param palette_name A string specifying the color palette to be used for the plot.
#' @param save_plots A logical value indicating whether to save the plot as a PDF file (default is TRUE).
#' @param plot_width A numeric value specifying the width of the saved plot in inches (default is 5).
#' @param plot_height A numeric value specifying the height of the saved plot in inches (default is 5).
#' @param save_dir A string specifying the directory where the plot will be saved (default is `here("StatObject", "correlation_info")`).
#'
#' @returns A scatter plot of the two specified variables, with the correlation coefficient and p-value annotated.
#'          If `save_plots` is TRUE, the plot is saved as a PDF file.
#'
#' @export
#'
#' @examples
#' # Perform cross-analysis of two variables 'var1' and 'var2'
#' cross_analysis(data = my_data, vars = c("var1", "var2"))
#'
#' # Perform cross-analysis with a custom plot size and save the plot
#' cross_analysis(data = my_data, vars = c("var1", "var2"), plot_width = 6, plot_height = 6)
#'

cross_analysis <- function(data,
                           vars,
                           group_col = "group",
                           palette_name = "Zissou1",
                           save_plots = TRUE,
                           plot_width = 5,
                           plot_height = 5,
                           save_dir = here("StatObject", "correlation_info")) {
  if (!is.data.frame(data)) {
    stop("Invalid input. Function expects a data frame.")
  }

  cat("Variables for analysis:", vars, "\n")

  if (length(vars) != 2) {
    stop("Invalid input. Function expects a character vector of length 2.")
  }

  var1 <- vars[1]
  var2 <- vars[2]

  if (is.numeric(data[[var1]]) && is.numeric(data[[var2]])) {
    plot_title <- paste("Scatter Plot of", var1, "vs", var2)

    cor_test <- cor.test(data[[var1]], data[[var2]])
    cor_coef <- round(cor_test$estimate, 2)

    if (cor_test$p.value < 0.00001) {
      p_value <- format(cor_test$p.value, scientific = TRUE, digits = 3)
    } else {
      p_value <- round(cor_test$p.value, 5)
    }

    gg <- ggplot(data, aes(x = !!sym(var1), y = !!sym(var2), color = !!sym(group_col))) +
      geom_point(size = 3, alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(subtitle = plot_title, y = var2, x = var1, title = "Scatter Plot") +
      scale_fill_manual(values = wes_palette(palette_name)) +
      theme_prism(base_size = 14) +
      theme(legend.position = "bottom") +
      annotate("text", x = Inf, y = Inf,
               label = paste("Correlation Coefficient (r) =", cor_coef, "\nP-value =", p_value),
               hjust = 1.1, vjust = 1.1, size = 5, color = "black",
               parse = FALSE)

    if (save_plots) {
      ggsave(filename = file.path(save_dir, paste0("scatter_plot_", var1, "_vs_", var2, ".pdf")),
             plot = gg,
             width = plot_width,
             height = plot_height,
             device = "pdf")
      cat("Plot saved to:", file.path(save_dir, paste0("scatter_plot_", var1, "_vs_", var2, ".pdf")), "\n")
    }

    print(gg)
  } else {
    stop("Unrecognized variable types. Unable to perform cross analysis.")
  }
}



#' Cross-Plot Analysis of Two Variables from a 'Stat' Object or Data Frame
#'
#' This function performs a cross-plot analysis of two variables by generating a scatter plot with a correlation test.
#' It handles input either as a 'Stat' object or a data frame, extracting the relevant data based on the specified `data_type`.
#' The plot will be saved as a PDF file if specified.
#' @import ggplot2
#' @import stats
#' @import  wesanderson
#' @importFrom ggprism theme_prism
#' @import here 
#' @importFrom rlang sym
#' @param object An object of class 'Stat' or a data frame containing the data for analysis.
#' @param vars A character vector of length 2 specifying the names of the two variables for the analysis.
#' @param palette_name A string specifying the color palette for the plot.
#' @param data_type A string indicating which type of data to extract from the 'Stat' object. Options are "clean" or "scale" (default is "scale").
#' @param group_col A string specifying the column name for the grouping variable (default is "group").
#' @param save_plots A logical value indicating whether to save the plot as a PDF file (default is TRUE).
#' @param plot_width A numeric value specifying the width of the saved plot in inches (default is 5).
#' @param plot_height A numeric value specifying the height of the saved plot in inches (default is 5).
#' @param save_dir A string specifying the directory where the plot will be saved (default is `here("StatObject", "correlation_info")`).
#'
#' @returns A scatter plot of the two specified variables, with the correlation coefficient and p-value annotated.
#'          If `save_plots` is TRUE, the plot is saved as a PDF file.
#'
#' @export
#'
#' @examples
#' # Perform cross-plot analysis on a 'Stat' object with specified variables
#' cross_plot_analysis(object = my_stat_object, vars = c("var1", "var2"))
#'
#' # Perform cross-plot analysis on a data frame with a custom plot size
#' cross_plot_analysis(object = my_data, vars = c("var1", "var2"), plot_width = 6, plot_height = 6)
#'
cross_plot_analysis <- function(object,
                                vars,
                                palette_name = "Zissou1",
                                data_type = "scale",
                                group_col = "group",
                                save_plots = TRUE,
                                plot_width = 5,
                                plot_height = 5,
                                save_dir = here("StatObject", "correlation_info")) {
  if (inherits(object, "Stat")) {
    data <- if (data_type == "clean") ExtractCleanData(object) else ExtractScaleData(object)
    group_col <- slot(object, "group_col")
    if (length(group_col) == 0 || !group_col %in% colnames(data)) {
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
  cross_analysis(data, vars, group_col = group_col, palette_name = palette_name, save_plots = save_plots, plot_width = plot_width, plot_height = plot_height, save_dir = save_dir)
}



