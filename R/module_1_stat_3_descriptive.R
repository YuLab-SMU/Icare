#' Compute Descriptive Statistics
#'
#' This function computes various descriptive statistics for the input dataset, including counts for categorical variables,
#' and summary statistics (mean, median, standard deviation, min, max) for numeric variables. It also checks for normality
#' of numeric variables and computes normality tests (Shapiro-Wilk or Anderson-Darling test), based on the number of unique
#' values and sample size.
#' @import stats
#' @importFrom nortest ad.test
#' @import methods
#' @param data A data frame containing the dataset to analyze.
#' @param count_feature A logical value indicating whether to compute counts for categorical variables (default is TRUE).
#' @param group_col The column name representing the grouping variable (default is "group").
#' @param max_unique_values The maximum number of unique values a variable can have to be considered categorical (default is 5).
#'
#' @returns A list containing the following components:
#' \item{Group_Counts}{The counts of each group if a `group_col` is provided.}
#' \item{Count_Results}{A list of counts for each categorical variable.}
#' \item{Num_Results}{A list of descriptive statistics for numeric variables, including separate statistics for normal and non-normal distributions.}
#' \item{Normality_Test}{A list of p-values and normality test results for numeric variables.}
#'
#' @export
#'
#' @examples
#' # Example 1: Compute descriptive statistics for a data frame
#' result <- compute_descriptive_stats(data = my_data, count_feature = TRUE, group_col = "group")
#'
#' # Example 2: Compute overall descriptive statistics without considering group
#' result <- compute_descriptive_stats(data = my_data, count_feature = FALSE)
#'
compute_descriptive_stats <- function(data,
                                      count_feature = TRUE,
                                      group_col = "group",
                                      max_unique_values = 5) {
  if (length(group_col) == 0) group_col <- NULL

  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)
  result <- list()

  if (count_feature && !is.null(variable_types) && !is.null(variable_types$categorical_vars)) {
    count_col <- variable_types$categorical_vars
    cat("Categorical variables identified:", count_col, "\n")

    if (group_col %in% colnames(data)) {
      data[[group_col]] <- as.factor(data[[group_col]])
      group_counts <- table(data[[group_col]])
      result$Group_Counts <- group_counts
      cat("Group counts computed for column:", group_col, "\n")
      num_cols <- setdiff(names(data), c(count_col, group_col))
    } else {
      num_cols <- setdiff(names(data), count_col)
    }

    cat("Counting categorical values...\n")
    count_results <- lapply(count_col, function(col) table(data[[col]]))
    names(count_results) <- count_col
    result$Count_Results <- count_results

    normal_vars <- character()
    non_normal_vars <- character()

    for (col in num_cols) {
      if (is.numeric(data[[col]])) {
        cat("Checking normality for numeric column:", col, "\n")
        if (length(unique(data[[col]])) <= 2) {
          non_normal_vars <- c(non_normal_vars, col)
        } else {
          if (sum(!is.na(data[[col]])) <= 5000) {
            is_normal <- shapiro.test(data[[col]])$p.value > 0.05
          } else {
            is_normal <- nortest::ad.test(data[[col]])$p.value > 0.05
          }
          if (is_normal) {
            normal_vars <- c(normal_vars, col)
          } else {
            non_normal_vars <- c(non_normal_vars, col)
          }
        }
      }
    }

    cat("Computing descriptive statistics for numeric columns...\n")
    all_stats <- apply(data[, num_cols, drop = FALSE], 2, function(x) {
      c(Mean = mean(x, na.rm = TRUE), Median = median(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), Max = max(x, na.rm = TRUE))
    })

    normal_stats <- if (length(normal_vars) > 0) {
      apply(data[, normal_vars, drop = FALSE], 2, function(x) {
        c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE))
      })
    } else {
      NULL
    }

    non_normal_stats <- if (length(non_normal_vars) > 0) {
      apply(data[, non_normal_vars, drop = FALSE], 2, function(x) {
        if (sum(!is.na(x)) <= 5000) {
          ad_test <- shapiro.test(x)
        } else {
          ad_test <- nortest::ad.test(x)
        }
        c(AD_p_value = ad_test$p.value, Median = median(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE))
      })
    } else {
      NULL
    }

    result$Num_Results <- list(All = all_stats, Normal = normal_stats, Non_Normal = non_normal_stats)
    result$Normality_Test <- lapply(data[, num_cols, drop = FALSE], function(x) {
      if (is.numeric(x)) {
        if (length(unique(x)) <= 2) {
          return(NULL)
        } else {
          if (sum(!is.na(x)) <= 5000) {
            test_result <- shapiro.test(x)
          } else {
            test_result <- nortest::ad.test(x)
          }
          return(list(p_value = test_result$p.value, is_normal = test_result$p.value > 0.05))
        }
      } else {
        return(NULL)
      }
    })

    cat("Descriptive statistics computation completed.\n")
  } else {
    cat("Computing overall descriptive statistics...\n")
    stats_compute <- apply(data, 2, function(x) c(Mean = mean(x, na.rm = TRUE), Median = median(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), Max = max(x, na.rm = TRUE)))
    result$Stats <- stats_compute
    cat("Overall statistics computed.\n")
  }

  return(result)
}




#' Compute Descriptive Statistics for 'Stat' Object or Data Frame
#'
#' This function computes descriptive statistics for numeric and categorical variables in the input data. It can handle both
#' `Stat` objects and standard data frames. The function computes counts for categorical variables, and summary statistics
#' (mean, median, standard deviation, min, max) for numeric variables, along with normality tests for numeric variables.
#'
#' @import stats
#' @importFrom nortest ad.test
#' @import methods
#' @param object An object of class 'Stat' or a data frame containing the dataset to analyze.
#' @param count_feature A logical value indicating whether to compute counts for categorical variables (default is TRUE).
#' @param group_col The column name representing the grouping variable (default is "group").
#' @param max_unique_values The maximum number of unique values a variable can have to be considered categorical (default is 5).
#'
#' @returns If the input is a 'Stat' object, the function returns the updated 'Stat' object with a new slot `compute.descriptive`
#'          containing the results. If the input is a data frame, it returns the descriptive statistics as a list.
#'
#' @export
#'
#' @examples
#' # Example 1: Compute descriptive statistics for a 'Stat' object
#' updated_stat <- stat_compute_descriptive(object = my_stat, count_feature = TRUE, group_col = "group")
#'
#' # Example 2: Compute descriptive statistics for a data frame
#' result <- stat_compute_descriptive(object = my_data, count_feature = TRUE)
stat_compute_descriptive <- function(
    object,
    count_feature = TRUE,
    group_col = "group",
    max_unique_values = 5
) {

  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    } 
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

  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)

  if (length(variable_types$numeric_vars) == 0 && length(variable_types$categorical_vars) == 0) {
    stop("No valid variables found after variable type diagnosis")
  }

  result <- list()

  result <- compute_descriptive_stats(data = data, count_feature = count_feature)

  if (inherits(object, "Stat")) {
    object@compute.descriptive <- result
    cat("Updating 'Stat' object...\n")
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'compute.descriptive' slot updated.\n")

    return(object)
  }

  return(result)
}




#' Plot Categorical Descriptive Statistics for 'Stat' Object or Data Frame
#'
#' This function generates bar plots displaying the distribution of categorical variables and their percentages across different groups.
#' It can handle both `Stat` objects and standard data frames. The function also supports saving plots to a specified directory.
#' The plots are created using `ggplot2` and the colors are customized using a specified color palette.
#'
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate
#' @import wesanderson
#' @importFrom ggprism theme_prism
#' @import here
#' @importFrom scales percent
#' @importFrom gridExtra grid.arrange
#' @import methods
#' @importFrom rlang sym
#' @param object An object of class 'Stat' or a data frame containing the dataset to analyze.
#' @param group_col The column name representing the grouping variable (default is "group").
#' @param palette_name The name of the color palette to use (default is "Royal1").
#' @param save_plots A logical value indicating whether to save the plots (default is TRUE).
#' @param save_dir The directory where the plots will be saved (default is "StatObject/categorical_descriptive").
#' @param sub_var A vector of column names to subset the variables to plot (default is NULL, meaning all categorical variables will be plotted).
#' @param plot_width The width of the saved plot (default is 5).
#' @param plot_height The height of the saved plot (default is 5).
#' @param base_size The base font size for the plots (default is 14).
#'
#' @returns If the input is a 'Stat' object, the function returns the updated 'Stat' object with a new slot `compute.descriptive$count_plots`
#'          containing the generated plots. If the input is a data frame, it returns a list of the generated plots.
#'
#' @export
#'
#' @examples
#' # Example 1: Plot categorical descriptive statistics for a 'Stat' object
#' updated_stat <- plot_categorical_descriptive(object = my_stat, group_col = "group")
#'
#' # Example 2: Plot categorical descriptive statistics for a data frame
#' plot_results <- plot_categorical_descriptive(object = my_data, group_col = "group")
plot_categorical_descriptive <- function(
    object,
    group_col = "group",
    palette_name = "Royal1",
    save_plots = TRUE,
    save_dir = here("StatObject", "categorical_descriptive"),
    sub_var = NULL,
    plot_width = 5,
    plot_height = 5,
    base_size = 14
) {
  tryCatch({
    if (inherits(object, "Stat")) {
      data <- slot(object, "clean.data")
      if (is.null(data) || nrow(data) == 0) {
        data <- slot(object, "raw.data")
      } 
      group_col <- slot(object, "group_col")

      if (length(group_col) == 0 || !group_col %in% colnames(data)) {
        group_col <- NULL
      }
    } else if (is.data.frame(object)) {
      data <- object
    } else {
      stop("Input must be an object of class 'Stat' or a data frame.\n")
    }

    if (is.null(data) || nrow(data) == 0) {
      stop("No valid data found in the input.\n")
    }

    if (inherits(object, "Stat")) {
      stats <- slot(object, "compute.descriptive")
    } else {
      stats <- compute_descriptive_stats(data, count_feature = TRUE, group_col = group_col)
    }

    if (is.null(stats$Count_Results)) {
      cat("No valid count results to display. Please check the data and parameters.\n")
      return(object)
    }

    count_plots <- list()

    for (col in names(stats$Count_Results)) {
      cat(paste("Processing column:", col, "\n"))

      if (!is.null(sub_var) && !(col %in% sub_var)) {
        cat(paste("Skipping column:", col, "not in sub_var.\n"))
        next
      }

      if (is.null(group_col)) {
        group_col <- "temp_group"
        data[[group_col]] <- "All"
      }

      data[[group_col]] <- as.factor(data[[group_col]])

      plot_data <- data %>%
        group_by(!!sym(col), !!sym(group_col)) %>%
        summarise(n = n(), .groups = 'drop') %>%
        mutate(perc = n / sum(n) * 100)

      if (nrow(plot_data) == 0) {
        cat(paste("No data available for column:", col, "\n"))
        next
      }

      p <- ggplot(plot_data, aes_string(x = col, y = "perc", fill = group_col)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = scales::percent(perc / 100)),
                  position = position_dodge(width = 0.9), size = 4, vjust = -0.5, hjust = 0.5) +
        labs(title = paste("Bar plot of", col, "with percentage"), y = "Percentage") +
        theme_prism(base_size = base_size)

      p <- p + scale_fill_manual(values = wes_palette(palette_name))

      count_plots[[col]] <- p

      if (save_plots) {
        ggsave(here(save_dir, paste0(col, "_plot.pdf")),
               plot = p,
               width = plot_width,
               height = plot_height,
               device = "pdf")
        cat(paste("Saved plot for column:", col, "\n"))
      }
    }

    plot_count <- length(count_plots)
    cat(paste("A total of", plot_count, "plots were generated.\n"))

    if (plot_count > 0) {
      grid.arrange(grobs = count_plots, nrow = 1)
      print(count_plots[[1]])
    } else {
      cat("No valid count plots to display.\n")
    }

    if (save_plots) {
      cat("Saved plots to:", save_dir, "\n")
    }

    if (inherits(object, "Stat")) {
      object@compute.descriptive[["count_plots"]] <- count_plots
      cat("Updating 'Stat' object...\n")
      cat("The 'Stat' object has been updated with the following slots:\n")
      cat("- 'compute.descriptive' slot updated.\n")
      return(object)
    }
    return(count_plots)
  }, error = function(e) {
    cat("An error occurred. Function terminated silently.\n")
    return(object)
  })
}

#' Violin Plots for Numeric Variables
#'
#' This function generates violin plots for numeric variables in the dataset,
#' with boxplots and jittered points for enhanced visualization. The plots
#' are grouped by a specified column and saved as PDF files, if desired.
#' It handles large numbers of variables by splitting them into multiple
#' plots, each containing a set of variables (defined by the `vars_per_plot`
#' argument). The function also allows for customization of plot appearance
#' and save locations.
#'
#'
#' @import ggplot2
#' @import wesanderson
#' @importFrom reshape2 melt
#' @importFrom ggpubr stat_compare_means
#' @import here
#' @import stats
#' @param data A data frame containing the data to plot.
#' @param vars_per_plot Integer. The number of variables (columns) to include
#'   in each individual plot. Default is 5.
#' @param save_dir String. The directory to save the plots in. Default is
#'   the directory `"StatObject"`.
#' @param palette_name String. The name of the color palette to use. The
#'   default is `"Royal1"`, which uses a predefined color scheme.
#' @param group_col String. The name of the column in the data used to group
#'   the samples. Default is `"group"`.
#' @param max_unique_values Integer. The maximum number of unique values
#'   allowed for categorical variables. Defaults to 5. Used for filtering
#'   variables when diagnosing their types.
#' @param sub_var Character vector. A subset of variable names (columns)
#'   to include in the plot. If `NULL`, all numeric variables will be included.
#' @param save_plots Logical. If `TRUE`, the plots will be saved as PDF
#'   files in the specified `save_dir`. Default is `TRUE`.
#' @param plot_width Numeric. The width of each saved plot (in inches).
#'   Default is 5.
#' @param plot_height Numeric. The height of each saved plot (in inches).
#'   Default is 5.
#' @param base_size Numeric. The base font size for the plot. Default is 14.
#'
#' @returns A list of ggplot objects containing the generated violin plots.
#'   If `save_plots` is `TRUE`, the plots are also saved as PDF files in
#'   the specified directory.
#' @export
#'
#' @examples
#' violin_plots(data = my_data,
#'              vars_per_plot = 4,
#'              save_dir = here::here("StatObject"),
#'              palette_name = "Zissou1",
#'              group_col = "group",
#'              save_plots = TRUE)
violin_plots <- function(data,
                         vars_per_plot = 1,
                         save_dir = here::here("StatObject"),
                         palette_name = "Zissou1",
                         group_col = "group",
                         max_unique_values = 5,
                         sub_var = NULL,
                         save_plots = TRUE,
                         plot_width = 5,
                         plot_height = 5,
                         base_size = 14,
                         stat_method = "wilcox.test",
                         paired_comparison = TRUE) {

  if (!is.null(group_col) && !group_col %in% names(data)) {
    stop(paste("Column", group_col, "not found in data"))
  }

  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)
  num_cols <- variable_types$numeric_vars

  if (!is.null(sub_var)) {
    num_cols <- intersect(num_cols, sub_var)
  }

  if (length(num_cols) == 0) {
    stop("No numeric variables found after filtering with sub_var")
  }

  if (is.null(group_col)) {
    melted_data <- melt(data, measure.vars = num_cols)
    melted_data$group <- "All"
  } else {
    melted_data <- melt(data, id.vars = group_col, measure.vars = num_cols)
    melted_data$group <- as.character(melted_data[[group_col]])
  }

  num_vars <- length(num_cols)
  num_groups <- ceiling(num_vars / vars_per_plot)
  var_groups <- split(num_cols, rep(1:num_groups, each = vars_per_plot, length.out = num_vars))

  plot_list <- list()

  for (i in seq_along(var_groups)) {
    group_vars <- var_groups[[i]]
    group_data <- melted_data[melted_data$variable %in% group_vars, ]

    if (nrow(group_data) > 0) {
      if (is.null(group_col)) {
        pal <- wes_palette(palette_name, n = 1, type = "continuous")
      } else {
        pal <- wes_palette(palette_name, n = length(unique(group_data$group)), type = "continuous")
      }

      p <- ggplot(group_data, aes(x = as.factor(group), y = value)) +
        geom_violin(aes(fill = group), scale = "area", alpha = 0.5) +
        geom_boxplot(width = 0.1, size = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3, color = "black", size = 0.5) +
        scale_fill_manual(values = pal) +
        facet_wrap(~variable, scales = "free_y", ncol = 1) +
        theme_minimal(base_size = base_size) +
        theme(legend.position = "bottom") +
        labs(title = paste("Violin Plots - Part", i), x = "Group", y = "Value")

      if (!is.null(group_col) && paired_comparison) {
        groups <- unique(group_data$group)
        if (length(groups) >= 3) {
          comparisons <- list(
            c(groups[1], groups[3]),
            c(groups[2], groups[3]),
            c(groups[1], groups[2])
          )
        } else {
          comparisons <- combn(groups, 2, simplify = FALSE)
        }

        if (length(comparisons) > 0) {
          p <- p + stat_compare_means(
            method = stat_method,
            comparisons = comparisons,
            label = "p.format",
            size = 4,
            tip.length = 0.02,
            step.increase = 0.08
          )
        }
      }

      plot_list[[i]] <- p

      if (save_plots) {
        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }
        ggsave(
          filename = paste0(save_dir, "/violin_plot_part_", i, ".pdf"),
          plot = p,
          width = plot_width,
          height = plot_height,
          device = "pdf"
        )
      }
    }
  }

  return(plot_list)
}

#' Density Ridge Plots for Numeric Variables
#'
#' This function generates density ridge plots for numeric variables in the dataset,
#' with optional grouping by a specified column. It handles large numbers of variables
#' by splitting them into multiple plots, each containing a set of variables (defined
#' by the `vars_per_plot` argument). The function also allows for customization of plot
#' appearance and save locations.
#'
#' @import ggplot2
#' @import wesanderson
#' @importFrom reshape2 melt
#' @importFrom ggpubr stat_compare_means
#' @import here
#' @import stats
#' @importFrom ggridges geom_density_ridges_gradient
#' @param data A data frame containing the data to plot.
#' @param vars_per_plot Integer. The number of variables (columns) to include
#'   in each individual plot. Default is 5.
#' @param save_dir String. The directory to save the plots in. Default is
#'   the directory `"StatObject"`.
#' @param palette_name String. The name of the color palette to use. The
#'   default is `"Royal1"`, which uses a predefined color scheme.
#' @param group_col String. The name of the column in the data used to group
#'   the samples. Default is `"group"`.
#' @param max_unique_values Integer. The maximum number of unique values
#'   allowed for categorical variables. Defaults to 5. Used for filtering
#'   variables when diagnosing their types.
#' @param sub_var Character vector. A subset of variable names (columns)
#'   to include in the plot. If `NULL`, all numeric variables will be included.
#' @param save_plots Logical. If `TRUE`, the plots will be saved as PDF
#'   files in the specified `save_dir`. Default is `TRUE`.
#' @param plot_width Numeric. The width of each saved plot (in inches).
#'   Default is 5.
#' @param plot_height Numeric. The height of each saved plot (in inches).
#'   Default is 5.
#' @param base_size Numeric. The base font size for the plot. Default is 14.
#'
#' @returns A list of ggplot objects containing the generated density ridge plots.
#'   If `save_plots` is `TRUE`, the plots are also saved as PDF files in
#'   the specified directory.
#' @export
#'
#' @examples
#' density_ridge_plots(data = my_data,
#'                     vars_per_plot = 4,
#'                     save_dir = here("StatObject"),
#'                     palette_name = "Zissou1",
#'                     group_col = "group",
#'                     save_plots = TRUE)
#'
density_ridge_plots <- function(data,
                                vars_per_plot = 1,
                                save_dir = here("StatObject"),
                                palette_name = "Zissou1",
                                group_col = "group",
                                max_unique_values = 5,
                                sub_var = NULL,
                                save_plots = TRUE,
                                plot_width = 5,
                                plot_height = 5,
                                base_size = 14) {

  if (length(group_col) == 0 || is.null(group_col) || !group_col %in% names(data)) {
    group_col <- NULL
  }

  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)
  num_cols <- variable_types$numeric_vars

  if (!is.null(sub_var)) {
    sub_var <- intersect(sub_var, num_cols)

    if (length(sub_var) == 0) {
      stop("No valid variables in sub_var for plotting.")
    }

    num_cols <- sub_var
  }

  if (length(num_cols) == 0) {
    stop("No numeric variables found after filtering with sub_var\n")
  }

  melted_data <- melt(data, id.vars = group_col, measure.vars = num_cols)

  if (!is.null(group_col)) {
    melted_data$group <- melted_data[[group_col]]
  } else {
    melted_data$group <- "Default Group"
  }

  num_vars <- length(num_cols)
  num_groups <- ceiling(num_vars / vars_per_plot)
  var_groups <- split(num_cols, rep(1:num_groups, each = vars_per_plot, length.out = num_vars))

  pal <- wes_palette(palette_name, 100, type = "continuous")
  plot_list <- list()

  for (i in seq_along(var_groups)) {
    group_vars <- var_groups[[i]]
    group_data <- melted_data[melted_data$variable %in% group_vars, ]
    group_data$group <- as.factor(group_data$group)
    if (nrow(group_data) > 0) {
      p <- ggplot(group_data, aes(x = value, y = group, fill = after_stat(density))) +
        geom_density_ridges_gradient(scale = 3, rel_min_height = 0.00) +
        scale_fill_gradientn(colours = pal) +
        facet_wrap(~variable, scales = "free_x", ncol = 1) +
        theme_prism(base_size = base_size) +
        labs(title = paste("Density Ridge Plots  - Part", i),
             x = "Value",
             y = ifelse(!is.null(group_col), group_col, "Group"))

      plot_list[[i]] <- p

      if (save_plots) {
        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }

        ggsave(paste0(save_dir, "/density_ridge_plot_part_", i, ".pdf"),
               plot = p,
               width = plot_width,
               height = plot_height,
               device = "pdf")
        cat("Saved plots to: ", save_dir, "\n")
      }
    }
  }

  return(plot_list)
}

#' Numeric Descriptive Plots (Violin or Ridge Density)
#'
#' This function generates either violin plots or ridge density plots for numeric variables
#' in the given dataset (either a `Stat` object or a data frame). The function allows
#' customizing the number of variables per plot, palette style, and grouping columns, and
#' saves the plots if required. It provides a flexible approach to visualizing the distribution
#' of numeric data and comparing different groups.
#' @import ggplot2
#' @import wesanderson
#' @importFrom reshape2 melt
#' @importFrom ggpubr stat_compare_means
#' @import here
#' @import stats
#' @param object An object of class `Stat` or a data frame containing numeric data.
#'   If the object is of class `Stat`, the clean data is extracted from the `Stat` object.
#' @param vars_per_plot Integer. The number of variables (columns) to include in each
#'   individual plot. Default is 5.
#' @param save_dir String. The directory where the plots will be saved. Default is
#'   `"StatObject/numeric_descriptive"`.
#' @param palette_name String. The name of the color palette to use. Default is
#'   `"Zissou1"`.
#' @param group_col String. The name of the column in the data to group the samples by.
#'   Default is `"group"`. If not provided, no grouping is performed.
#' @param max_unique_values Integer. The maximum number of unique values allowed for
#'   categorical variables. Default is 5. Used for filtering variables when diagnosing
#'   their types.
#' @param plot_type String. The type of plot to generate. Options are `"violin"` or
#'   `"ridge"`. Default is `"violin"`.
#' @param save_plots Logical. If `TRUE`, the plots will be saved as PDF files in
#'   the specified `save_dir`. Default is `TRUE`.
#' @param plot_width Numeric. The width of each saved plot (in inches). Default is 5.
#' @param plot_height Numeric. The height of each saved plot (in inches). Default is 5.
#' @param base_size Numeric. The base font size for the plot. Default is 14.
#' @param sub_var Character vector. A subset of variable names (columns) to include in
#'   the plot. If `NULL`, all numeric variables are included.
#'
#' @returns The input object (either a `Stat` object or a data frame) with the
#'   updated plots added to the appropriate slot. If the input is a `Stat` object,
#'   the `compute.descriptive` slot is updated with the generated plots.
#' @export
#'
#' @examples
#' plot_numeric_descriptive(object = my_stat_object,
#'                           vars_per_plot = 4,
#'                           save_dir = here("StatObject", "numeric_descriptive"),
#'                           palette_name = "Zissou1",
#'                           group_col = "group",
#'                           plot_type = "violin",
#'                           save_plots = TRUE)
#'
plot_numeric_descriptive <- function(
    object,
    vars_per_plot = 1,
    save_dir = here("StatObject", "numeric_descriptive"),
    palette_name = "Zissou1",
    group_col = "group",
    max_unique_values = 5,
    plot_type = "violin",
    save_plots = TRUE,
    plot_width = 5,
    plot_height = 5,
    base_size = 14,
    sub_var = NULL) {


  if (inherits(object, "Stat")) {
    cat("Input is an object of class 'Stat'. Extracting data...\n")
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    } 
    group_col <- object@group_col
    if (length(group_col) == 0) {
      group_col <- NULL
      cat("No group_col found in Stat object. Using NULL.\n")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using it directly.\n")
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame\n")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input\n")
  }


  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)


  if (length(variable_types$numeric_vars) == 0) {
    stop("No valid numeric variables found after diagnosis\n")
  }

  if (!is.null(sub_var)) {
    cat("Filtering numeric variables based on sub_var:", sub_var, "\n")
    numeric_vars <- intersect(variable_types$numeric_vars, sub_var)
  } else {
    numeric_vars <- variable_types$numeric_vars
  }

  if (length(numeric_vars) == 0) {
    stop("No valid numeric variables found after filtering with sub_var\n")
  }


  plots_list <- list()

  if (plot_type == "violin") {
    cat("Generating violin plots...\n")
    plots_list <- violin_plots(data,
                               vars_per_plot = vars_per_plot,
                               save_dir = save_dir,
                               palette_name = palette_name,
                               group_col = group_col,
                               max_unique_values = max_unique_values,
                               save_plots = save_plots,
                               plot_width = plot_width,
                               plot_height = plot_height,
                               base_size = base_size,
                               sub_var = numeric_vars)
  } else if (plot_type == "ridge") {
    cat("Generating ridge density plots...\n")
    plots_list <- density_ridge_plots(data,
                                      vars_per_plot = vars_per_plot,
                                      save_dir = save_dir,
                                      palette_name = palette_name,
                                      group_col = group_col,
                                      max_unique_values = max_unique_values,
                                      save_plots = save_plots,
                                      plot_width = plot_width,
                                      plot_height = plot_height,
                                      base_size = base_size,
                                      sub_var = numeric_vars)
  } else {
    stop("Invalid plot type. Choose either 'violin' or 'ridge'.\n")
  }

  if (length(plots_list) > 0) {
    print(plots_list[[1]])
  }

  total_plots <- length(plots_list)
  cat("A total of", total_plots, "plots were generated.\n")

  if (inherits(object, "Stat")) {
    object@compute.descriptive[[paste0(plot_type, "_plots")]] <- plots_list
    cat("Updating 'Stat' object...\n")
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'compute.descriptive' slot updated.\n")
  }


  return(object)
}


#' Convert Variables to Numeric or Factor Based on Their Types
#'
#' This function converts the variables (columns) in the input data frame to either
#' numeric or factor types based on the information provided in `variable_types`.
#' Numeric variables are converted to numeric data type, and non-numeric variables
#' are converted to factors. The function uses `variable_types` to determine which
#' variables should be treated as numeric and which should be treated as factors.
#'
#' @param data A data frame containing the variables to be converted.
#' @param variable_types A list containing the variable types. This should have a
#'   component `numeric_vars`, which is a character vector of variable names that
#'   should be converted to numeric.
#'
#' @returns A data frame with the variables converted to the appropriate types
#'   (numeric or factor).
#' @export
#'
#' @examples
#' # Example of usage
#' data <- data.frame(a = c(1, 2, 3), b = c("low", "medium", "high"))
#' variable_types <- list(numeric_vars = c("a"))
#' converted_data <- convert_variables(data, variable_types)
#'
convert_variables <- function(data, 
                              variable_types,
                              save_dir = here::here("StatObject","Data"),
                              save_data = F,
                              csv_filename = "clean_data.csv") {
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
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  
  return(data)
}

#' Convert Variables in a 'Stat' Object or Data Frame
#'
#' This function converts variables in a given object (either of class 'Stat' or a
#' data frame) to numeric or factor types based on the information provided by
#' `diagnose_variable_type`. If the input is a 'Stat' object, the conversion will
#' update its `clean.data` slot. If the input is a data frame, it will directly
#' return the converted data frame.
#'
#' @param object An object of class 'Stat' or a data frame. If a 'Stat' object is
#'   provided, the function will update its `clean.data` slot.
#' @param group_col A string representing the column name that groups the data.
#'   Default is "group". This column is used to determine the type of variables.
#' @param max_unique_values The maximum number of unique values allowed for a variable
#'   to be considered as numeric. Default is 5.
#'
#' @returns If the input is a 'Stat' object, it returns the updated 'Stat' object.
#'   If the input is a data frame, it returns the converted data frame.
#' @export
#'
#' @examples
#' # Example of usage
#' data <- data.frame(a = c(1, 2, 3), b = c("low", "medium", "high"))
#' stat_object <- Stat$new(clean.data = data, group_col = "group")
#' updated_stat_object <- stat_convert_variables(stat_object)
#'
stat_convert_variables <- function(object,
                                   group_col = "group",
                                   max_unique_values = 5,
                                   save_dir = here::here("StatObject","Data"),
                                   save_data = TRUE,
                                   csv_filename = "clean_data.csv") {
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

  variable_types <- diagnose_variable_type(data, group_col = group_col, max_unique_values = max_unique_values)

  row_names <- rownames(data)
  data <- convert_variables(data, 
                            variable_types,
                            save_dir =save_dir,
                            save_data = save_data,
                            csv_filename = csv_filename)
  rownames(data) <- row_names

  if (inherits(object, "Stat")) {
    object@clean.data <- data
    cat("Updating 'Stat' object...\n")
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'clean.data' slot updated.\n")

    return(object)
  }

  return(data)
}





#' One-Hot Encode Categorical Variables in a Data Frame
#'
#' This function performs one-hot encoding on categorical variables with limited unique values,
#' as determined by the `diagnose_variable_type` function. It skips encoding the `group_col`
#' and moves it to the last column of the resulting data frame if provided.
#'
#' @param data A data frame to be processed.
#' @param group_col The column name for grouping (e.g., "group"). This column will not be encoded
#'        and will be moved to the end of the result. Default is `"group"`.
#' @param max_unique_values The maximum number of unique values a variable can have to be
#'        considered for one-hot encoding. Default is 5.
#'
#' @returns A new data frame with one-hot encoded variables.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' one_hot_encode(iris, group_col = "Species", max_unique_values = 5)
#' }
one_hot_encode <- function(data, 
                           
                           group_col = "group",
                           max_unique_values = 5,
                           save_dir = here::here("StatObject","Data"),
                           save_data = TRUE,
                           csv_filename = "clean_data.csv"
) {
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
    cat("Encoding variable:", var, "with unique values:", unique_values, "\n")
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
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(encoded_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  return(encoded_data)
}

#' Apply One-Hot Encoding to a Stat Object or Data Frame
#'
#' This function applies one-hot encoding to categorical variables within a `Stat` object
#' or a regular `data.frame`, using `one_hot_encode()`. If a `Stat` object is provided,
#' the encoded result will be updated in its `clean.data` slot.
#'
#' @import methods
#' @param object An object of class `Stat` or a `data.frame`.
#' @param method Reserved for future use. Currently not used. Default is 1.
#' @param group_col A string specifying the name of the grouping column. This column will be
#'        excluded from one-hot encoding and moved to the end of the result. Default is `"group"`.
#' @param max_unique_values The maximum number of unique values a variable can have to be
#'        considered for one-hot encoding. Default is 5.
#'
#' @returns A `Stat` object with updated one-hot encoded `clean.data` slot, or a
#'          one-hot encoded `data.frame`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' stat_onehot_encode(my_data, group_col = "status")
#' stat_onehot_encode(my_stat_object)
#' }
#'
stat_onehot_encode <- function(object, 
                               method = 1, 
                               group_col = "group", 
                               max_unique_values = 5,
                               save_dir = here::here("StatObject","Data"),
                               save_data = TRUE,
                               csv_filename = "clean_data.csv") {
  cat("Input object class:", class(object), "\n")

  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    group_col <- slot(object, "group_col")

    if (length(group_col) == 0 || is.null(group_col)) {
      cat("Group column is not valid, setting to NULL.\n")
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

  cat("Starting one-hot encoding on data...\n")

  onehot_data <- one_hot_encode(data, 
                                
                                group_col = group_col, 
                                max_unique_values = max_unique_values,
                                save_dir =save_dir,
                                save_data = save_data,
                                csv_filename = csv_filename)

  if (inherits(object, "Stat")) {
    if (!is.null(slotNames(object)) && "clean.data" %in% slotNames(object)) {
      object@clean.data <- onehot_data
      cat("Updating 'Stat' object...\n")
      cat("The 'Stat' object has been updated with the following slots:\n")
      cat("- 'clean.data' slot updated.\n")
    } else {
      stop("The 'Stat' object does not have a 'clean.data' slot.")
    }
    return(object)
  }

  cat("One-hot encoding complete, returning encoded data frame.\n")
  return(onehot_data)
}
