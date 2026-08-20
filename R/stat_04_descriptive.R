#' Compute Descriptive Statistics
#'
#' This function computes various descriptive statistics for the input dataset, including counts for categorical variables,
#' and summary statistics (mean, median, standard deviation, min, max) for numeric variables. It also checks for normality
#' of numeric variables and computes normality tests (Shapiro-Wilk or Anderson-Darling test), based on the number of unique
#' values and sample size.
#' @importFrom stats sd median quantile var shapiro.test p.adjust
#' @importFrom nortest ad.test
#' @import methods
#' @param data A data frame containing the dataset to analyze.
#' @param count_feature A logical value indicating whether to compute counts for categorical variables (default is TRUE).
#' @param group_col The column name representing the grouping variable (default is "group").
#' @param max_unique_values The maximum number of unique values a variable can have to be considered categorical (default is 5).
#'
#' @returns A list containing the following components:
#' \describe{
#' \item{Group_Counts}{The counts of each group if a `group_col` is provided.}
#' \item{Count_Results}{A list of counts for each categorical variable.}
#' \item{Num_Results}{A list of descriptive statistics for numeric variables, including separate statistics for normal and non-normal distributions.}
#' \item{Normality_Test}{A list of p-values and normality test results for numeric variables.}
#' }
#'
#' @export
#'
#' @examples
#' data(iris)
#' compute_descriptive_stats(iris, group_col = "Species")
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
#' @importFrom stats sd median quantile var shapiro.test p.adjust
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
#' \dontrun{
#' data("mtcars")
#' mtcars$cyl <- factor(mtcars$cyl)
#' mtcars$vs  <- factor(mtcars$vs, labels = c("V-shaped", "Straight"))
#' # Create a Stat object with "vs" as the grouping variable
#' stat_obj <- CreateStatObject(raw.data = mtcars, group_col = "vs")
#' #Compute descriptive statistics for a 'Stat' object
#' updated_stat <- stat_compute_descriptive(object = stat_obj, 
#' count_feature = TRUE, group_col = "group")
#' }
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
#' This function generates bar plots displaying the distribution of categorical
#' variables and their percentages across different groups. It handles both
#' `Stat` objects and data frames. Plots can optionally be saved.
#'
#' @param object A Stat object or a data frame.
#' @param group_col Grouping column name (default "group").
#' @param palette_name Color palette name (default "Royal1").
#' @param save_plots Logical; whether to save the plots (default TRUE).
#' @param save_dir Directory to save plots (if NULL and save_plots=TRUE,
#'   a default figure directory is used).
#' @param sub_var Character vector of variable names to plot (NULL = all).
#' @param plot_width Width of saved plot.
#' @param plot_height Height of saved plot.
#' @param base_size Base font size.
#'
#' @returns If object is a Stat object, the updated Stat object with
#'   `compute.descriptive$count_plots` slot. Otherwise, a list of plots.
#' @export
#' @examples
#' if (interactive()) {
#'   stat <- CreateStatObject(raw.data = mtcars, group_col = "cyl")
#'   plot_categorical_descriptive(stat, save_plots = FALSE)
#' }
plot_categorical_descriptive <- function(
    object,
    group_col = "group",
    palette_name = "Royal1",
    save_plots = TRUE,
    save_dir = NULL,
    sub_var = NULL,
    plot_width = 5,
    plot_height = 5,
    base_size = 14
) {
  tryCatch({
    # ---- 1. Extract data and group column ----
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
    
    # ---- 2. Get descriptive statistics ----
    if (inherits(object, "Stat")) {
      stats <- slot(object, "compute.descriptive")
    } else {
      stats <- compute_descriptive_stats(data, count_feature = TRUE, group_col = group_col)
    }
    
    if (is.null(stats$Count_Results)) {
      cat("No valid count results to display.\n")
      return(object)
    }
    
    count_plots <- list()
    
    # ---- 3. Loop through categorical variables ----
    for (col in names(stats$Count_Results)) {
      cat("Processing column:", col, "\n")
      
      # Skip if not in sub_var
      if (!is.null(sub_var) && !(col %in% sub_var)) {
        cat("Skipping column:", col, "not in sub_var.\n")
        next
      }
      
      # Skip if the categorical column is the same as the grouping column
      if (!is.null(group_col) && col == group_col) {
        cat("Skipping column:", col, "(same as group column)\n")
        next
      }
      
      # Create a dummy group if grouping column is absent
      if (is.null(group_col)) {
        group_col <- "temp_group"
        data[[group_col]] <- "All"
      }
      
      data[[group_col]] <- as.factor(data[[group_col]])
      
      # Compute percentages within groups
      plot_data <- data %>%
        dplyr::group_by(!!rlang::sym(col), !!rlang::sym(group_col)) %>%
        dplyr::summarise(n = dplyr::n(), .groups = 'drop') %>%
        dplyr::group_by(!!rlang::sym(group_col)) %>%
        dplyr::mutate(perc = n / sum(n) * 100) %>%
        dplyr::ungroup()
      
      if (nrow(plot_data) == 0) {
        cat("No data available for column:", col, "\n")
        next
      }
      
      # Create bar plot
      p <- ggplot2::ggplot(plot_data,
                           ggplot2::aes_string(x = col, y = "perc", fill = group_col)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::geom_text(ggplot2::aes(label = scales::percent(perc / 100)),
                           position = ggplot2::position_dodge(width = 0.9),
                           size = 4, vjust = -0.5, hjust = 0.5) +
        ggplot2::labs(title = paste("Bar plot of", col, "with percentage"), y = "Percentage") +
        ggprism::theme_prism(base_size = base_size) +
        ggplot2::scale_fill_manual(values = wesanderson::wes_palette(palette_name))
      
      count_plots[[col]] <- p
      
      # Save plot if requested
      if (save_plots) {
        if (is.null(save_dir)) {
          # Use a default figure directory (from viz_functions helpers)
          save_dir <- .get_viz_output_dir("Stat")
        }
        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }
        ggplot2::ggsave(filename = file.path(save_dir, paste0(col, "_plot.pdf")),
                        plot = p, width = plot_width, height = plot_height, device = "pdf")
        cat("Saved plot for column:", col, "\n")
      }
    }
    
    plot_count <- length(count_plots)
    cat("A total of", plot_count, "plots were generated.\n")
    
    if (plot_count > 0) {
      # Use gridExtra::grid.arrange to avoid dependency on global library
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(grobs = count_plots, nrow = 1)
      } else {
        # fallback: print the first plot
        print(count_plots[[1]])
      }
      print(count_plots[[1]])
    } else {
      cat("No valid count plots to display.\n")
    }
    
    if (save_plots) {
      cat("Saved plots to:", save_dir, "\n")
    }
    
    # ---- 4. Update Stat object if applicable ----
    if (inherits(object, "Stat")) {
      object@compute.descriptive[["count_plots"]] <- count_plots
      cat("Updating 'Stat' object...\n")
      cat("- 'compute.descriptive' slot updated.\n")
      return(object)
    }
    return(count_plots)
  }, error = function(e) {
    cat("An error occurred:\n")
    cat("Error message:", conditionMessage(e), "\n")
    return(object)
  })
}

#' Generate Violin Plots for Numeric Variables
#'
#' This function creates faceted violin plots for selected numeric variables,
#' grouped by a categorical variable. It supports custom group colors via
#' `group_colors` or a palette name via `palette_name` (RColorBrewer or wesanderson).
#'
#' @param data A data frame containing the variables to plot.
#' @param vars_per_plot Number of variables to display per facet page.
#'   Default is 1.
#' @param save_dir Directory to save the generated plots. If NULL and
#'   save_plots is TRUE, a default directory is used.
#' @param palette_name Name of the color palette to use when `group_colors`
#'   is not provided. Supports RColorBrewer palettes (e.g., "Set1") or
#'   wesanderson palettes (e.g., "Zissou1"). Default is "Zissou1".
#' @param group_col Name of the column defining the grouping variable.
#'   If NULL, no grouping is applied. Default is "group".
#' @param max_unique_values Maximum number of unique values to treat a
#'   variable as categorical (used internally). Default is 5.
#' @param sub_var Character vector of specific numeric variables to plot.
#'   If NULL, all numeric variables are used.
#' @param save_plots Logical; if TRUE, saves each plot as a PDF.
#'   Default is TRUE.
#' @param plot_width Width of the saved plot (inches). Default is 5.
#' @param plot_height Height of the saved plot (inches). Default is 5.
#' @param base_size Base font size for the plot theme. Default is 14.
#' @param stat_method Statistical test for pairwise comparisons.
#'   Default is "wilcox.test".
#' @param paired_comparison Logical; if TRUE, adds pairwise significance
#'   labels. Default is TRUE.
#' @param group_colors An optional named character vector specifying custom
#'   colors for each group. Names must match the group levels. If provided,
#'   it overrides `palette_name`. Default is NULL.
#'
#' @return A list of ggplot objects (one per facet page).
#' @export
#'
#' @examples
#' if (interactive()) {
#'   data(iris)
#'   violin_plots(iris, group_col = "Species", save_plots = FALSE)
#' }
violin_plots <- function(data,
                         vars_per_plot = 1,
                         save_dir = NULL,
                         palette_name = "Zissou1",
                         group_col = "group",
                         max_unique_values = 5,
                         sub_var = NULL,
                         save_plots = TRUE,
                         plot_width = 5,
                         plot_height = 5,
                         base_size = 14,
                         stat_method = "wilcox.test",
                         paired_comparison = TRUE,
                         group_colors = NULL) {
  
  # ---- Handle default save_dir ----
  if (save_plots && is.null(save_dir)) {
    if (exists(".get_viz_output_dir")) {
      save_dir <- .get_viz_output_dir("Stat")
    } else {
      save_dir <- file.path(".", "figures", "Stat")
    }
  }
  
  # ---- Validate group column ----
  if (!is.null(group_col) && !group_col %in% names(data)) {
    stop(paste("Column", group_col, "not found in data"))
  }
  
  # ---- Identify numeric variables ----
  variable_types <- diagnose_variable_type(data, group_col = group_col,
                                           max_unique_values = max_unique_values)
  num_cols <- variable_types$numeric_vars
  
  if (!is.null(sub_var)) {
    num_cols <- intersect(num_cols, sub_var)
  }
  
  if (length(num_cols) == 0) {
    stop("No numeric variables found after filtering with sub_var")
  }
  
  # ---- Reshape data ----
  if (is.null(group_col)) {
    melted_data <- reshape2::melt(data, measure.vars = num_cols)
    melted_data$group <- "All"
  } else {
    melted_data <- reshape2::melt(data, id.vars = group_col, measure.vars = num_cols)
    melted_data$group <- as.character(melted_data[[group_col]])
  }
  
  # ---- Split variables into groups for faceting ----
  num_vars <- length(num_cols)
  num_groups <- ceiling(num_vars / vars_per_plot)
  var_groups <- split(num_cols, rep(1:num_groups, each = vars_per_plot,
                                    length.out = num_vars))
  
  plot_list <- list()
  
  for (i in seq_along(var_groups)) {
    group_vars <- var_groups[[i]]
    group_data <- melted_data[melted_data$variable %in% group_vars, ]
    
    if (nrow(group_data) > 0) {
      # ---- Determine group levels ----
      groups <- unique(group_data$group)
      n_groups <- length(groups)
      
      # ---- Generate colors ----
      if (!is.null(group_colors)) {
        # Check if named vector and match groups
        if (!is.null(names(group_colors))) {
          matched <- intersect(names(group_colors), groups)
          if (length(matched) == n_groups) {
            pal <- group_colors[groups]   # reorder to match group order
          } else {
            warning("Names in 'group_colors' do not match all group levels. Falling back to palette.")
            pal <- NULL   # will be generated from palette_name
          }
        } else {
          # Unnamed: assume order matches group levels
          if (length(group_colors) == n_groups) {
            pal <- group_colors
            names(pal) <- groups
            message("Using 'group_colors' in the order of group levels: ",
                    paste(groups, collapse = ", "))
          } else {
            warning("Length of 'group_colors' does not match number of groups. Falling back to palette.")
            pal <- NULL
          }
        }
      } else {
        pal <- NULL
      }
      
      # If pal is NULL (either because group_colors not provided or mismatch),
      # generate from palette_name
      if (is.null(pal)) {
        # Try RColorBrewer first
        if (requireNamespace("RColorBrewer", quietly = TRUE) &&
            palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
          max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
          if (n_groups <= max_colors) {
            pal <- RColorBrewer::brewer.pal(n_groups, palette_name)
          } else {
            base_cols <- RColorBrewer::brewer.pal(max_colors, palette_name)
            pal <- grDevices::colorRampPalette(base_cols)(n_groups)
          }
        } else if (requireNamespace("wesanderson", quietly = TRUE)) {
          # Fallback to wesanderson
          tryCatch({
            pal <- wesanderson::wes_palette(palette_name, n = n_groups, type = "discrete")
          }, error = function(e) {
            warning("wesanderson palette '", palette_name, "' not found. Using rainbow colors.")
            pal <- grDevices::rainbow(n_groups)
          })
        } else {
          # Ultimate fallback
          warning("No known palette found. Using rainbow colors.")
          pal <- grDevices::rainbow(n_groups)
        }
        names(pal) <- groups
      }
      
      # ---- Build plot ----
      p <- ggplot2::ggplot(group_data,
                           ggplot2::aes(x = as.factor(group), y = value)) +
        ggplot2::geom_violin(ggplot2::aes(fill = group), scale = "area", alpha = 0.5) +
        ggplot2::geom_boxplot(width = 0.1, size = 0.7, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.3, color = "black", size = 0.5) +
        ggplot2::scale_fill_manual(values = pal) +
        ggplot2::facet_wrap(~variable, scales = "free_y", ncol = 1) +
        ggprism::theme_prism(base_size = base_size) +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(title = paste("Violin Plots - Part", i),
                      x = "Group", y = "Value")
      
      # ---- Add pairwise comparisons ----
      if (!is.null(group_col) && paired_comparison && n_groups >= 2) {
        comparisons <- if (n_groups >= 3) {
          list(c(groups[1], groups[3]), c(groups[2], groups[3]), c(groups[1], groups[2]))
        } else {
          utils::combn(groups, 2, simplify = FALSE)
        }
        if (length(comparisons) > 0) {
          p <- p + ggpubr::stat_compare_means(
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
      
      # ---- Save plot ----
      if (save_plots) {
        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }
        ggplot2::ggsave(
          filename = paste0(save_dir, "/violin_plot_part_", i, ".pdf"),
          plot = p,
          width = plot_width,
          height = plot_height,
          device = "pdf"
        )
        cat("Saved violin plot part", i, "\n")
      }
    }
  }
  
  return(plot_list)
}

#' Density Ridge Plots for Numeric Variables
#'
#' This function generates density ridge plots for numeric variables in the dataset,
#' with optional grouping by a specified column. It handles large numbers of variables
#' by splitting them into multiple plots.
#'
#' @param data A data frame containing the data to plot.
#' @param vars_per_plot Integer. The number of variables (columns) to include
#'   in each individual plot. Default is 1.
#' @param save_dir String. The directory to save the plots in. Default is `NULL`.
#'   When `save_plots = TRUE` and this is `NULL`, a default figure directory is used.
#' @param palette_name String. The name of the color palette to use. Default `"Zissou1"`.
#' @param group_col String. The name of the column in the data used to group
#'   the samples. Default is `"group"`.
#' @param max_unique_values Integer. The maximum number of unique values
#'   allowed for categorical variables. Defaults to 5.
#' @param sub_var Character vector. A subset of variable names (columns)
#'   to include in the plot. If `NULL`, all numeric variables will be included.
#' @param save_plots Logical. If `TRUE`, the plots will be saved as PDF
#'   files in the specified `save_dir`. Default is `TRUE`.
#' @param plot_width Numeric. The width of each saved plot (in inches). Default 5.
#' @param plot_height Numeric. The height of each saved plot (in inches). Default 5.
#' @param base_size Numeric. The base font size for the plot. Default 14.
#'
#' @returns A list of ggplot objects containing the generated density ridge plots.
#'   If `save_plots` is `TRUE`, the plots are also saved as PDF files.
#' @export
#' @examples
#' if (interactive()) {
#'   data(iris)
#'   density_ridge_plots(iris[, 1:4], group_col = NULL, save_plots = FALSE)
#' }
density_ridge_plots <- function(data,
                                vars_per_plot = 1,
                                save_dir = NULL,
                                palette_name = "Zissou1",
                                group_col = "group",
                                max_unique_values = 5,
                                sub_var = NULL,
                                save_plots = TRUE,
                                plot_width = 5,
                                plot_height = 5,
                                base_size = 14) {
  
  # ---- Handle default save_dir ----
  if (save_plots && is.null(save_dir)) {
    if (exists(".get_viz_output_dir")) {
      save_dir <- .get_viz_output_dir("Stat")
    } else {
      save_dir <- file.path(".", "figures", "Stat")
    }
  }
  
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
  
  melted_data <- reshape2::melt(data, id.vars = group_col, measure.vars = num_cols)
  
  if (!is.null(group_col)) {
    melted_data$group <- melted_data[[group_col]]
  } else {
    melted_data$group <- "Default Group"
  }
  
  num_vars <- length(num_cols)
  num_groups <- ceiling(num_vars / vars_per_plot)
  var_groups <- split(num_cols, rep(1:num_groups, each = vars_per_plot, length.out = num_vars))
  
  pal <- wesanderson::wes_palette(palette_name, 100, type = "continuous")
  plot_list <- list()
  
  for (i in seq_along(var_groups)) {
    group_vars <- var_groups[[i]]
    group_data <- melted_data[melted_data$variable %in% group_vars, ]
    group_data$group <- as.factor(group_data$group)
    
    if (nrow(group_data) > 0) {
      p <- ggplot2::ggplot(group_data,
                           ggplot2::aes(x = value, y = group, fill = ggplot2::after_stat(density))) +
        ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.00) +
        ggplot2::scale_fill_gradientn(colours = pal) +
        ggplot2::facet_wrap(~variable, scales = "free_x", ncol = 1) +
        ggprism::theme_prism(base_size = base_size) +
        ggplot2::labs(title = paste("Density Ridge Plots - Part", i),
                      x = "Value",
                      y = ifelse(!is.null(group_col), group_col, "Group"))
      
      plot_list[[i]] <- p
      
      if (save_plots) {
        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }
        ggplot2::ggsave(paste0(save_dir, "/density_ridge_plot_part_", i, ".pdf"),
                        plot = p,
                        width = plot_width,
                        height = plot_height,
                        device = "pdf")
        cat("Saved ridge plot part", i, "\n")
      }
    }
  }
  
  return(plot_list)
}


#' Descriptive Numeric Variable Plots (Violin or Ridge)
#'
#' This function generates descriptive plots for numeric variables in a dataset,
#' either violin plots or ridge density plots, grouped by a categorical variable.
#' It supports custom color specification for groups.
#'
#' @param object A `Stat` object or a data frame.
#' @param vars_per_plot Number of variables to display per plot (for faceting). Default is 1.
#' @param save_dir Directory to save plots. If NULL and save_plots is TRUE, default directory is used.
#' @param palette_name Name of the palette to use if `group_colors` is not provided.
#'   Supports RColorBrewer (e.g., "Set1") or wesanderson (e.g., "Zissou1"). Default is "Zissou1".
#' @param group_col Column name defining the grouping variable. If NULL, attempts to extract from Stat object.
#' @param max_unique_values Maximum unique values to treat a variable as categorical. Default is 5.
#' @param plot_type Type of plot: "violin" or "ridge". Default is "violin".
#' @param save_plots Logical; if TRUE, saves plots to disk. Default is TRUE.
#' @param plot_width Width of saved plots (inches). Default is 5.
#' @param plot_height Height of saved plots (inches). Default is 5.
#' @param base_size Base font size for plots. Default is 14.
#' @param sub_var Character vector of specific numeric variables to plot. If NULL, all numeric variables are used.
#' @param group_colors An optional named character vector specifying custom colors for each group.
#'   Names must match group levels. If provided, overrides `palette_name`.
#'
#' @return The updated object (if `object` is a Stat) with plots stored in `compute.descriptive` slot.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   stat <- CreateStatObject(raw.data = mtcars, group_col = "cyl")
#'   plot_numeric_descriptive(stat, save_plots = FALSE,group_col = "cyl")
#' }
plot_numeric_descriptive <- function(object,
                                     vars_per_plot = 1,
                                     save_dir = NULL,
                                     palette_name = "Zissou1",
                                     group_col = "group",
                                     max_unique_values = 5,
                                     plot_type = "violin",
                                     save_plots = TRUE,
                                     plot_width = 5,
                                     plot_height = 5,
                                     base_size = 14,
                                     sub_var = NULL,
                                     group_colors = NULL) {
  
  if (save_plots && is.null(save_dir)) {
    if (exists(".get_viz_output_dir")) {
      save_dir <- .get_viz_output_dir("Stat")
    } else {
      save_dir <- file.path(".", "figures", "Stat")
    }
  }
  
  if (inherits(object, "Stat")) {
    cat("Input is an object of class 'Stat'. Extracting data...\n")
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    }
    if (is.null(group_col)) {
      group_col <- object@group_col
    }
    if (length(group_col) == 0) {
      group_col <- NULL
      cat("No group_col found in Stat object. Using NULL.\n")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using it directly.\n")
    data <- object
    if (is.null(group_col)) {
      stop("When object is a data frame, 'group_col' must be specified.")
    }
  } else {
    stop("Input must be an object of class 'Stat' or a data frame\n")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input\n")
  }
  
  # diagnose_variable_type is an internal helper (assumed available)
  variable_types <- diagnose_variable_type(data, group_col = group_col,
                                           max_unique_values = max_unique_values)
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
                               sub_var = numeric_vars,
                               group_colors = group_colors)   # pass to internal
  } else if (plot_type == "ridge") {
    cat("Generating ridge density plots...\n")
    # density_ridge_plots also needs to support group_colors; if not, ignore or adapt.
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
    # Note: group_colors not passed to ridge as per request
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
#' data(iris)
#' types <- diagnose_variable_type(iris, group_col = "Species")
#' converted <- convert_variables(iris, types, save_data = FALSE)
#' str(converted)
convert_variables <- function(data, 
                              variable_types,
                              save_dir = NULL,
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
#' \code{\link{diagnose_variable_type}}. If the input is a 'Stat' object, the
#' conversion will update its \code{clean.data} slot. If the input is a data frame,
#' it will directly return the converted data frame.
#'
#' @param object An object of class \code{Stat} or a data frame. If a \code{Stat}
#'   object is provided, the function will update its \code{clean.data} slot.
#' @param group_col A string representing the column name that groups the data.
#'   Default is \code{"group"}. This column is used to determine the type of variables.
#' @param max_unique_values The maximum number of unique values allowed for a
#'   variable to be considered as numeric. Default is \code{5}.
#' @param save_dir A character string specifying the directory where the converted
#'   data will be saved as a CSV file. If \code{NULL}, the default output directory
#'   is used. Default is \code{NULL}.
#' @param save_data Logical. If \code{TRUE}, saves the converted data to a CSV file.
#'   Default is \code{TRUE}.
#' @param csv_filename A character string specifying the name of the CSV file to
#'   save the converted data. Default is \code{"clean_data.csv"}.
#'
#' @return If the input is a \code{Stat} object, returns the updated \code{Stat}
#'   object with the converted data stored in the \code{clean.data} slot.
#'   If the input is a data frame, returns the converted data frame.
#'
#' @seealso \code{\link{diagnose_variable_type}}, \code{\link{convert_variables}}
#'
#' @export
#'
#' @examples
#' stat <- CreateStatObject(raw.data = mtcars, group_col = "cyl")
#' stat <- stat_convert_variables(stat, save_data = FALSE)
#' str(stat@clean.data)
stat_convert_variables <- function(object,
                                   group_col = "group",
                                   max_unique_values = 5,
                                   save_dir = NULL,
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
#' data(iris)
#' encoded <- one_hot_encode(iris, group_col = "Species", save_data = FALSE)
#' head(encoded)
one_hot_encode <- function(data, 
                           group_col = "group",
                           max_unique_values = 5,
                           save_dir = NULL,
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
#' @param object An object of class `Stat` or a `data.frame`.
#' @param method Reserved for future use. Currently not used. Default is `1`.
#' @param group_col A string specifying the name of the grouping column. This column will be
#'        excluded from one-hot encoding and moved to the end of the result. Default is `"group"`.
#' @param max_unique_values The maximum number of unique values a variable can have to be
#'        considered for one-hot encoding. Default is `5`.
#' @param save_dir A character string specifying the directory where the encoded data
#'        will be saved as a CSV file. If `NULL`, uses the default output directory
#'        from `get_output_dir("StatObject", "Data")`. Default is `NULL`.
#' @param save_data Logical. If `TRUE`, saves the encoded data to a CSV file.
#'        Default is `TRUE`.
#' @param csv_filename A character string specifying the name of the CSV file to
#'        save the encoded data. Default is `"clean_data.csv"`.
#'
#' @return If the input is a `Stat` object, returns the updated `Stat` object with
#'         the one-hot encoded data stored in the `clean.data` slot. If the input is
#'         a data frame, returns the one-hot encoded data frame.
#'
#' @export
#'
#' @examples
#' stat <- CreateStatObject(clean.data = mtcars, group_col = "cyl")
#' stat <- stat_onehot_encode(stat, save_data = FALSE,group_col = "cyl")
#' head(stat@clean.data)
stat_onehot_encode <- function(object, 
                               method = 1, 
                               group_col = "group", 
                               max_unique_values = 5,
                               save_dir = NULL,
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
