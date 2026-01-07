#' Perform Batch Wilcoxon Test for Multiple Variables
#'
#' This function performs pairwise Wilcoxon tests between two groups for each numeric variable in the input matrix.
#' It returns a summary table with the test statistics (W), p-values, mean and median values for each group,
#' standard deviations, log2 fold change, and the significance of the tests based on the p-values.
#' It also adjusts the p-values using the Bonferroni correction.
#'
#' @param mat A data frame or matrix containing the data. The first column should be the grouping variable (`group_col`),
#' and the remaining columns should contain numeric variables for the Wilcoxon tests.
#' @param group_col A string representing the column name in `mat` that defines the grouping variable.
#' It is assumed that the grouping variable has exactly two unique levels (e.g., treatment vs control).
#'
#' @returns A data frame with the following columns:
#' - `W`: The Wilcoxon test statistic.
#' - `p`: The p-value from the Wilcoxon test.
#' - `mean_x`: The mean of the variable for the first group.
#' - `mean_y`: The mean of the variable for the second group.
#' - `median_x`: The median of the variable for the first group.
#' - `median_y`: The median of the variable for the second group.
#' - `p.adjust`: The p-value adjusted using the Bonferroni correction.
#' - `logFC`: The log2 fold change between the two groups.
#' - `change`: A categorical variable indicating whether the log2 fold change is "Up", "Down", or "Stable",
#' based on a threshold of logFC > 0.5 or < -0.5, and p < 0.05.
#'
#' @export
#'
#' @examples
#' # Example of using the batch_Wilcoxon function
#' # Assume 'data' is a data frame with numeric variables and a grouping variable 'group'
#' result <- batch_Wilcoxon(data, group_col = "group")
#' print(result)
batch_Wilcoxon <- function(mat, 
                           group_col = "group",
                           save_dir = here::here("StatObject","deg_info"),
                           save_data = TRUE,
                           csv_filename = "last_test_sig.csv"
) {
  test.fun <- function(dat, col) {
    index <- unique(dat[[group_col]])
    sigs <- wilcox.test(
      dat[dat[[group_col]] == index[1], col],
      dat[dat[[group_col]] == index[2], col]
    )
    tests <- data.frame(
      W = sigs$statistic,
      p = sigs$p.value,
      mean_x = mean(dat[dat[[group_col]] == index[1], col]),
      mean_y = mean(dat[dat[[group_col]] == index[2], col]),
      median_x = median(dat[dat[[group_col]] == index[1], col]),
      median_y = median(dat[dat[[group_col]] == index[2], col])
    )
    return(tests)
  }

  mat[[group_col]] <- as.factor(as.character(mat[[group_col]]))
  numeric_cols <- sapply(mat, is.numeric)
  numeric_cols[group_col] <- TRUE
  mat <- mat[, numeric_cols]

  tests <- do.call(rbind, lapply(colnames(mat)[-which(colnames(mat) == group_col)], function(x) test.fun(mat, x)))

  rownames(tests) <- colnames(mat)[-which(colnames(mat) == group_col)]
  test_sig <- tests[tests$p < 0.01, ]
  test_sig$p.adjust <- p.adjust(test_sig$p, method = "bonferroni")
  test_sig <- test_sig[order(test_sig$p), ]

  sd_file <- mat %>%
    group_by_at(group_col) %>%
    dplyr::summarise_all(sd) %>%
    t(.)
  colnames(sd_file) <- sd_file[1, ]
  sd_file <- as.data.frame(sd_file[-1, ])
  sd_file$id <- rownames(sd_file)
  colnames(sd_file)[1:2] <- paste("sd", colnames(sd_file)[1:2], sep = "_")
  mean_file <- mat %>%
    group_by_at(group_col) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::select(-group_col) %>%
    log2(.)
  logFC <- mean_file[2, ] - mean_file[1, ]
  logFC <- as.data.frame(t(logFC))
  logFC$id <- rownames(logFC)
  colnames(logFC)[1] <- "logFC"

  test_sig$id <- rownames(test_sig)
  last_test_sig <- merge(test_sig, sd_file, by = "id")
  last_test_sig <- merge(last_test_sig, logFC, by = "id")
  last_test_sig <- last_test_sig[order(last_test_sig$p), ]
  last_test_sig$change <- as.factor(
    ifelse(last_test_sig$p < 0.05, ifelse(last_test_sig$logFC > 0.5, "Up",
                                          ifelse(last_test_sig$logFC < -0.5, "Down", "Stable")
    ))
  )
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(last_test_sig, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  
  return(last_test_sig)
}

#' Perform Feature Selection Using Batch Wilcoxon Test
#'
#' This function performs a batch Wilcoxon test for feature selection, comparing different groups based on the specified
#' group column. The function accepts an object of class `Stat` or a data frame, extracts the relevant data based on the
#' `data_type` argument, and applies the batch Wilcoxon test to select significant features.
#' If the input is an object of class `Stat`, the results will be updated in the object. If the input is a data frame,
#' the result will be returned as a data frame.
#' @importFrom dplyr group_by_at
#' @param object An object of class `Stat` or a data frame. If an object of class `Stat`, the relevant data will be
#' extracted using the `data_type` argument. If a data frame, it is used directly.
#' @param group_col A string representing the column name in `object` or the data frame that defines the grouping variable.
#' This column should contain at least two distinct non-missing values for comparison.
#' @param data_type A string indicating the type of data to extract from the `Stat` object. Options are `"clean"` for clean
#' data or `"scale"` for scaled data. Default is `"scale"`.
#'
#' @returns If the input is an object of class `Stat`, it returns the updated `Stat` object with the `var.result` slot
#' updated with the results of the feature selection. If the input is a data frame, it returns a data frame containing the
#' significant features based on the Wilcoxon test.
#'
#' @export
#'
#' @examples
#' # Example of using the stat_var_feature function with a data frame
#' result <- stat_var_feature(data, group_col = "group")
#' print(result)
#'
#' # Example of using the stat_var_feature function with a 'Stat' object
#' stat_object <- stat_var_feature(stat_object, group_col = "group", data_type = "clean")
#' print(stat_object)
stat_var_feature <- function(object,
                             group_col = "group",
                             save_dir = here::here("StatObject","deg_info"),
                             data_type = "scale",
                             save_data = TRUE,
                             csv_filename = "last_test_sig.csv") {
  if (inherits(object, "Stat")) {
    data <- switch(data_type,
                   clean = ExtractCleanData(object),
                   scale = ExtractScaleData(object),
                   stop("Invalid data_type. Choose either 'clean' or 'scale'."))

    group_col <- slot(object, "group_col")
    if (length(group_col) == 0) {
      group_col <- NULL
    }

  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (!group_col %in% colnames(data)) {
    stop(paste("The specified group column", group_col, "does not exist in the data."))
  }

  group_values <- data[[group_col]]

  if (any(is.na(group_values)) || length(unique(group_values)) < 2) {
    stop(paste("The group column", group_col, "must have at least two distinct non-missing values for comparison."))
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input.")
  }

  data <- data.frame(lapply(data, function(x) if (is.numeric(x)) as.numeric(x) else x))
  cat("Starting batch Wilcoxon test...\n")
  last_test_sig <- batch_Wilcoxon(data,
                                  group_col=group_col,
                                  save_dir =save_dir,
                                  save_data = save_data,
                                  csv_filename = csv_filename)

  if (inherits(object, 'Stat')) {
    object@var.result <- list(last_test_sig = last_test_sig)

    cat("Feature selection completed. Number of significant features: ", nrow(last_test_sig), "\n")
    cat("Updating 'Stat' object...\n")
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'var.result' slot updated.\n")

    return(object)
  } else {
    cat("Number of significant features: ", nrow(last_test_sig), "\n")
    return(last_test_sig)
  }
}





#' Extract Last Significant Test Results
#'
#' This function extracts the last significant test results stored in the `var.result` slot of a `Stat` object. If
#' the `last_test_sig` result is not available, it returns `NULL`.
#'
#' @param object An object of class `Stat`. The function attempts to extract the `last_test_sig` from the `var.result` slot.
#'
#' @returns Returns the last significant test results stored in the `last_test_sig` slot of the `var.result` list.
#' If the `last_test_sig` is not found, it returns `NULL`.
#'
#' @export
#'
#' @examples
#' # Example of extracting the last significant test results from a 'Stat' object
#' last_sig <- ExtractLastTestSig(stat_object)
#' print(last_sig)
ExtractLastTestSig <- function(object) {
  last_test_sig <- tryCatch(object@var.result[["last_test_sig"]],
                            error = function(e) NULL)
  return(last_test_sig)
}

#' Plot DEG Radar Chart
#'
#' This function generates a radar chart visualizing the differential expression (DEG) results from a data frame.
#' It uses a set of custom parameters to customize the appearance, such as colors, size, and grouping. The chart can be
#' saved to a specified directory as a PDF file.
#'
#'
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr pivot_longer
#' @import wesanderson
#' @import here
#' @param df A data frame containing the DEG results with columns for ID, means, log-fold change, p-values,
#' adjusted p-values, and other necessary data.
#' @param palette_name A string indicating the name of the color palette to use for the chart .
#' @param x_col The column name to use as the x-axis (default: "id").
#' @param y_cols A vector of column names to be plotted on the y-axis (default: c("mean_x", "mean_y")).
#' @param size_col The column name for the size of the points (default: "logFC").
#' @param color_col The column name for the color of the points, typically for p-values (default: "logp").
#' @param fill_col The column name for the fill color, which defines the grouping (default: "change").
#' @param p_adjust_col The column name for the adjusted p-values (default: "p.adjust").
#' @param plot_width The width of the plot (default: 5).
#' @param plot_height The height of the plot (default: 5).
#' @param save_dir The directory where the radar chart image will be saved (default: "here('StatObject', 'deg_info')").
#' @param base_size The base font size for the plot (default: 14).
#' @param title The title of the radar chart (default: "Radar Chart Title").
#'
#' @returns The ggplot object representing the radar chart.
#' @export
#'
#' @examples
#' # Example usage of the function
#' plot_deg_radarchart(df, palette_name = "Zissou1", title = "Differential Expression Radar")
plot_deg_radarchart <- function(df,
                                palette_name = "Zissou1",
                                x_col = "id",
                                y_cols = c("mean_x", "mean_y"),
                                size_col = "logFC",
                                color_col = "logp",
                                fill_col = "change",
                                p_adjust_col = "p.adjust",
                                plot_width = 5,
                                plot_height = 5,
                                save_dir = here("StatObject", "deg_info"),
                                base_size = 14,
                                title = "Radar Chart Title") {

  colors <- wes_palette(n = 3, name = palette_name, type = "continuous")

  colors <- as.list(colors)


  primary_color <- colors[[1]]
  secondary_color <- colors[[2]]


  df <- df %>%
    mutate(
      logp = -log10(get(p_adjust_col)),
      group = factor(get(fill_col), levels = c("Down", "Stable"))
    )

  df_long <- df %>%
    pivot_longer(cols = all_of(y_cols), names_to = "var", values_to = "value")

  plot <- ggplot(df_long, aes(x = get(x_col), y = value, fill = group)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.7) +
    geom_point(aes(size = abs(get(size_col)), color = get(color_col)), position = position_dodge(width = 0.9)) +
    coord_polar(start = 0) +
    scale_fill_manual(values = c("Down" = primary_color, "Stable" = secondary_color)) +
    scale_color_gradient(low = primary_color, high = secondary_color) +
    scale_size_continuous(range = c(2, 10)) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_line(color = "gray50"),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
      legend.position = "right",
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 12, color = "black"),
      panel.grid = element_line(color = "gray80")
    ) +
    labs(fill = "Change", color = "-log10(p.adjust)", size = "LogFC", title = title) +
    geom_hline(yintercept = 0, color = "gray60", linetype = "dashed")


  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(save_dir, "radar_chart.pdf"), plot = plot, width = plot_width, height = plot_height,
         device = "pdf")
  cat("Radar chart saved to: ", file.path(save_dir, "radar_chart.pdf"))

  return(plot)
}

#' Generate Variable Feature Radar Chart
#'
#' This function generates a radar chart for variable features from either a "Stat" object or a data frame.
#' It uses the `plot_deg_radarchart` function to visualize the data and saves the plot as an image.
#' The chart is customized based on various input parameters such as the color palette, plot size, and title.
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr pivot_longer
#' @import wesanderson
#' @import here
#' @param object An object of class "Stat" or a data frame containing the variable features to plot.
#' @param group_col The name of the column representing the groups in the data (default: "group").
#' @param palette_name A string indicating the color palette for the plot .
#' @param plot_width The width of the plot (default: 5).
#' @param plot_height The height of the plot (default: 5).
#' @param save_dir The directory where the plot image will be saved (default: "here('StatObject', 'deg_info')").
#' @param base_size The base font size for the plot (default: 14).
#' @param title The title for the radar chart (default: "Variable Feature Radar Chart").
#'
#' @returns If the input is a "Stat" object, the updated "Stat" object with the plot saved in `var.result`.
#' If the input is a data frame, it returns the generated radar chart plot.
#' @export
#'
#' @examples
#' # Example usage with a "Stat" object:
#' VarFeature_radarchart(stat_object, palette_name = "Zissou1", title = "Variable Feature Radar")
#'
#' # Example usage with a data frame:
#' VarFeature_radarchart(df, palette_name = "Zissou1", title = "Variable Feature Radar")
VarFeature_radarchart <- function(object,
                                  group_col = "group",
                                  palette_name = "Zissou1",
                                  plot_width = 5,
                                  plot_height = 5,
                                  save_dir = here("StatObject", "deg_info"),
                                  base_size = 14,
                                  title = "Variable Feature Radar Chart") {
  cat("Generating variable feature radar chart...\n")

  if (inherits(object, "Stat")) {
    mat <- ExtractLastTestSig(object)
    group_col = slot(object, "group_col")
    if (length(group_col) == 0) {
      group_col <- NULL
    }
  } else if (is.data.frame(object)) {
    mat <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (is.null(mat) || nrow(mat) == 0) {
    stop("No valid data found in the input.")
  }

  plot <- plot_deg_radarchart(mat,
                              palette_name = palette_name,
                              plot_width = plot_width,
                              plot_height = plot_height,
                              save_dir = save_dir,
                              base_size = base_size,
                              title = title)

  print(plot)

  if (inherits(object, "Stat")) {
    object@var.result[["VarFeaturePlot"]] <- plot
    return(object)
  }

  return(plot)
}


#' Generate Volcano Plot for Differential Expression
#'
#' This function generates a volcano plot based on the log-fold change and adjusted p-values of differential expression results.
#' The plot highlights significant and non-significant genes using different colors and includes labels for the significant genes.
#' The volcano plot is saved as an image in the specified directory.
#' @importFrom dplyr mutate %>%
#' @param last_test_sig A data frame containing the results of differential expression analysis. It should have columns
#'                      for log-fold change and adjusted p-values.
#' @param logFC_col The name of the column representing log-fold change values (default: "logFC").
#' @param p_adjust_col The name of the column representing adjusted p-values (default: "p.adjust").
#' @param title The title of the volcano plot (default: "Volcano Plot").
#' @param palette_name A string indicating the color palette to use for the plot.
#' @param save_dir The directory to save the volcano plot image (default: "here('StatObject', 'deg_info')").
#' @param plot_width The width of the plot in inches (default: 5).
#' @param plot_height The height of the plot in inches (default: 5).
#' @param base_size The base font size for the plot (default: 14).
#'
#' @returns A `ggplot` object representing the volcano plot.
#' @export
#'
#' @examples
#' # Example usage:
#' plot_deg_volcano(last_test_sig = deg_results, logFC_col = "logFC", p_adjust_col = "p.adjust")
#'
#' plot_deg_volcano(last_test_sig = deg_results, title = "Custom Volcano Plot", palette_name = "Zissou1")
#'
plot_deg_volcano <- function(last_test_sig,
                             logFC_col = "logFC",
                             p_adjust_col = "p.adjust",
                             title = "Volcano Plot",
                             palette_name = "Zissou1",
                             save_dir = here("StatObject", "deg_info"),
                             plot_width = 5,
                             plot_height = 5,
                             base_size = 14) {

  if (is.null(last_test_sig) || nrow(last_test_sig) == 0) {
    stop("No valid data found in last_test_sig.")
  }

  last_test_sig <- last_test_sig %>%
    mutate(log_p = -log10(get(p_adjust_col)))


  colors <- wes_palette(n = 3, name = palette_name, type = "continuous")

  colors <- as.list(colors)
  significant_color <- colors[[1]]
  not_significant_color <- colors[[2]]

  volcano_plot <- ggplot(last_test_sig, aes_string(x = logFC_col, y = "log_p")) +
    geom_point(aes(color = ifelse(get(p_adjust_col) < 0.05, "Significant", "Not Significant")),
               alpha = 0.6, size = 3) +
    scale_color_manual(values = c("Significant" = significant_color, "Not Significant" = not_significant_color)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") +
    labs(title = title, x = "Log Fold Change", y = "-Log10 Adjusted P-value") +
    theme_minimal(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      panel.grid = element_line(color = "gray80"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2)
    ) +
    geom_text(aes(label = id), vjust = -1, hjust = 0.5, size = 3, check_overlap = TRUE, color = "black")


  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    cat("Created directory: ", save_dir, "\n")
  }

  ggsave(filename = file.path(save_dir, "volcano_plot.pdf"), plot = volcano_plot,
         width = plot_width, height = plot_height,
         device = "pdf")

  cat("Volcano plot saved to: ", file.path(save_dir, "volcano_plot.pdf"), "\n")

  return(volcano_plot)
}


#' Generate Volcano Plot for Variable Features
#'
#' This function generates a volcano plot based on the log-fold change and adjusted p-values of variable features from
#' a 'Stat' object or a data frame. It then saves the volcano plot as an image in the specified directory.
#'
#' @param object An object of class 'Stat' or a data frame containing the results of differential expression analysis.
#'               It should include columns for log-fold change and adjusted p-values.
#' @param logFC_col The name of the column representing log-fold change values (default: "logFC").
#' @param p_adjust_col The name of the column representing adjusted p-values (default: "p.adjust").
#' @param title The title of the volcano plot (default: "Volcano Plot").
#' @param palette_name A string indicating the color palette to use for the plot.
#' @param plot_width The width of the plot in inches (default: 5).
#' @param plot_height The height of the plot in inches (default: 5).
#' @param save_dir The directory to save the volcano plot image (default: "here('StatObject', 'deg_info')").
#' @param base_size The base font size for the plot (default: 14).
#'
#' @returns A `ggplot` object representing the volcano plot. If the input is a 'Stat' object, it also updates the
#'          'Stat' object with the volcano plot.
#' @export
#'
#' @examples
#' # Example usage:
#' VarFeature_volcano(object = stat_object, logFC_col = "logFC", p_adjust_col = "p.adjust")
#'
#' VarFeature_volcano(object = stat_object, title = "Custom Volcano Plot", palette_name = "Zissou1")
VarFeature_volcano <- function(object,
                               logFC_col = "logFC",
                               p_adjust_col = "p.adjust",
                               title = "Volcano Plot",
                               palette_name = "Zissou1",
                               plot_width = 5,
                               plot_height = 5,
                               save_dir = here("StatObject", "deg_info"),
                               base_size = 14) {

  if (inherits(object, "Stat")) {
    last_test_sig <- ExtractLastTestSig(object)
  } else if (is.data.frame(object)) {
    last_test_sig <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (is.null(last_test_sig) || nrow(last_test_sig) == 0) {
    stop("No valid data found in last_test_sig.")
  }

  volcano_plot <- plot_deg_volcano(last_test_sig,
                                   logFC_col = logFC_col,
                                   p_adjust_col = p_adjust_col,
                                   title = title,
                                   palette_name = palette_name,
                                   plot_width = plot_width,
                                   plot_height = plot_height,
                                   base_size = base_size)

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(save_dir, "volcano_plot.pdf"), plot = volcano_plot,
         width = plot_width, height = plot_height,device = "pdf")

  if (inherits(object, "Stat")) {
    object@var.result[["VolcanoPlot"]] <- volcano_plot
    return(object)
  }

  return(volcano_plot)
}





#' Generate Violin Plot for Differentially Expressed Genes
#'
#' This function generates a violin plot for the top differentially expressed genes based on their log-fold change
#' and statistical significance. It also calculates the Wilcoxon test p-values for group comparisons and annotates
#' the plot with significance stars. The plot can be saved to a specified directory.
#' @importFrom dplyr %>% group_by arrange desc top_n filter mutate case_when
#' @importFrom dplyr summarise left_join
#' @importFrom tidyr pivot_longer
#'
#' @import stats
#' @param last_test_sig A data frame containing the results of differential expression analysis with columns
#'                      such as 'change' (e.g., 'Stable', 'Upregulated', 'Downregulated') and 'logFC' (log-fold change).
#' @param data A data frame containing the expression data with rows as samples and columns as features.
#' @param control The label of the control group (default: "health").
#' @param case The label of the case group (default: "cancer").
#' @param top_n The number of top significant features to display (default: 5).
#' @param palette_name The name of the color palette for the plot.
#' @param name_identity The identity for the type of analysis, default is "deg".
#' @param save_plots Logical value to indicate whether to save the plot (default: TRUE).
#' @param save_dir The directory to save the plot (default: "here('StatObject', 'deg_info')").
#' @param plot_width The width of the plot in inches (default: 5).
#' @param plot_height The height of the plot in inches (default: 5).
#' @param base_size The base font size for the plot (default: 14).
#' @param title The title of the plot (default: "Violin Plot").
#' @param group_col The column name used to group the data (default: 'group').
#'
#' @returns A ggplot object representing the violin plot.
#' @export
#'
#' @examples
#' # Example usage:
#' plot_deg_violinplot(last_test_sig = sig_results, data = expression_data)
#' plot_deg_violinplot(last_test_sig = sig_results, data = expression_data, top_n = 10)
plot_deg_violinplot <- function(last_test_sig,
                                data,
                                control = 'health',
                                case = 'cancer',
                                top_n = 5,
                                palette_name = "Royal1",
                                name_identity = 'deg',
                                save_plots = TRUE,
                                save_dir = here("StatObject", "deg_info"),
                                plot_width = 5,
                                plot_height = 5,
                                base_size = 14,
                                title = "Violin Plot",
                                group_col = 'group') {

  cat("Data columns:", colnames(data), "\n")

  cat("Filtering significant features...\n")
  left <- last_test_sig[last_test_sig$change != 'Stable', ]

  for_label <- left %>%
    group_by(change) %>%
    arrange(desc(abs(logFC)), .by_group = TRUE) %>%
    top_n(n = top_n, wt = abs(logFC)) %>%
    dplyr::filter(!is.na(id))

  cat("Filtered features (for_label):\n")
  print(for_label)

  if (nrow(for_label) == 0) {
    cat("No significant features found. Exiting the function.\n")
    return(NULL)
  }

  available_genes <- nrow(left)
  n <- min(top_n * 2, available_genes)

  selected_columns <- c(for_label$id[1:n], group_col)
  cat("Selected columns:", selected_columns, "\n")
  missing_columns <- setdiff(selected_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop(paste("Missing columns in data:", paste(missing_columns, collapse = ", ")))
  }

  box_test <- data[, selected_columns]
  box_test <- pivot_longer(box_test, cols = 1:n, values_to = 'value', names_to = 'id')
  box_test[[group_col]] <- as.factor(box_test[[group_col]])

  cat("Calculating significance differences...\n")
  p_values <- box_test %>%
    group_by(id) %>%
    summarise(p_value = wilcox.test(value ~ .data[[group_col]])$p.value, .groups = 'drop') %>%
    mutate(significance = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ))

  y_positions <- box_test %>%
    group_by(id) %>%
    summarise(y_position = max(log10(abs(value) + 1e-10)) + 0.1)

  p_values <- p_values %>%
    left_join(y_positions, by = "id")

  box_test <- left_join(box_test, p_values)
  index <- box_test$id[order(box_test$y_position, decreasing = FALSE)]
  box_test$id <- factor(box_test$id, levels = unique(index))

  cat("Drawing violin plot...\n")
  p <- ggplot(data = box_test, aes(x = id, y = log10(value), fill = .data[[group_col]])) +
    geom_violin(trim = FALSE, linewidth = 0.3) +
    geom_text(data = box_test, aes(x = id, y = y_position, label = significance),
              size = 3, color = "black") +
    scale_fill_manual(values = wes_palette(palette_name)) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.key.size = unit(0.4, 'cm'),
      legend.text = element_text(size = base_size * 0.5),
      legend.title = element_text(size = base_size * 0.6, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(size = base_size * 0.8, face = "bold")
    ) +
    labs(title = title) +
    labs(x = "Features", y = "Log10(Value)")

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    output_path <- file.path(save_dir, "violinplot.pdf")
    ggsave(filename = output_path, plot = p, height = plot_height, width = plot_width, device = "pdf")

    cat("Plot saved to:", output_path, "\n")
  }

  return(p)
}

#' Generate Violin Plot for Differential Features
#'
#' This function generates a violin plot to visualize the distribution of the top differentially expressed (DE) features between two groups (e.g., control vs case) based on statistical significance. The function allows users to customize the number of features, color palette, and save the plot to a specified directory.
#'
#' @param object The input object, which can either be of class 'Stat' or a data frame.
#' @param control The label for the control group in the data (default is 'health').
#' @param case The label for the case group in the data (default is 'cancer').
#' @param top_n The number of top features to plot based on differential expression (default is 5).
#' @param palette_name The color palette to use for the plot.
#' @param name_identity The identifier for the features, used in plot titles or labels (default is 'deg').
#' @param data_type Specifies whether to use "clean" or "scaled" data (default is "clean").
#' @param save_dir The directory where the plot will be saved (default is `here("StatObject", "deg_info")`).
#' @param plot_width The width of the plot in centimeters (default is 5).
#' @param plot_height The height of the plot in centimeters (default is 5).
#' @param base_size The base font size for the plot (default is 14).
#'
#' @returns If the input is a 'Stat' object, the function returns the updated 'Stat' object with the violin plot stored in the 'var.result' slot. If the input is a data frame, the function returns the generated violin plot.
#' @export
#'
#' @examples
#' # Example usage with a 'Stat' object:
#' var_plot <- VarFeature_violinplot(object = stat_object, control = 'healthy', case = 'cancer', top_n = 10)
#'
#' # Example usage with a data frame:
#' var_plot_df <- VarFeature_violinplot(object = df, control = 'healthy', case = 'cancer', top_n = 5)
#'
VarFeature_violinplot <- function(object,
                                  control='health',
                                  case='cancer',
                                  top_n = 5,
                                  palette_name = "Royal1",
                                  name_identity = 'deg',
                                  data_type = "clean",
                                  save_dir = here("StatObject", "deg_info"),
                                  plot_width = 5,
                                  plot_height = 5,
                                  base_size = 10) {

  if (inherits(object, "Stat")) {
    last_test_sig <- ExtractLastTestSig(object)

    if (inherits(object, "Stat")) {
      group_col <- slot(object, "group_col")
      if (length(group_col) == 0 || !group_col %in% colnames(data)) {
        group_col <- NULL
      }
      cat("Extracting data from 'Stat' object...\n")
      data <- if (data_type == "clean") {
        ExtractCleanData(object)
      } else {
        ExtractScaleData(object)
      }
    } else if (is.data.frame(object)) {
      data <- object
    } else {
      stop("Input must be an object of class 'Stat' or a data frame.")
    }
  } else if (is.data.frame(object)) {
    last_test_sig <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (is.null(last_test_sig) || nrow(last_test_sig) == 0) {
    stop("No valid data found in last_test_sig.")
  }

  violinplot <- plot_deg_violinplot(last_test_sig = last_test_sig,
                                    data = data,
                                    control = control,
                                    case = case,
                                    top_n = top_n,
                                    name_identity = name_identity,
                                    save_plots = TRUE,
                                    save_dir = save_dir,
                                    plot_width = plot_width,
                                    plot_height = plot_height,
                                    base_size = base_size,
                                    palette_name = palette_name)

  if (inherits(object, "Stat")) {
    object@var.result[["violinplot"]] <- violinplot
    return(object)
  }

  return(violinplot)
}

#' Plot ROC Curve for Differential Expression Gene Analysis
#'
#' This function generates a ROC curve for the top features identified in differential expression analysis.
#' It filters the significant features, calculates the ROC curve, and visualizes the AUC (Area Under the Curve)
#' for each feature. The plot can be saved to a specified directory.
#'
#' @importFrom dplyr %>% group_by arrange desc top_n filter
#' @importFrom pROC roc ggroc
#' @param deg_test A data frame containing the results of differential expression analysis.
#'                 It should have a column `change` indicating whether a feature is "Stable" or changed.
#' @param mat_test A data frame containing the expression values of features across samples.
#'                 It should have columns representing features (e.g., genes) and rows representing samples.
#' @param group_col The column in `mat_test` that indicates the grouping of samples (e.g., control vs case).
#'                  Default is `'group'`.
#' @param control The label used to indicate the control group in the data. Default is `'health'`.
#' @param case The label used to indicate the case group in the data. Default is `'lung cancer'`.
#' @param top_n The number of top features to consider based on log fold change. Default is 5.
#' @param palette_name The name of the palette to be used for the ROC curves.
#' @param name_identity The name or identifier used for the features in the `deg_test` data frame. Default is `'deg'`.
#' @param save_plots A logical value indicating whether the plot should be saved. Default is `TRUE`.
#' @param save_dir The directory path where the plot will be saved. Default is `here("StatObject", "deg_info")`.
#' @param plot_width The width of the saved plot in cm. Default is 5.
#' @param plot_height The height of the saved plot in cm. Default is 5.
#' @param base_size The base font size for the plot. Default is 14.
#' @param title The title of the plot. Default is `'ROC Curve'`.
#'
#' @returns A `ggplot` object representing the ROC curve for the top features.
#'          If `save_plots` is `TRUE`, the plot is saved to the specified directory.
#' @export
#'
#' @examples
#' plot_deg_Roc_plot(deg_test = deg_results, mat_test = expression_data)
plot_deg_Roc_plot <- function(deg_test,
                              mat_test,
                              group_col = 'group',
                              control = 'health',
                              case = 'lung cancer',
                              top_n = 5,
                              palette_name = "Royal1",
                              name_identity = 'deg',
                              save_plots = TRUE,
                              save_dir = here("StatObject", "deg_info"),
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 10,
                              title = 'ROC Curve') {

  left <- deg_test[deg_test$change != 'Stable', ]

  for_label <- left %>%
    group_by(change) %>%
    arrange(desc(abs(logFC)), .by_group = TRUE) %>%
    top_n(n = top_n, wt = abs(logFC)) %>%
    dplyr::filter(!is.na(id))

  cat("Filtered features (for_label):\n")
  print(for_label)

  if (nrow(for_label) == 0) {
    cat("No significant features found. Exiting the function.\n")
    return(NULL)
  }

  if (!all(for_label$id %in% colnames(mat_test))) {
    stop("Some features in for_label$id are not found in mat_test.")
  }

  roc_tem <- list()
  auc_tem <- numeric()

  for (j in 1:min(nrow(for_label), top_n * 2)) {
    id <- for_label$id[j]
    cat("Processing feature:", id, "\n")

    if (!is.numeric(mat_test[[id]])) {
      stop(paste("Feature", id, "is not numeric."))
    }

    roc_tem[[j]] <- roc(mat_test[[group_col]], mat_test[[id]])
    auc_tem[j] <- round(roc_tem[[j]]$auc, 2)
  }

  colors <- as.vector(wes_palette(n = length(roc_tem), name = palette_name, type = "discrete"))

  p <- ggroc(roc_tem) +
    theme_classic() +
    theme(legend.position = 'none',
          text = element_text(size = base_size),
          plot.title = element_text(size = base_size * 0.8, face = "bold")) +  # 标题字体大小
    geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
    scale_color_manual(values = wes_palette(palette_name)) +
    labs(title = title)

  for (j in 1:length(auc_tem)) {
    p <- p + annotate("text", x = 0.25, y = 0.01 + (j - 1) * 0.05,
                      label = paste(for_label$id[j], "AUC:", auc_tem[j]),
                      size = 2, color = colors[j])
  }

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    output_path <- file.path(save_dir, "roc_plot.pdf")
    ggsave(filename = output_path, plot = p, height = plot_height, width = plot_width, device = "pdf")
    cat("Saved plots to: ", output_path, "\n")
  }

  return(p)
}


#' Plot ROC Curve for Variable Features in Statistical Object
#'
#' This function generates a ROC curve for the top variable features identified in a statistical analysis.
#' It uses the results from a `Stat` object or a data frame, computes the ROC curves for the significant features,
#' and visualizes the AUC (Area Under the Curve) for each feature. The plot can be saved to a specified directory.
#'
#' @param object An object of class `'Stat'` or a data frame containing the statistical test results.
#'               If the input is a `Stat` object, it must contain the `last_test_sig` and data (either clean or scaled).
#' @param group_col The column in `data` or `object` representing the grouping variable (e.g., control vs case).
#'                  Default is `'group'`.
#' @param control The label used to indicate the control group in the data. Default is `'health'`.
#' @param case The label used to indicate the case group in the data. Default is `'cancer'`.
#' @param top_n The number of top features to consider based on log fold change. Default is 5.
#' @param palette_name The name of the palette to be used for the ROC curves.
#' @param name_identity The name or identifier used for the features in the `last_test_sig` data frame. Default is `'deg'`.
#' @param data_type The type of data to extract from the `Stat` object. Can be either `'clean'` or `'scale'`. Default is `'clean'`.
#' @param save_dir The directory path where the plot will be saved. Default is `here("StatObject", "deg_info")`.
#' @param plot_width The width of the saved plot in cm. Default is 5.
#' @param plot_height The height of the saved plot in cm. Default is 5.
#' @param base_size The base font size for the plot. Default is 14.
#'
#' @returns The updated `Stat` object with the ROC plot stored in the `var.result` slot if the input is a `Stat` object.
#'          If the input is a data frame, the ROC plot is returned directly.
#' @export
#'
#' @examples
#' VarFeature_ROC(object = stat_object)
#' VarFeature_ROC(object = df_data)
VarFeature_ROC <- function(object,
                           group_col = 'group',
                           control = 'health',
                           case = 'cancer',
                           top_n = 5,
                           palette_name = "Royal1",
                           name_identity = 'deg',
                           data_type = "clean",
                           save_dir = here("StatObject", "deg_info"),
                           plot_width = 5,
                           plot_height = 5,
                           base_size = 10) {

  if (inherits(object, "Stat")) {
    last_test_sig <- ExtractLastTestSig(object)

    if (missing(group_col) || is.null(group_col)) {
      group_col <- slot(object, "group_col")
    }

    data <- if (data_type == "clean") {
      ExtractCleanData(object)
    } else {
      ExtractScaleData(object)
    }

    if (!group_col %in% colnames(data)) {
      stop(paste("Group column", group_col, "not found in data. Available columns:",
                 paste(colnames(data), collapse = ", ")))
    }

  } else if (is.data.frame(object)) {
    if (!group_col %in% colnames(object)) {
      stop(paste("Group column", group_col, "not found in data frame."))
    }
    last_test_sig <- object$last_test_sig
    data <- object$data
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }
  if (is.null(last_test_sig) || nrow(last_test_sig) == 0) {
    stop("No valid data found in last_test_sig.")
  }

  roc_plot <- plot_deg_Roc_plot(deg_test = last_test_sig,
                                mat_test = data,
                                group_col = group_col,
                                control = control,
                                case = case,
                                top_n = top_n,
                                palette_name = palette_name,
                                name_identity = name_identity,
                                save_plots = TRUE,
                                save_dir = save_dir,
                                plot_width = plot_width,
                                plot_height = plot_height,
                                base_size = base_size)

  if (inherits(object, "Stat")) {
    object@var.result[["Rocplot"]] <- roc_plot
    return(object)
  }

  return(roc_plot)
}
