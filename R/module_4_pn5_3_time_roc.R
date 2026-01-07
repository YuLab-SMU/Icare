#' Generate Time-Dependent ROC Curve
#'
#' Generates time-dependent ROC curves for survival data.
#'
#' @param data Data frame.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param var_col Variable column (marker).
#' @param palette_name Palette name.
#' @param time_unit Time unit ("days", "years", "months").
#' @param save_dir Save directory.
#' @param plot_width Plot width.
#' @param plot_height Plot height.
#' @param base_size Base font size.
#' @param save_plots Logical.
#' @param save_format File format ("pdf", "png", etc.).
#'
#' @return List containing ROC results and the plot.
#' @export
generate_time_dependent_roc <- function(data,
                                        time_col,
                                        status_col,
                                        var_col,
                                        palette_name = "AsteroidCity1",
                                        time_unit = "months",
                                        save_dir = here::here('PrognosiX', "univariate_analysis"),
                                        plot_width = 5,
                                        plot_height = 5,
                                        base_size = 14,
                                        save_plots = TRUE,
                                        save_format = "pdf"){


  time <- data[[time_col]]
  status <- data[[status_col]]

  marker_data <- data[[var_col]]

  if (time_unit == "days") {
    time <- time * 365
  } else if (time_unit == "years") {
    time <- time / 12
  }

  time_points <- quantile(time, probs = c(0.25, 0.5, 0.75))
  time_points <- round(time_points, 2)

  roc_results <- timeROC::timeROC(
    T = time,
    delta = status,
    marker = marker_data,
    cause =1,
    weighting = "marginal",
    times = time_points
  )


  roc_plot <- ggplot2::ggplot()

  for (i in 1:length(time_points)) {
    time_point <- time_points[i]

    if (time_unit == "days") {
      time_label <- paste("t=", round(time_point * 365, 0), " days", sep = "")
    } else if (time_unit == "years") {
      time_label <- paste("t=", round(time_point / 12, 2), " years", sep = "")
    } else {
      time_label <- paste("t=", time_point, " months", sep = "")
    }

    FPR_1 <- roc_results[["FP_1"]][,i]
    FPR_2 <- roc_results[["FP_2"]][,i]

    TPR <- 1 - FPR_2

    Specificity <- 1 - FPR_1

    temp_data <- data.frame(Specificity = Specificity, Sensitivity = TPR, TimePoint = time_label)

    if (nrow(temp_data) > 0) {
      roc_plot <- roc_plot +
        geom_path(data = temp_data,
                  aes(x = Specificity, y = Sensitivity, color = TimePoint),
                  size = 1.5)
    }
  }

  roc_plot <- roc_plot +
    scale_color_manual(values = wesanderson::wes_palette(palette_name)) +
    labs(
      title = paste("Time-dependent ROC Curves by", var_col),
      subtitle = paste("AUCs:", paste0("   t=", time_points, ": ", round(roc_results$AUC_1, 2), collapse = "; ")),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Time Points"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.25, "cm"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, face = "italic"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10)
    )
  print(roc_plot)
  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat(paste("Created directory:", save_dir, "\n"))
    }
    ggsave(filename = file.path(save_dir, paste0("time_dependent_roc_by_", var_col,".",save_format)),
           plot = roc_plot,
           device = save_format,
           width = plot_width,
           height = plot_height)
    cat(paste("Plot saved as:", file.path(save_dir, paste0("time_dependent_roc_by_", var_col,".",save_format)), "\n"))
  }



  list(
    ROC_Results = roc_results,
    Plot = roc_plot
  )
}


#' Wrapper for Plotting Time-Dependent ROC
#'
#' Wrapper function to plot time-dependent ROC for PrognosiX object or data frame.
#'
#' @param object PrognosiX object or data frame.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param var_col Variable column.
#' @param binwidth Bin width (unused here but kept for signature).
#' @param palette_name Palette name.
#' @param save_plots Logical.
#' @param save_dir Save directory.
#' @param file_name Filename (unused here).
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param use_subgroup_data Logical.
#'
#' @return Updated PrognosiX object or plot object.
#' @export
plot_var_roc_plot <- function(object,
                              time_col = "time",
                              status_col = "status",
                              var_col,
                              binwidth = 10,
                              palette_name = "AsteroidCity1",
                              save_plots = TRUE,
                              save_dir = here::here('PrognosiX', "univariate_analysis"),
                              file_name = "survival_distribution_plot.pdf",
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 14,
                              use_subgroup_data = FALSE) {

  if (inherits(object, 'PrognosiX')) {
    if (use_subgroup_data) {
      data <- methods::slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
      data <- convert_to_numeric(data)
      data <- convert_two_class_to_binary(data)

    } else {
      data <- methods::slot(object, "survival.data")
      cat("Using original survival data...\n")
    }
    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")

    if (is.null(data) || nrow(data) == 0) {
      stop("The survival.data in the PrognosiX object is empty.")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  result <- generate_time_dependent_roc(
    data = data,
    time_col = time_col,
    status_col = status_col,
    var_col = var_col,
    palette_name = palette_name,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    save_plots = save_plots
  )

  result_string <- paste0("roc_time_", var_col)

  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object...\n")

    object@univariate.analysis[["roc_time_results"]][[result_string]]<- result

    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }

  cat("Plotting function execution completed.\n")
  return(result$Plot)
}
