#' Visualize Classification Accuracy Across Probability Thresholds
#'
#' This function analyzes and visualizes how classification accuracy and other metrics
#' vary across different probability thresholds for binary classification. It identifies
#' the optimal threshold that maximizes accuracy and creates an informative plot showing
#' the relationship between thresholds and performance metrics.
#'
#' @param true_labels A vector of true class labels (factor or character)
#' @param pred_prob A vector of predicted probabilities for the positive class
#' @param positive_class The class to consider as positive (defaults to first level if NULL)
#' @param palette_name Name of Wes Anderson color palette to use (default: "AsteroidCity1")
#' @param base_size Base font size for the plot (default: 14)
#' @param save_plots Logical indicating whether to save the plot (default: TRUE)
#' @param save_dir Directory path to save plots (default: here("ModelData", "best_cutoff"))
#' @param plot_width Width of the saved plot in inches (default: 5)
#' @param plot_height Height of the saved plot in inches (default: 5)
#'
#' @return A list containing:
#' \itemize{
#'   \item plot - The ggplot object showing metrics across thresholds
#'   \item best_threshold - The probability threshold that maximizes accuracy
#'   \item accuracy - The maximum accuracy achieved
#'   \item metrics_data - Data frame containing all calculated metrics at each threshold
#' }
#'
#' @import ggplot2
#' @importFrom wesanderson wes_palette
#' @importFrom here here
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data <- data.frame(
#'   true = factor(sample(c("Disease", "Healthy"), 100, replace = TRUE)),
#'   prob = runif(100)
#' )
#' 
#' results <- visualize_accuracy(
#'   true_labels = data$true,
#'   pred_prob = data$prob,
#'   positive_class = "Disease"
#' )
#' 
#' # Access results:
#' print(results$best_threshold)
#' print(results$accuracy)
#' plot(results$plot)
#' }
visualize_accuracy <- function(true_labels, 
                               pred_prob, 
                               positive_class = NULL,
                               palette_name = "AsteroidCity1",
                               base_size = 14,
                               save_plots = TRUE,
                               save_dir = here("ModelData", "best_cutoff"),
                               plot_width = 5,
                               plot_height = 5) {
  
  if (is.null(positive_class)) {
    positive_class <- levels(factor(true_labels))[1]
    cat(paste("Positive class not specified, using first level:", positive_class))
  }
  
  true_labels <- factor(true_labels)
  negative_class <- levels(true_labels)[levels(true_labels) != positive_class]
  
  thresholds <- seq(min(pred_prob), max(pred_prob), length.out = 100)
  
  metrics <- lapply(thresholds, function(thresh) {
    pred_class <- ifelse(pred_prob >= thresh, positive_class, negative_class)
    pred_class <- factor(pred_class, levels = levels(true_labels))
    
    cf <- table(factor(true_labels, levels = levels(true_labels)), pred_class)
    
    TN <- cf[1, 1]
    FP <- cf[1, 2]
    FN <- cf[2, 1]
    TP <- cf[2, 2]
    
    data.frame(
      Threshold = thresh,
      Accuracy = (TP + TN) / sum(cf),
      PPV = TP / (TP + FP),
      NPV = TN / (TN + FN),
      Sensitivity = TP / (TP + FN),
      Specificity = TN / (TN + FP)
    )
  }) %>% dplyr::bind_rows()
  
  best_row <- metrics[which.max(metrics$Accuracy), ]
  best_threshold <- best_row$Threshold
  best_accuracy <- best_row$Accuracy
  
  palette_colors <- wesanderson::wes_palette(palette_name)
  accuracy_color <- palette_colors[1]
  point_color <- "red"
  
  x_mid <- mean(range(metrics$Threshold))
  
  metrics_plot <- ggplot(metrics, aes(x = Threshold)) +
    geom_line(aes(y = Accuracy, color = "Accuracy"), size = 1.2) +
    geom_vline(xintercept = best_threshold, linetype = "dashed", color = "gray40") +
    annotate("point", x = best_threshold, y = best_accuracy, color = point_color, size = 3) +
    scale_color_manual(values = c("Accuracy" = accuracy_color)) +
    labs(
      title = "Performance Metrics Across Thresholds",
      subtitle = "Optimizing for Accuracy", 
      x = "Threshold", 
      y = "Metric Value",
      color = "Metric") +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1.1),face = "bold"),
          legend.position = "top",
          plot.subtitle = element_text(hjust = 0.5,
                                       size = rel(1.0),
                                       face = "bold")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    annotate("text", 
             x = x_mid,  
             y = 0.25,
             label = paste0("Best Threshold: ", round(best_threshold, 3),
                            "\nAccuracy: ", round(best_accuracy, 3)),
             hjust = 0.5, vjust = 0, size = 4, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = 0.1))
  
  print(metrics_plot)
  
  if (save_plots) {
    ggplot2::ggsave(filename = file.path(save_dir, "accuracy_metrics_plot.pdf"),
                    plot = metrics_plot,
                    width = plot_width, height = plot_height, 
                    device = "pdf")
    cat("\nPlot saved to: ", file.path(save_dir, "accuracy_metrics_plot.pdf"))
  }
  
  return(list(
    plot = metrics_plot,
    best_threshold = best_threshold,
    accuracy = best_accuracy,
    metrics_data = metrics
  ))
}

#' Visualize Positive Predictive Value (PPV) Across Probability Thresholds
#'
#' This function analyzes and visualizes how Positive Predictive Value (PPV) and other 
#' classification metrics vary across different probability thresholds for binary 
#' classification. It can either find the threshold that maximizes PPV or the one that 
#' gets closest to a target PPV value.
#'
#' @param true_labels A vector of true class labels (factor or character)
#' @param pred_prob A vector of predicted probabilities for the positive class
#' @param positive_class The class to consider as positive (defaults to first level if NULL)
#' @param target_value Optional target PPV value to approach (if NULL, maximizes PPV)
#' @param palette_name Name of Wes Anderson color palette to use (default: "AsteroidCity1")
#' @param base_size Base font size for the plot (default: 14)
#' @param save_plots Logical indicating whether to save the plot (default: TRUE)
#' @param save_dir Directory path to save plots (default: here("ModelData", "best_cutoff"))
#' @param plot_width Width of the saved plot in inches (default: 5)
#' @param plot_height Height of the saved plot in inches (default: 5)
#'
#' @return A list containing:
#' \itemize{
#'   \item plot - The ggplot object showing metrics across thresholds
#'   \item best_threshold - The optimal probability threshold
#'   \item best_ppv - The PPV value at optimal threshold
#'   \item metrics_data - Data frame containing all calculated metrics at each threshold
#'   \item optimization_note - Text describing optimization approach used
#' }
#'
#' @details
#' The function calculates the following metrics at each threshold:
#' \itemize{
#'   \item PPV (Positive Predictive Value/Precision): TP/(TP+FP)
#'   \item Accuracy: (TP+TN)/(TP+FP+TN+FN)
#'   \item NPV (Negative Predictive Value): TN/(TN+FN)
#'   \item Sensitivity (Recall/TPR): TP/(TP+FN)
#'   \item Specificity (TNR): TN/(TN+FP)
#' }
#' When target_value is provided, the function finds the threshold that gets closest to
#' the target PPV value. Otherwise, it finds the threshold that maximizes PPV.
#'
#' @import ggplot2 
#' @importFrom wesanderson wes_palette
#' @importFrom here here
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' data <- data.frame(
#'   true = factor(sample(c("Disease", "Healthy"), 100, replace = TRUE, 
#'                       prob = c(0.3, 0.7))),
#'   prob = runif(100)
#' )
#'
#' # Example 1: Find threshold that maximizes PPV
#' result_max <- visualize_ppv(
#'   true_labels = data$true,
#'   pred_prob = data$prob,
#'   positive_class = "Disease"
#' )
#' 
#' # Example 2: Find threshold closest to target PPV of 0.8
#' result_target <- visualize_ppv(
#'   true_labels = data$true,
#'   pred_prob = data$prob,
#'   positive_class = "Disease",
#'   target_value = 0.8
#' )
#' 
#' # Access results
#' print(result_max$best_threshold)
#' print(result_target$best_ppv)
#' 
#' # View plots
#' print(result_max$plot)
#' print(result_target$plot)
#' }
visualize_ppv <- function(true_labels, 
                          pred_prob, 
                          positive_class = NULL,
                          target_value = NULL,
                          palette_name = "AsteroidCity1",
                          base_size = 14,
                          save_plots = TRUE,
                          save_dir = here("ModelData", "best_cutoff"),
                          plot_width = 5,
                          plot_height = 5) {
  
  if (is.null(positive_class)) {
    positive_class <- levels(factor(true_labels))[1]
    cat(paste("Positive class not specified, using first level:", positive_class))
  }
  
  true_labels <- factor(true_labels)
  negative_class <- levels(true_labels)[levels(true_labels) != positive_class]
  
  thresholds <- seq(min(pred_prob), max(pred_prob), length.out = 100)
  
  metrics <- lapply(thresholds, function(thresh) {
    pred_class <- ifelse(pred_prob >= thresh, positive_class, negative_class)
    pred_class <- factor(pred_class, levels = levels(true_labels))
    
    cf <- table(factor(true_labels, levels = levels(true_labels)), pred_class)
    
    TN <- cf[1, 1]
    FP <- cf[1, 2]
    FN <- cf[2, 1]
    TP <- cf[2, 2]
    
    data.frame(
      Threshold = thresh,
      Accuracy = (TP + TN) / sum(cf),
      PPV = TP / (TP + FP),
      NPV = TN / (TN + FN),
      Sensitivity = TP / (TP + FN),
      Specificity = TN / (TN + FP)
    )
  }) %>% dplyr::bind_rows()
  
  if (!is.null(target_value)) {
    best_row <- metrics[which.min(abs(metrics$PPV - target_value)), ]
    optimization_note <- paste0("Closest to PPV=", target_value)
  } else {
    best_row <- metrics[which.max(metrics$PPV), ]
    optimization_note <- "Maximum PPV"
  }
  
  best_threshold <- best_row$Threshold
  best_ppv <- best_row$PPV
  
  palette_colors <- wesanderson::wes_palette(palette_name)
  ppv_color <- palette_colors[1]
  clost_color <- palette_colors[2]
  point_color <- "red"
  
  x_mid <- mean(range(metrics$Threshold))
  
  metrics_plot <- ggplot(metrics, aes(x = Threshold)) +
    geom_line(aes(y = PPV, color = "PPV"), size = 1.2) +
    geom_vline(xintercept = best_threshold, linetype = "dashed", color = "gray40") +
    annotate("point", x = best_threshold, y = best_ppv, color = point_color, size = 3) +
    scale_color_manual(values = c("PPV" = ppv_color)) +
    labs(
      title = "Performance Metrics Across Thresholds",
      subtitle = "Optimizing for PPV", 
      x = "Threshold", 
      y = "Metric Value",
      color = "Metric") +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                    size = rel(1.1)),
          legend.position = "top",
          plot.subtitle = element_text(hjust = 0.5, face = "bold",
                                       size = rel(1.0))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    annotate("text", 
             x = x_mid, 
             y = 0.25,
             label = paste0("Best Threshold: ", round(best_threshold, 3),
                            "\nPPV: ", round(best_ppv, 3),
                            "\n(", optimization_note, ")"),
             hjust = 0.5, vjust = 0, size = 4, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = 0.1))
  
  if (!is.null(target_value)) {
    metrics_plot <- metrics_plot +
      geom_hline(yintercept = target_value, 
                 linetype = "dotted",
                 color = clost_color,
                 size = 1)
  }
  
  print(metrics_plot)
  
  if (save_plots) {
    ggplot2::ggsave(filename = file.path(save_dir, "ppv_metrics_plot.pdf"),
                    plot = metrics_plot,
                    width = plot_width, height = plot_height, 
                    device = "pdf")
    cat("\nPlot saved to: ", file.path(save_dir, "ppv_metrics_plot.pdf"))
  }
  
  return(list(
    plot = metrics_plot,
    best_threshold = best_threshold,
    best_ppv = best_ppv,
    metrics_data = metrics,
    optimization_note = optimization_note
  ))
}


#' Visualize Negative Predictive Value (NPV) Across Probability Thresholds
#'
#' This function analyzes and visualizes how Negative Predictive Value (NPV) varies across 
#' different probability thresholds for binary classification. It can either find the 
#' threshold that maximizes NPV or the one that gets closest to a target NPV value.
#'
#' @param true_labels A vector of true class labels (factor or character)
#' @param pred_prob A vector of predicted probabilities for the positive class
#' @param positive_class The class to consider as positive (defaults to first level if NULL)
#' @param target_value Optional target NPV value to approach (if NULL, maximizes NPV)
#' @param palette_name Name of Wes Anderson color palette to use (default: "AsteroidCity1")
#' @param base_size Base font size for the plot (default: 14)
#' @param save_plots Logical indicating whether to save the plot (default: TRUE)
#' @param save_dir Directory path to save plots (default: here("ModelData", "best_cutoff"))
#' @param plot_width Width of the saved plot in inches (default: 5)
#' @param plot_height Height of the saved plot in inches (default: 5)
#'
#' @return A list containing:
#' \itemize{
#'   \item plot - The ggplot object showing NPV across thresholds
#'   \item best_threshold - The optimal probability threshold
#'   \item best_npv - The NPV value at optimal threshold
#'   \item metrics_data - Data frame containing NPV at each threshold
#'   \item optimization_note - Text describing optimization approach used
#' }
#'
#' @details
#' The function calculates Negative Predictive Value (NPV) at each threshold:
#' \itemize{
#'   \item NPV = TN / (TN + FN)
#' }
#' When target_value is provided, the function finds the threshold that gets closest to
#' the target NPV value. Otherwise, it finds the threshold that maximizes NPV.
#'
#' @import ggplot2 
#' @importFrom wesanderson wes_palette
#' @importFrom here here
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' data <- data.frame(
#'   true = factor(sample(c("Disease", "Healthy"), 100, replace = TRUE, 
#'                       prob = c(0.3, 0.7))),
#'   prob = runif(100)
#' )
#'
#' # Example 1: Find threshold that maximizes NPV
#' result_max <- visualize_npv(
#'   true_labels = data$true,
#'   pred_prob = data$prob,
#'   positive_class = "Disease"
#' )
#' 
#' # Example 2: Find threshold closest to target NPV of 0.9
#' result_target <- visualize_npv(
#'   true_labels = data$true,
#'   pred_prob = data$prob,
#'   positive_class = "Disease",
#'   target_value = 0.9
#' )
#' 
#' # Access results
#' print(result_max$best_threshold)
#' print(result_target$best_npv)
#' 
#' # View plots
#' print(result_max$plot)
#' print(result_target$plot)
#' }
visualize_npv <- function(true_labels, 
                          pred_prob, 
                          positive_class = NULL,
                          target_value = NULL,
                          palette_name = "AsteroidCity1",
                          base_size = 14,
                          save_plots = TRUE,
                          save_dir = here("ModelData", "best_cutoff"),
                          plot_width = 5,
                          plot_height = 5) {
  
  if (is.null(positive_class)) {
    positive_class <- levels(factor(true_labels))[1]
    cat(paste("Positive class not specified, using first level:", positive_class))
  }
  
  true_labels <- factor(true_labels)
  negative_class <- levels(true_labels)[levels(true_labels) != positive_class]
  
  thresholds <- seq(min(pred_prob), max(pred_prob), length.out = 100)
  
  metrics <- lapply(thresholds, function(thresh) {
    pred_class <- ifelse(pred_prob >= thresh, positive_class, negative_class)
    pred_class <- factor(pred_class, levels = levels(true_labels))
    
    cf <- table(factor(true_labels, levels = levels(true_labels)), pred_class)
    
    TN <- cf[1, 1]
    FP <- cf[1, 2]
    FN <- cf[2, 1]
    TP <- cf[2, 2]
    
    data.frame(
      Threshold = thresh,
      NPV = TN / (TN + FN)
    )
  }) %>% dplyr::bind_rows()
  
  if (!is.null(target_value)) {
    best_row <- metrics[which.min(abs(metrics$NPV - target_value)), ]
    optimization_note <- paste0("Closest to NPV=", target_value)
  } else {
    best_row <- metrics[which.max(metrics$NPV), ]
    optimization_note <- "Maximum NPV"
  }
  
  best_threshold <- best_row$Threshold
  best_npv <- best_row$NPV
  
  if (!requireNamespace("wesanderson", quietly = TRUE)) {
    stop("Package 'wesanderson' required for color palette. Please install it.")
  }
  palette_colors <- wesanderson::wes_palette(palette_name)
  npv_color <- palette_colors[1]
  point_color <- "red"
  clost_color <- palette_colors[2]
  
  x_mid <- mean(range(metrics$Threshold))
  
  metrics_plot <- ggplot(metrics, aes(x = Threshold)) +
    geom_line(aes(y = NPV, color = "NPV"), size = 1.2) +
    geom_vline(xintercept = best_threshold, linetype = "dashed", color = "gray40") +
    annotate("point", x = best_threshold, y = best_npv, color = point_color, size = 3) +
    scale_color_manual(values = c("NPV" = npv_color)) +
    labs(
      title = "Performance Metrics Across Thresholds",
      subtitle = "Optimizing for NPV", 
      x = "Threshold", 
      y = "Metric Value",
      color = "Metric") +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                    size = rel(1.1)),
          plot.subtitle = element_text(hjust = 0.5, face = "bold",
                                       size = rel(1.0)),
          axis.text = element_text(color = "black"),
          legend.position = "top") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    annotate("text", 
             x = x_mid,  
             y = 0.25,
             label = paste0("Best Threshold: ", round(best_threshold, 3),
                            "\nNPV: ", round(best_npv, 3),
                            "\n(", optimization_note, ")"),
             hjust = 0.5, vjust = 0, size = 4, fontface = "bold")
  
  if (!is.null(target_value)) {
    metrics_plot <- metrics_plot +
      geom_hline(yintercept = target_value, 
                 linetype = "dotted",
                 color = clost_color,
                 size = 1) 
    
  }
  
  print(metrics_plot)
  
  if (save_plots) {
    ggplot2::ggsave(filename = file.path(save_dir, "npv_metrics_plot.pdf"),
                    plot = metrics_plot,
                    width = plot_width, height = plot_height, 
                    device = "pdf")
    cat("\nPlot saved to: ", file.path(save_dir, "npv_metrics_plot.pdf"))
  }
  
  return(list(
    plot = metrics_plot,
    best_threshold = best_threshold,
    best_npv = best_npv,
    metrics_data = metrics,
    optimization_note = optimization_note
  ))
}



#' Visualize ROC Curve with Youden's Optimal Threshold
#'
#' This function creates a ROC curve plot and identifies the optimal threshold using
#' Youden's Index (J = Sensitivity + Specificity - 1). The optimal threshold maximizes
#' the difference between true positive rate and false positive rate.
#'
#' @param true_labels A vector of true class labels (factor or character)
#' @param pred_prob A vector of predicted probabilities for the positive class
#' @param positive_class The class to consider as positive (defaults to first level if NULL)
#' @param palette_name Name of Wes Anderson color palette to use (default: "AsteroidCity1")
#' @param base_size Base font size for the plot (default: 14)
#' @param save_plots Logical indicating whether to save the plot (default: TRUE)
#' @param save_dir Directory path to save plots (default: here("ModelData", "best_cutoff"))
#' @param plot_width Width of the saved plot in inches (default: 5)
#' @param plot_height Height of the saved plot in inches (default: 5)
#' @param target_value Optional target value for Youden's Index (currently unused)
#'
#' @return A list containing:
#' \itemize{
#'   \item plot - The ggplot object showing ROC curve with optimal point
#'   \item best_threshold - The optimal probability threshold
#'   \item sensitivity - Sensitivity at optimal threshold
#'   \item specificity - Specificity at optimal threshold
#'   \item youden_index - Youden's Index value at optimal threshold
#'   \item roc_object - The full ROC curve object from pROC
#' }
#'
#' @details
#' The function calculates the following:
#' \itemize{
#'   \item ROC curve using pROC package
#'   \item Optimal threshold using Youden's Index (maximizes Sensitivity + Specificity - 1)
#'   \item Plots ROC curve with optimal point marked
#'   \item Includes performance metrics in plot annotation
#' }
#'
#' @import ggplot2 
#' @importFrom wesanderson wes_palette
#' @importFrom here here
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' data <- data.frame(
#'   true = factor(sample(c("Disease", "Healthy"), 100, replace = TRUE, 
#'                       prob = c(0.3, 0.7))),
#'   prob = runif(100)
#' )
#'
#' # Generate ROC plot with Youden's optimal threshold
#' result <- visualize_youden(
#'   true_labels = data$true,
#'   pred_prob = data$prob,
#'   positive_class = "Disease"
#' )
#' 
#' # Access results
#' print(result$best_threshold)
#' print(result$youden_index)
#' 
#' # View plot
#' print(result$plot)
#' }
visualize_youden <- function(true_labels, 
                             pred_prob, 
                             positive_class = NULL,
                             palette_name = "AsteroidCity1",
                             base_size = 14,
                             save_plots = TRUE,
                             save_dir = here::here("ModelData", "best_cutoff"),
                             plot_width = 5,
                             plot_height = 5,
                             target_value =NULL) {
  
  
  true_labels <- factor(true_labels)
  if (is.null(positive_class)) {
    positive_class <- 1
    cat("\nPositive class not specified, using first level: ", positive_class)
  }
  negative_class <- levels(true_labels)[levels(true_labels) != positive_class]
  
  roc_obj <- safe_roc(
    true_labels = true_labels,
    pred_prob = pred_prob,
    positive_class = positive_class,
    negative_class = negative_class
  )
  
  best_coords <- pROC::coords(
    roc_obj, 
    "best", 
    best.method = "youden",
    ret = c("threshold", "specificity", "sensitivity")
  )
  
  best_threshold <- as.numeric(best_coords["threshold"])
  best_spec <- as.numeric(best_coords["specificity"])
  best_sens <- as.numeric(best_coords["sensitivity"])
  best_youden <- best_spec +best_sens-1
  
  roc_data <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Threshold = roc_obj$thresholds
  )
  
  palette_colors <- wesanderson::wes_palette(palette_name)
  roc_color <- palette_colors[1]
  youden_color <- "red"
  
  x_mid <- mean(range(roc_data$FPR))
  
  roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = roc_color, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    annotate(
      "point",
      x = 1 - best_spec,
      y = best_sens,
      color = youden_color,
      size = 4
    ) +
    annotate(
      "text",
      x = x_mid, 
      y = 0.25,
      label = sprintf(
        "Best Threshold = %.3f\nSensitivity = %.3f\nSpecificity = %.3f\nYouden Index = %.3f",
        best_threshold,
        best_sens,
        best_spec,
        best_youden
      ),
      hjust = 0.5,  
      vjust = 0,    
      size = 4,
      fontface = "bold"
    ) +
    labs(
      title = "ROC Curve with Youden's Optimal Threshold",
      subtitle = "Optimizing for Youden Index (Sensitivity + Specificity - 1)",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5, 
        face = "bold", 
        color = "black",
        size = rel(1.1)  
      ),
      plot.subtitle = element_text(
        hjust = 0.5, 
        face = "bold", 
        size = base_size * 0.9  
      ),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    ) +
    scale_x_continuous(
      limits = c(0, 1), 
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = 0.1)  
    ) +
    scale_y_continuous(
      limits = c(0, 1), 
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = 0.1)  
    )
  
  print(roc_plot)
  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    ggplot2::ggsave(
      filename = file.path(save_dir, "roc_youden_plot.pdf"),
      plot = roc_plot,
      width = plot_width,
      height = plot_height,
      device = "pdf"
    )
    cat("\nPlot saved to: ", file.path(save_dir, "roc_youden_plot.pdf"))
  }
  
  return(list(
    plot = roc_plot,
    best_threshold = best_threshold,
    sensitivity = best_sens,
    specificity = best_spec,
    youden_index = best_youden,
    roc_object = roc_obj
  ))
}

#' Safe ROC Curve Calculation with Direction Auto-Correction
#'
#' Calculates ROC curve while automatically detecting and correcting for cases where
#' the probability direction needs to be reversed (when AUC < 0.5). This ensures
#' the ROC curve is always calculated in the optimal direction.
#'
#' @param true_labels A vector of true class labels (factor or character)
#' @param pred_prob A vector of predicted probabilities for the positive class
#' @param positive_class The class to consider as positive (must match levels in true_labels)
#' @param negative_class The class to consider as negative (must match levels in true_labels)
#'
#' @return An object of class "roc" from the pROC package, representing the ROC curve
#' calculated in the optimal direction (either ">" or "<").
#'
#' @details
#' The function:
#' \itemize{
#'   \item First calculates ROC curve with default direction (">")
#'   \item If AUC < 0.5, recalculates with reversed direction ("<")
#'   \item Returns the ROC curve with higher AUC
#'   \item Issues a message if direction was corrected
#' }
#'
#' @importFrom pROC roc auc
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with correct direction
#' roc1 <- safe_roc(
#'   true_labels = c("Disease", "Healthy", "Disease"),
#'   pred_prob = c(0.9, 0.2, 0.8),
#'   positive_class = "Disease",
#'   negative_class = "Healthy"
#' )
#'
#' # Example needing direction correction
#' roc2 <- safe_roc(
#'   true_labels = c("Disease", "Healthy", "Disease"),
#'   pred_prob = c(0.1, 0.8, 0.2),  # Lower probabilities predict Disease
#'   positive_class = "Disease",
#'   negative_class = "Healthy"
#' )
#' }
safe_roc <- function(true_labels, pred_prob, positive_class, negative_class) {
  roc_auto <- roc(
    response = true_labels,
    predictor = pred_prob,
    levels = c(negative_class, positive_class),
    direction = ">" 
  )
  
  if (auc(roc_auto) < 0.5) {
    roc_reversed <- roc(
      response = true_labels,
      predictor = pred_prob,
      levels = c(negative_class, positive_class),
      direction = "<"  
    )
    
    if (auc(roc_reversed) > auc(roc_auto)) {
      message("Warning: Direction was auto-corrected to '<' (low prob = positive class)")
      return(roc_reversed)
    }
  }
  
  return(roc_auto)
}


#' Select Optimal Classification Threshold Using Specified Method
#'
#' This function provides a unified interface for selecting optimal classification thresholds
#' using different methods (accuracy, PPV, NPV, or Youden's index). It wraps around specialized
#' visualization functions and returns comprehensive results.
#'
#' @param true_labels A vector of true class labels (factor or character)
#' @param pred_prob A vector of predicted probabilities for the positive class
#' @param positive_class The class to consider as positive (defaults to first level if NULL)
#' @param method The optimization method to use: "accuracy", "ppv", "npv", or "youden" 
#'        (default: "accuracy")
#' @param target_value Optional target value for PPV/NPV methods (default: NULL)
#' @param palette_name Name of Wes Anderson color palette to use (default: "AsteroidCity1")
#' @param base_size Base font size for the plot (default: 14)
#' @param save_plots Logical indicating whether to save the plot (default: TRUE)
#' @param save_dir Directory path to save plots (default: here("ModelData", "best_cutoff"))
#' @param plot_width Width of the saved plot in inches (default: 5)
#' @param plot_height Height of the saved plot in inches (default: 5)
#' @param seed Random seed for reproducibility (default: 1234)
#'
#' @return A list containing:
#' \itemize{
#'   \item All elements returned by the underlying visualization function
#'   \item method - The optimization method used
#'   \item target_value - The target value (if specified)
#' }
#'
#' @details
#' The function supports four threshold selection methods:
#' \itemize{
#'   \item "accuracy" - Maximizes overall accuracy (default)
#'   \item "ppv" - Maximizes or approaches target Positive Predictive Value
#'   \item "npv" - Maximizes or approaches target Negative Predictive Value
#'   \item "youden" - Maximizes Youden's Index (Sensitivity + Specificity - 1)
#' }
#' When target_value is provided for PPV/NPV methods, finds the threshold closest to the target.
#'
#' @importFrom here here
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Maximize accuracy
#' result1 <- select_best_threshold(
#'   true_labels = sample(c("Disease", "Healthy"), 100, replace = TRUE),
#'   pred_prob = runif(100),
#'   method = "accuracy"
#' )
#'
#' # Example 2: Approach target PPV of 0.8
#' result2 <- select_best_threshold(
#'   true_labels = sample(c("Disease", "Healthy"), 100, replace = TRUE),
#'   pred_prob = runif(100),
#'   method = "ppv",
#'   target_value = 0.8
#' )
#'
#' # Access results
#' print(result1$best_threshold)
#' plot(result2$plot)
#' }
select_best_threshold <- function(true_labels, 
                                  pred_prob, 
                                  positive_class = NULL,
                                  method = c("accuracy", "ppv", "npv", "youden"),
                                  target_value = NULL,
                                  palette_name = "AsteroidCity1",
                                  base_size = 14,
                                  save_plots = TRUE,
                                  save_dir = here("ModelData", "best_cutoff"),
                                  plot_width = 5,
                                  plot_height = 5,
                                  seed =1234) {
  
  set.seed(seed)
  method <- match.arg(method)
  if (is.null(positive_class)) {
    positive_class <- levels(true_labels)[1]
  }
  
  result <- switch(method,
                   "accuracy" = visualize_accuracy(
                     true_labels = true_labels,
                     pred_prob = pred_prob,
                     positive_class = positive_class,
                     palette_name = palette_name,
                     base_size = base_size,
                     save_plots = save_plots,
                     save_dir = save_dir,
                     plot_width = plot_width,
                     plot_height = plot_height
                   ),
                   "ppv" = visualize_ppv(
                     true_labels = true_labels,
                     pred_prob = pred_prob,
                     positive_class = positive_class,
                     target_value = target_value,
                     palette_name = palette_name,
                     base_size = base_size,
                     save_plots = save_plots,
                     save_dir = save_dir,
                     plot_width = plot_width,
                     plot_height = plot_height
                   ),
                   "npv" = visualize_npv(
                     true_labels = true_labels,
                     pred_prob = pred_prob,
                     positive_class = positive_class,
                     target_value = target_value,
                     palette_name = palette_name,
                     base_size = base_size,
                     save_plots = save_plots,
                     save_dir = save_dir,
                     plot_width = plot_width,
                     plot_height = plot_height
                   ),
                   "youden" = visualize_youden(
                     true_labels = true_labels,
                     pred_prob = pred_prob,
                     positive_class = positive_class,
                     target_value = target_value,
                     palette_name = palette_name,
                     base_size = base_size,
                     save_plots = save_plots,
                     save_dir = save_dir,
                     plot_width = plot_width,
                     plot_height = plot_height
                   )
  )
  
  result$method <- method
  result$target_value <- target_value
  
  cat("\nThreshold Selection Summary:")
  cat("\nMethod: ", toupper(method))
  cat("\nBest Threshold: ", round(result$best_threshold, 4),"\n")
  
  if (!is.null(target_value)) {
    cat("\nTarget Value: ", target_value)
  }
  
  return(result)
}



#' Determine Optimal Threshold for Best Model
#'
#' This function calculates and visualizes the optimal classification threshold for the best model
#' stored in a Best_Model object, using specified optimization criteria. It updates the Best_Model
#' object with threshold information.
#'
#' @param object A Best_Model object containing a trained model and data
#' @param group_col Name of the column containing class labels (default: "group")
#' @param method Optimization method: "accuracy", "ppv", "npv", or "youden" (default: "accuracy")
#' @param target_value Target value for PPV/NPV methods (default: NULL)
#' @param palette_name Name of color palette for plots (default: "AsteroidCity1")
#' @param base_size Base font size for plots (default: 14)
#' @param save_plots Logical to save output plots (default: TRUE)
#' @param save_dir Directory to save plots (default: here("ModelData", "best_cutoff"))
#' @param plot_width Plot width in inches (default: 5)
#' @param plot_height Plot height in inches (default: 5)
#' @param positive_class Name of positive class (default: NULL, uses first factor level)
#'
#' @return The input Best_Model object with updated slots:
#' \itemize{
#'   \item best_threshold - Optimal threshold value
#'   \item threshold_info - List containing method, target_value, and plot
#' }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts the best model and training data from the Best_Model object
#'   \item Calculates predicted probabilities on training data
#'   \item Determines optimal threshold using specified method
#'   \item Updates the Best_Model object with threshold information
#'   \item Optionally saves threshold selection plot
#' }
#'
#' @importFrom here here
#' @export
#'
#' @examples
#' \dontrun{
#' # After creating a Best_Model object
#' best_model <- ModelThreshold(
#'   object = best_model,
#'   method = "youden",
#'   save_plots = TRUE
#' )
#' 
#' # Access updated threshold information
#' print(best_model@best_threshold$best_threshold)
#' plot(best_model@best_threshold$threshold_info$plot)
#' }
ModelThreshold <- function(object,
                           group_col = "group",
                           method = c("accuracy"),
                           target_value = NULL,
                           palette_name = "AsteroidCity1",
                           base_size = 14,
                           save_plots = TRUE,
                           save_dir = here::here("ModelData", "best_cutoff"),
                           plot_width = 5,
                           plot_height = 5,
                           positive_class =NULL) {
  
  if (!inherits(object, "Best_Model")) {
    stop("Input must be an object of class 'Best_Model'")
  }
  
  best_model<-object@best.model[[1]]
  best_model_type<-object@best.model.type
  group_col  <-  object@group_col
  data_sets <- object@filtered.set
  
  train_data <- data_sets$training
  
  predicted_values <- predict(best_model, newdata = train_data, type = "prob")[, 2]
  true_labels <- train_data[[group_col]]
  
  if (!is.factor(true_labels)) {
    true_labels <- factor(true_labels)
  }
  
  result <- select_best_threshold(
    true_labels = true_labels,
    pred_prob = predicted_values,
    positive_class = positive_class, 
    method = method,
    target_value = target_value,
    palette_name = palette_name,
    base_size = base_size,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height
  )
  
  object@best_threshold[["best_threshold"]] <- result$best_threshold
  object@best_threshold[["threshold_info"]] <- list(
    method = method,
    target_value = target_value,
    plot = result$plot
  )
  cat("Updating 'Best_Model' object...\n")
  cat("The 'Best_Model' object has been updated with the following slots:\n")
  cat("- 'best_threshold' slot updated.\n")
  
  return(object)
}



