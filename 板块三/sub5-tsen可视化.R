library(Rtsne)
library(ggplot2)
perform_tsne <- function(data,
                         perplexity = 30,
                         dims = 2,
                         iterations = 1000) {
  cat("Starting t-SNE with perplexity:", perplexity, "\n")
  
  if (!is.data.frame(data)) {
    stop("Input data must be of data frame type.")
  }
  
  if (nrow(data) <= perplexity) {
    stop("The number of rows in the data must be greater than perplexity.")
  }
  
  data <- data[sapply(data, is.numeric)]
  
  if (ncol(data) == 0) {
    stop("Input data must be entirely numeric. Please check and process non-numeric data.")
  }
  
  unique_rows <- !duplicated(data)
  data <- data[unique_rows, ]
  row_names <- rownames(data)
  
  cat("Removed duplicates, remaining rows:", nrow(data), "\n")
  
  
  if (perplexity < 5 || perplexity > 50) {
    stop("Perplexity should be between 5 and 50.")
  }
  
  if (iterations < 250) {
    stop("Iterations should be greater than 250.")
  }
  
  cat("Performing t-SNE...\n")
  tsne_result <- Rtsne(data, perplexity = perplexity, dims = dims, max_iter = iterations)
  tsne_df <- as.data.frame(tsne_result$Y)
  rownames(tsne_df) <- row_names
  
  colnames(tsne_df) <- c("Dimension 1", "Dimension 2")
  cat("t-SNE completed successfully.\n")
  re<-list(tsne_result=tsne_result,
           tsne_df=tsne_df)
  return(re)
}



Sub_tsne_analyse <- function(object,
                             perplexity = 30,
                             dims = 2,
                             iterations = 1000,
                             use_scaled_data = TRUE) {
  cat("Starting Sub_tsne_analyse...\n")
  if (inherits(object, "Subtyping")) {
    if (use_scaled_data) {
      data <- slot(object, "scale.data")
    } else {
      data <- slot(object, "clean.data")
    }
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  # Ensure the data is numeric
  numeric_data <- data[sapply(data, is.numeric)]
  
  if (ncol(numeric_data) == 0) {
    stop("Input data must be entirely numeric. Please check and process non-numeric data.")
  }
  
  re <- perform_tsne(numeric_data,
                     perplexity = perplexity,
                     dims = dims,
                     iterations = iterations)
  
  
  if (inherits(object, "Subtyping")) {
    object@visualization.results[["tsne.result"]] <- re$tsne_result
    object@visualization.results[["tsne.df"]] <- re$tsne_df
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'visualization.results' slot updated.\n")
    return(object)
  } else {
    cat("Returning t-SNE results as data frame.\n")
    return(tsne_df)
  }
}


Sub_plot_tsne <- function(object,
                          palette_name = "AsteroidCity1",
                          save_plots = TRUE,
                          save_dir = here('Subtyping', "visualization_results","tsne"),
                          plot_width = 5,
                          plot_height = 5,
                          base_size = 14,
                          seed = 123) {
  
  set.seed(seed)
  if (inherits(object, "Subtyping")) {
    clustered.data <- object@clustered.data
    tsne_df <- object@visualization.results[["tsne.df"]]
  } else if (is.list(object)) {
    clustered.data <- object$clustered.data
    tsne_df <- object$tsne_df
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame.")
  }
  
  if (is.null(tsne_df) || nrow(tsne_df) == 0) {
    stop("No valid t-SNE data found in the input.")
  }
  clustered.data <- clustered.data[!duplicated(clustered.data), ]
  cat("Removed duplicates, remaining rows:", nrow(clustered.data), "\n")
  
  cluster_labels<-clustered.data$group
  
  tsne_df$cluster <- factor(cluster_labels)
  
  
  p <- ggplot(tsne_df, aes(x = `Dimension 1`, y = `Dimension 2`, color = cluster)) +
    geom_point() +
    scale_color_manual(values = wes_palette(palette_name)) +
    theme_minimal(base_size = base_size) +
    theme_prism(base_size = base_size) +
    labs(title = "T-SNE Visualization with Clustering",
         x = "Dimension 1",
         y = "Dimension 2") +
    theme(legend.position = "right")
  
  print(p)
  
  if (save_plots) {
    ggsave(filename = file.path(save_dir, "tsne_clustering_plot.pdf"),
           plot = p,
           width = plot_width,
           height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "tsne_clustering_plot.pdf"), "\n")
  }
  
  if (inherits(object, "Subtyping")) {
    object@visualization.results[["tsne.plot"]] <- p
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'visualization.results' slot updated.\n")
    return(object)
  } else {
    return(p)
  }
}
