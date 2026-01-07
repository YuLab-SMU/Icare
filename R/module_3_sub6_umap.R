#' Perform UMAP
#'
#' @param data Data.
#' @param dims Dims.
#' @param n_neighbors Neighbors.
#' @param min_dist Min dist.
#' @param metric Metric.
#' @importFrom umap umap
#' @export
perform_umap <- function(data, 
                         dims = 2, 
                         n_neighbors = 15,
                         min_dist = 0.1,
                         metric = "euclidean") {
  cat("Starting UMAP...\n")
  
  if (!is.data.frame(data)) {
    stop("Input data must be of data frame type.")
  }
  
  if (nrow(data) <= n_neighbors) {
    stop("The number of rows in the data must be greater than the number of neighbors.")
  }
  
  if (!all(sapply(data, is.numeric))) {
    stop("Input data must be entirely numeric. Please check and process non-numeric data.")
  }
  
  
  data <- data[sapply(data, is.numeric)]
  
  if (ncol(data) == 0) {
    stop("Input data must be entirely numeric. Please check and process non-numeric data.")
  }
  
  unique_rows <- !duplicated(data)
  data <- data[unique_rows, ]
  row_names <- rownames(data)
  
  cat("Removed duplicates, remaining rows:", nrow(data), "\n")
  
  umap_result <- umap(data,
                      dims = dims,
                      n_neighbors = n_neighbors,
                      min_dist = min_dist,
                      metric = metric)
  umap_df <- as.data.frame(umap_result$layout)
  cat("UMAP completed successfully.\n")
  rownames(umap_df) <- row_names
  
  colnames(umap_df) <- c("Dimension 1", "Dimension 2")
  re<-list(umap_result=umap_result,
           umap_df=umap_df)
  return(re)
}

#' Sub UMAP Analyse
#'
#' @param object Subtyping object.
#' @param dims Dims.
#' @param n_neighbors Neighbors.
#' @param min_dist Min dist.
#' @param metric Metric.
#' @export
Sub_umap_analyse <- function(object,
                             dims = 2,
                             n_neighbors = 15,
                             min_dist = 0.1,
                             metric = "euclidean") {
  cat("Starting Sub_umap_analyse...\n")
  
  if (inherits(object, "Subtyping")) {
    data <- slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  umap_result <- perform_umap(data,
                              dims = dims,
                              n_neighbors = n_neighbors,
                              min_dist = min_dist,
                              metric = metric)
  
  
  if (inherits(object, "Subtyping")) {
    object@visualization.results[["umap.result"]] <- umap_result$umap_result
    object@visualization.results[["umap.df"]]<- umap_result$umap_df
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'visualization.results' slot updated.\n")
    return(object)
  } else {
    cat("Returning UMAP results as data frame.\n")
    return(umap_result$umap_df)
  }
}

#' Sub Plot UMAP
#'
#' @param object Subtyping object.
#' @param palette_name Palette.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param seed Seed.
#' @export
Sub_plot_umap <- function(object,
                          palette_name = "AsteroidCity1",
                          save_plots = TRUE,
                          save_dir = here('Subtyping', "visualization_results","umap"),
                          plot_width = 5,
                          plot_height = 5,
                          base_size = 14,
                          seed = 123) {
  
  set.seed(seed)
  
  if (inherits(object, "Subtyping")) {
    clustered.data <- object@clustered.data
    umap_df <- object@visualization.results[["umap.df"]]
  } else if (is.list(object)) {
    cluster_labels <- object$cluster_labels
    umap_df <- object$umap_df
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame.")
  }
  
  if (is.null(umap_df) || nrow(umap_df) == 0) {
    stop("No valid UMAP data found in the input.")
  }
  clustered.data <- clustered.data[!duplicated(clustered.data), ]
  cat("Removed duplicates, remaining rows:", nrow(clustered.data), "\n")
  
  cluster_labels<-clustered.data$group
  
  umap_df$cluster <- factor(cluster_labels)
  cat("UMAP data updated with cluster labels. Cluster levels:", levels(umap_df$cluster), "\n")
  
  p <- ggplot(umap_df, aes(x = `Dimension 1`, y = `Dimension 2`, color = cluster)) +
    geom_point() +
    theme_minimal(base_size = base_size) +
    theme_prism(base_size = base_size) +
    labs(title = "UMAP Visualization with Clustering",
         x = "Dimension 1",
         y = "Dimension 2") +
    scale_color_manual(values = wes_palette(palette_name)) +
    theme(legend.position = "right")
  
  print(p)
  cat("UMAP plot created and printed.\n")
  
  if (save_plots) {
    plot_file <- file.path(save_dir, "umap_clustering_plot.pdf")
    ggsave(filename = plot_file, plot = p, width = plot_width, height = plot_height,
           device = "pdf")
    cat("Plot saved to:", plot_file, "\n")
  }
  
  if (inherits(object, "Subtyping")) {
    object@visualization.results[["umap.plot"]] <- p
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'visualization.results' slot updated.\n")
    return(object)
  } else {
    cat("Returning UMAP plot.\n")
    return(p)
  }
}
