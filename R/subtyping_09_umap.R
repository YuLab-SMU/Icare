#' Perform UMAP
#'
#' @param data Data.
#' @param dims Dims.
#' @param n_neighbors Neighbors.
#' @param min_dist Min dist.
#' @param metric Metric.
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
#'   result <- perform_umap(df, n_neighbors = 10, dims = 2)
#'   plot(result$umap_df)
#' }
#' @export
perform_umap <- function(data, 
                         dims = 2, 
                         n_neighbors = 15,
                         min_dist = 0.1,
                         metric = "euclidean") {
  .require_pkgs("umap")
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
  
  cfg              <- umap::umap.defaults
  cfg$n_components <- dims
  cfg$n_neighbors  <- n_neighbors
  cfg$min_dist     <- min_dist
  cfg$metric       <- metric

  umap_result <- umap::umap(data, config = cfg)
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
#' @param use_scaled_data Logical.
#' @param seed Numeric.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object
#'   obj <- Sub_umap_analyse(obj, n_neighbors = 10, dims = 2)
#' }
#' @export
Sub_umap_analyse <- function(object,
                             dims = 2,
                             n_neighbors = 15,
                             min_dist = 0.1,
                             metric = "euclidean",
                             use_scaled_data = TRUE,
                             seed = 123) {
  cat("Starting Sub_umap_analyse...\n")
  
  if (inherits(object, "Subtyping")) {
    data <- if (use_scaled_data) slot(object, "scale.data") else slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  set.seed(seed)
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
#' @examples
#' \dontrun{
#' set.seed(1)
#' demo_df <- data.frame(
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2),
#'   feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#' )
#' rownames(demo_df) <- paste0("S", 1:60)
#'
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' sub_obj <- Sub_normalize_process(sub_obj, normalize_method = "min_max")
#'
#' # Sub_plot_umap colors points by cluster assignment, so a clustering step
#' # must run first to populate info.data with cluster labels.
#' sub_obj <- Sub_kmeans_with_optimal_k(sub_obj, k.max = 3, save_plots = FALSE, seed = 1)
#' sub_obj <- Sub_umap_analyse(sub_obj, n_neighbors = 10, use_scaled_data = TRUE, seed = 1)
#'
#' p <- Sub_plot_umap(sub_obj)
#' print(p)
#' }
#' @export
Sub_plot_umap <- function(object,
                          palette_name = "AsteroidCity1",
                          save_plots = TRUE,
                          save_dir = file.path(get_output_dir("m3", "visualization_results"), "umap"),
                          plot_width = 5,
                          plot_height = 5,
                          base_size = 14) {
  
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
  common_rows <- intersect(rownames(umap_df), rownames(clustered.data))
  umap_df <- umap_df[common_rows, , drop = FALSE]
  clustered.data <- clustered.data[common_rows, , drop = FALSE]
  cat("Removed duplicates, remaining rows:", nrow(clustered.data), "\n")
  
  cluster_labels<-clustered.data$group
  
  umap_df$cluster <- factor(cluster_labels)
  cat("UMAP data updated with cluster labels. Cluster levels:", levels(umap_df$cluster), "\n")
  
  n_clusters <- length(unique(cluster_labels))
  pal <- tryCatch(
    wesanderson::wes_palette(palette_name, n_clusters, type = "continuous"),
    error = function(e) scales::hue_pal()(n_clusters)
  )

  p <- ggplot2::ggplot(umap_df,
                       ggplot2::aes(x = `Dimension 1`, y = `Dimension 2`,
                                    color = cluster)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal, name = "Cluster") +
    ggplot2::labs(title = "UMAP Visualization with Clustering",
                  x = "Dimension 1", y = "Dimension 2") +
    .pub_theme(base_size)

  if (save_plots) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    plot_file <- file.path(save_dir, "umap_clustering_plot.pdf")
    ggplot2::ggsave(plot_file, p, width = plot_width, height = plot_height, dpi = 300)
    ggplot2::ggsave(sub("\\.pdf$", ".png", plot_file), p,
                    width = plot_width, height = plot_height, dpi = 150)
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
