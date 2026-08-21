#' Perform t-SNE
#'
#' @param data Data.
#' @param perplexity Perplexity.
#' @param dims Dims.
#' @param iterations Iterations.
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
#'   result <- perform_tsne(df, perplexity = 10, dims = 2)
#'   plot(result$tsne_df)
#' }
#' @export
perform_tsne <- function(data,
                         perplexity = 30,
                         dims = 2,
                         iterations = 1000) {
  .require_pkgs("Rtsne")
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
  tsne_result <- Rtsne::Rtsne(data, perplexity = perplexity, dims = dims, max_iter = iterations)
  tsne_df <- as.data.frame(tsne_result$Y)
  rownames(tsne_df) <- row_names
  
  colnames(tsne_df) <- c("Dimension 1", "Dimension 2")
  cat("t-SNE completed successfully.\n")
  re<-list(tsne_result=tsne_result,
           tsne_df=tsne_df)
  return(re)
}

#' Sub t-SNE Analyse
#'
#' @param object Subtyping object.
#' @param perplexity Perplexity.
#' @param dims Dims.
#' @param iterations Iterations.
#' @param use_scaled_data Logical.
#' @param seed Numeric.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object
#'   obj <- Sub_tsne_analyse(obj, perplexity = 20, dims = 2)
#' }
#' @export
Sub_tsne_analyse <- function(object,
                             perplexity = 30,
                             dims = 2,
                             iterations = 1000,
                             use_scaled_data = TRUE,
                             seed = 123) {
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
  set.seed(seed)
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
    return(re$tsne_df)
  }
}

#' Sub Plot t-SNE
#'
#' @param object Subtyping object.
#' @param palette_name Palette.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' demo_df <- data.frame(
#'  group = rep(c(0, 1), each = 30),
#'  feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'  feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'  feat3 = rnorm(60, 10, 2),
#'  feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#')
#' rownames(demo_df) <- paste0("S", 1:60)

#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                             group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' sub_obj <- Sub_normalize_process(sub_obj, normalize_method = "min_max")

#'  #Sub_plot_tsne colors points by cluster assignment, so a clustering step
#'  #must run first to populate info.data with cluster labels.
#' sub_obj <- Sub_kmeans_with_optimal_k(sub_obj, k.max = 3, save_plots = FALSE, seed = 1)

#' sub_obj <- Sub_tsne_analyse(sub_obj, perplexity = 10, use_scaled_data = TRUE, seed = 1)

#' p <- Sub_plot_tsne(sub_obj)
#' p
#' }
Sub_plot_tsne <- function(object,
                          palette_name = "AsteroidCity1",
                          save_plots = TRUE,
                          save_dir = file.path(get_output_dir("m3", "visualization_results"), "tsne"),
                          plot_width = 5,
                          plot_height = 5,
                          base_size = 14) {
  
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

  if (!is.null(clustered.data) && nrow(clustered.data) > 0) {
    common_rows <- intersect(rownames(tsne_df), rownames(clustered.data))
  } else {
    common_rows <- character(0)
  }

  if (length(common_rows) == 0 && inherits(object, "Subtyping")) {
    info <- object@info.data
    cl_cols <- grep("^cluster_", colnames(info), value = TRUE)
    if (length(cl_cols) > 0) {
      cl_col      <- cl_cols[1]
      common_rows <- intersect(rownames(tsne_df), rownames(info))
      tsne_df     <- tsne_df[common_rows, , drop = FALSE]
      cluster_labels <- info[common_rows, cl_col]
      cat("Note: clustered.data row names did not match tsne_df;",
          "using info.data$", cl_col, "instead.\n", sep = "")
    } else {
      stop("No cluster labels found in clustered.data or info.data.")
    }
  } else {
    tsne_df        <- tsne_df[common_rows, , drop = FALSE]
    clustered.data <- clustered.data[common_rows, , drop = FALSE]
    cluster_labels <- clustered.data$group
  }
  cat("Remaining rows after alignment:", length(common_rows), "\n")
  
  tsne_df$cluster <- factor(cluster_labels)
  
  
  n_clusters <- length(unique(cluster_labels))
  pal <- tryCatch(
    wesanderson::wes_palette(palette_name, n_clusters, type = "continuous"),
    error = function(e) scales::hue_pal()(n_clusters)
  )

  p <- ggplot2::ggplot(tsne_df,
                       ggplot2::aes(x = `Dimension 1`, y = `Dimension 2`,
                                    color = cluster)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal, name = "Cluster") +
    ggplot2::labs(title = "t-SNE Visualization with Clustering",
                  x = "Dimension 1", y = "Dimension 2") +
    .pub_theme(base_size)

  if (save_plots) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tsne_clustering_plot.pdf"),
                    p, width = plot_width, height = plot_height, dpi = 300)
    ggplot2::ggsave(file.path(save_dir, "tsne_clustering_plot.png"),
                    p, width = plot_width, height = plot_height, dpi = 150)
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

