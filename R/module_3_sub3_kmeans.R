#' K-means Clustering
#'
#' @param data Numeric data.
#' @param k Number of clusters.
#' @param nstart nstart.
#' @export
kmeans_clustering <- function(data, 
                              k,
                              nstart = 50) {
  kmeans_result <- kmeans(data, centers = k, nstart = nstart)
  return(kmeans_result$cluster)
}

#' K-means with Optimal K
#'
#' @param data Data frame.
#' @param palette_name Palette.
#' @param save_plots Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param nstart nstart.
#' @param seed Seed.
#' @importFrom factoextra fviz_nbclust
#' @export
kmeans_with_optimal_k <- function(data,
                                  palette_name = "Zissou1",
                                  save_plots = TRUE,
                                  save_dir = here('Subtyping', "kmeans_result"),
                                  plot_width = 5,
                                  plot_height = 5,
                                  base_size = 14,
                                  nstart = 50,
                                  seed = 123) {
  set.seed(seed)

  if (!is.data.frame(data) || !all(sapply(data, is.numeric))) {
    stop("Input data must be a numeric data frame.")
  }


  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }



  silhouette <- fviz_nbclust(data, kmeans, method = "silhouette", k.max = 15) +
    labs(title = "Silhouette Method for Optimal K") +
    scale_color_manual(values = wes_palette(palette_name)) +
    theme_minimal() +
    theme_prism(base_size = base_size)

  print(silhouette)

  if (save_plots) {
    ggsave(filename = file.path(save_dir, "silhouette_plot.pdf"),
           plot = silhouette,
           width = plot_width,
           height = plot_height,
           device = "pdf")
    cat("Silhouette plot saved to:", file.path(save_dir, "silhouette_plot.pdf"), "\n")
  }

  optimal_k <- as.numeric(silhouette$data[which.max(silhouette$data$y), "clusters"])
  cat("Optimal K value is:", optimal_k, "\n")

  cluster_labels <- kmeans(data, centers = optimal_k, nstart = nstart)$cluster

  clustered_data <- cbind(data, group = factor(cluster_labels))

  return(list(clustered_data = clustered_data,
              cluster_labels = cluster_labels,
              optimal_k = optimal_k,
              plot_k = silhouette))
}

#' Sub K-means with Optimal K
#'
#' @param object Subtyping object.
#' @param use_scaled_data Logical.
#' @param palette_name Palette.
#' @param save_plots Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param seed Seed.
#' @export
Sub_kmeans_with_optimal_k <- function(object,
                                      use_scaled_data = TRUE,  
                                      palette_name = "Zissou1",
                                      save_plots = TRUE,
                                      save_dir = here('Subtyping', "cluster_results", "kmeans_result"),
                                      plot_width = 5,
                                      plot_height = 5,
                                      base_size = 14,
                                      seed = 123) {

  if (inherits(object, "Subtyping")) {
    if (use_scaled_data) {
      data <- slot(object, "scale.data")
    } else {
      data <- slot(object, "clean.data")
    }
    raw_data<-slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }

  numeric_data <- data[sapply(data, is.numeric)]

  if (ncol(numeric_data) == 0) {
    stop("No numeric columns found in the data. Please provide numeric data.")
  }

  cat("Starting K-means clustering analysis...\n")

  kmeans_result <- kmeans_with_optimal_k(numeric_data,
                                         palette_name = palette_name,
                                         save_plots = save_plots,
                                         save_dir = save_dir,
                                         plot_width = plot_width,
                                         plot_height = plot_height,
                                         base_size = base_size,
                                         seed = seed)
  
  clustered_data<-kmeans_result$clustered_data
  raw_data_with_rownames <- raw_data %>% 
    tibble::rownames_to_column("sample_id")
  
  clustered_data_with_rownames <- clustered_data %>% 
    tibble::rownames_to_column("sample_id") %>% 
    dplyr::select(sample_id, group)
  
  combined_data <- raw_data_with_rownames %>% 
    dplyr::left_join(clustered_data_with_rownames, by = "sample_id") %>% 
    tibble::column_to_rownames("sample_id") 
  
  
  if (inherits(object, "Subtyping")) {
    object@cluster.results[["kmeans.result"]] <- kmeans_result
    object@Optimal.cluster <- kmeans_result$optimal_k
    object@clustered.data <- combined_data
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'cluster.results' slot updated.\n")
    cat("- 'Optimal.cluster' slot updated.\n")
    cat("- 'clustered.data' slot updated.\n")

    return(object)
  } else {
    return(kmeans_result)
  }
}
