#' Calinski-Harabasz Index
#'
#' @param clustered.data Clustered data.
#' @param label_col label column.
#' @return Numeric CH index value.
#' @export
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(30), y = rnorm(30), group = sample(1:3, 30, replace = TRUE))
#'   ch <- calinski_harabasz(df, label_col = "group")
#'   cat("Calinski-Harabasz Index:", ch, "\n")
#' }
calinski_harabasz <- function(clustered.data,label_col = "group") {
  X <- clustered.data[, !names(clustered.data) %in% label_col, drop = FALSE]
  labels <- clustered.data[[label_col]]
  
  n_samples <- nrow(X)
  unique_labels <- unique(labels)
  n_clusters <- length(unique_labels)
  
  if (n_clusters == 1) {
    cat("Warning: All samples belong to a single cluster. Returning 0.\n")
    return(0)
  }
  
  global_center <- colMeans(X)
  cat("Global center calculated.\n")
  
  cluster_centers <- matrix(NA, nrow = n_clusters, ncol = ncol(X))
  intra_disp <- 0
  label_to_index <- match(labels, unique_labels)
  cat("Cluster labels mapped to continuous indices.\n")
  
  for (k in 1:n_clusters) {
    cluster_points <- X[label_to_index == k, , drop = FALSE]
    if (nrow(cluster_points) == 0) {
      cat(sprintf("Warning: Cluster %d is empty. Using global center.\n", k))
      cluster_centers[k, ] <- global_center  
      next
    }
    
    center_k <- colMeans(cluster_points)
    cluster_centers[k, ] <- center_k
    intra_disp <- intra_disp + sum(rowSums((cluster_points - center_k)^2))
  }
  cat("Cluster centers and within-cluster dispersion calculated.\n")
  
  n_points_per_cluster <- as.vector(table(label_to_index))
  inter_disp <- sum(n_points_per_cluster * rowSums((cluster_centers - global_center)^2))
  cat("Between-cluster dispersion calculated.\n")
  
  if (intra_disp == 0) {
    cat("Warning: Zero within-cluster dispersion. Returning 0.\n")
    return(0)
  }
  
  chi_score <- (inter_disp / intra_disp) * ((n_samples - n_clusters) / (n_clusters - 1))
  
  cat("\nCalinski-Harabasz Index Calculation Complete:\n")
  cat("Calinski-Harabasz Index:", chi_score, "\n")
  
  
  return(chi_score)
}


#' Davies-Bouldin Index
#'
#' @param clustered.data Clustered data.
#' @param label_col Column for label.
#' @return Numeric DB index value.
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(30), y = rnorm(30), group = sample(1:3, 30, replace = TRUE))
#'   db <- davies_bouldin(df)
#'   cat("Davies-Bouldin Index:", db, "\n")
#' }
#' @export
davies_bouldin <- function(clustered.data, label_col = "group") {
  cat("Calculating Davies-Bouldin Index...\n")
  
  X <- clustered.data[, !names(clustered.data) %in% label_col, drop = FALSE]
  labels <- clustered.data[[label_col]]
  
  
  n_clusters <- length(unique(labels))
  
  if (n_clusters == 1) {
    cat("Warning: All samples belong to a single cluster. Returning 0.\n")
    return(0)
  }
  
  cluster_centers <- matrix(NA, nrow = n_clusters, ncol = ncol(X))
  cluster_distances <- numeric(n_clusters)
  
  label_to_index <- match(labels, unique(labels))
  cat("Cluster labels mapped to continuous indices.\n")
  
  for (k in 1:n_clusters) {
    cluster_points <- X[label_to_index == k, , drop = FALSE]
    if (nrow(cluster_points) == 0) {
      cat(sprintf("Warning: Cluster %d is empty. Skipping.\n", k))
      cluster_centers[k, ] <- NA
      cluster_distances[k] <- 0
      next
    }
    
    center_k <- colMeans(cluster_points)
    cluster_centers[k, ] <- center_k
    
    S_k <- mean(sqrt(rowSums((cluster_points - center_k)^2)))
    cluster_distances[k] <- S_k
  }
  cat("Cluster centers and within-cluster distances calculated.\n")
  
  valid_clusters <- !is.na(cluster_centers[, 1])
  cluster_centers <- cluster_centers[valid_clusters, , drop = FALSE]
  cluster_distances <- cluster_distances[valid_clusters]
  n_clusters <- sum(valid_clusters)
  
  if (n_clusters <= 1) {
    cat("Warning: Only one valid cluster found. Returning 0.\n")
    return(0)
  }
  
  center_distances <- as.matrix(dist(cluster_centers))
  diag(center_distances) <- Inf  
  cat("Between-cluster distances calculated.\n")
  
  db_values <- numeric(n_clusters)
  
  for (k in 1:n_clusters) {
    ratios <- (cluster_distances[k] + cluster_distances) / center_distances[k, ]
    db_values[k] <- max(ratios)
  }
  
  dbi_score <- mean(db_values)
  
  cat("\nDavies-Bouldin Index Calculation Complete:\n")
  cat("Davies-Bouldin Index:", dbi_score, "\n")
  
  return(dbi_score)
}



#' Silhouette Score
#'
#' @param clustered.data Clustered data.
#' @param label_col Column for label.
#' @return Numeric average silhouette score.
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(30), y = rnorm(30), group = sample(1:3, 30, replace = TRUE))
#'   sil <- silhouette_score(df)
#'   cat("Silhouette Score:", sil, "\n")
#' }
#' @export
silhouette_score <- function(clustered.data, label_col = "group") {
  cat("Calculating Silhouette Score...\n")
  
  X <- clustered.data[, !names(clustered.data) %in% label_col, drop = FALSE]
  labels <- clustered.data[[label_col]]
  
  n_clusters <- length(unique(labels))
  
  if (n_clusters < 2) {
    cat("Warning: All samples belong to a single cluster. Returning 0.\n")
    return(0)
  }
  
  int_labels <- as.integer(factor(labels))
  sil <- cluster::silhouette(int_labels, dist(X))
  overall_silhouette <- mean(sil[, 3], na.rm = TRUE)
  
  cat("\nSilhouette Score Calculation Complete:\n")
  cat("Silhouette Score:", overall_silhouette, "\n")
  
  return(overall_silhouette)
}


#' Sub Evaluation Results
#'
#' @param object Subtyping object.
#' @param seed Seed.
#' @return The updated \code{Subtyping} object (if input was Subtyping),
#'   or a named list of evaluation metrics.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object with clustered.data
#'   obj <- Sub_evaluation_results(obj, seed = 123)
#'   print(obj@evaluation_results)
#' }
#' @export
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
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' sub_obj <- Sub_normalize_process(sub_obj, normalize_method = "min_max")
#' sub_obj <- Sub_kmeans_with_optimal_k(sub_obj, k.max = 4, save_plots = FALSE, seed = 1)
#'
#' sub_obj@clustered.data$group <- sub_obj@info.data$cluster_kmeans
#' eval_res <- Sub_evaluation_results(sub_obj, seed = 1)
#' eval_res@evaluation_results
#' }
Sub_evaluation_results <- function(object ,
                                   seed = 123) {
  
  set.seed(seed)
  if (inherits(object, "Subtyping")) {
    clustered.data <- object@clustered.data
  } else if (is.data.frame(object)) {
    clustered.data <- object
  } else {
    stop("Input must be a Subtyping object or a data frame.")
  }
  set.seed(seed)
  
  if (ncol(clustered.data) < 2) {
    stop("clustered.data must contain features and one cluster label column.")
  }
  
  # Convert cluster labels to consecutive integers regardless of format
  # (handles "S1"/"S2", "1"/"2", factor levels, etc.)
  raw_labels <- clustered.data[, ncol(clustered.data)]
  int_labels  <- as.integer(factor(as.character(raw_labels)))
  clustered.data[, ncol(clustered.data)] <- int_labels

  labels     <- int_labels
  n_clusters <- length(unique(labels))
  n_samples  <- nrow(clustered.data)
  
  if (n_clusters < 2) {
    warning("Only one cluster detected. All evaluation metrics will be set to 0.")
  }
  
  label_col <- colnames(clustered.data)[ncol(clustered.data)]
  dbi_score <- davies_bouldin(clustered.data, label_col = label_col)
  chi_score <- calinski_harabasz(clustered.data, label_col = label_col)
  sil_score <- silhouette_score(clustered.data, label_col = label_col)
  
  evaluation_results <- list(
    Calinski_Harabasz = chi_score,
    Davies_Bouldin    = dbi_score,
    Silhouette        = sil_score,
    n_clusters        = n_clusters,
    n_samples         = n_samples
  )
  
  if (inherits(object, "Subtyping")) {
    object@evaluation_results <- evaluation_results
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'evaluation_results' slot updated.\n")
    return(object)
    
  } else {
    return(evaluation_results)
  }
}
