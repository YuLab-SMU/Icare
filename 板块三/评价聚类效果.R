calinski_harabasz <- function(clustered.data) {
  X <- clustered.data[, -ncol(clustered.data), drop = FALSE]
  labels <- clustered.data[, ncol(clustered.data)]
  
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
      cluster_centers[k, ] <- global_center  # 空簇用全局中心填充
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


davies_bouldin <- function(clustered.data) {
  cat("Calculating Davies-Bouldin Index...\n")
  
  # 提取数据和标签
  X <- clustered.data[, -ncol(clustered.data), drop = FALSE]
  labels <- clustered.data[, ncol(clustered.data)]
  
 
  n_clusters <- length(unique(labels))
  
  if (n_clusters == 1) {
    cat("Warning: All samples belong to a single cluster. Returning 0.\n")
    return(0)
  }
  
  # 初始化簇中心和平均簇内距离
  cluster_centers <- matrix(NA, nrow = n_clusters, ncol = ncol(X))
  cluster_distances <- numeric(n_clusters)
  
  # 将标签映射到 1:n_clusters 的连续索引
  label_to_index <- match(labels, unique(labels))
  cat("Cluster labels mapped to continuous indices.\n")
  
  # 计算每个簇的中心和平均簇内距离
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
    
    # 计算平均簇内距离 S_k
    S_k <- mean(sqrt(rowSums((cluster_points - center_k)^2)))
    cluster_distances[k] <- S_k
  }
  cat("Cluster centers and within-cluster distances calculated.\n")
  
  # 移除空簇（如果有）
  valid_clusters <- !is.na(cluster_centers[, 1])
  cluster_centers <- cluster_centers[valid_clusters, , drop = FALSE]
  cluster_distances <- cluster_distances[valid_clusters]
  n_clusters <- sum(valid_clusters)
  
  if (n_clusters <= 1) {
    cat("Warning: Only one valid cluster found. Returning 0.\n")
    return(0)
  }
  
  # 计算簇间距离矩阵
  center_distances <- as.matrix(dist(cluster_centers))
  diag(center_distances) <- Inf  # 避免除以零
  cat("Between-cluster distances calculated.\n")
  
  # 计算每个簇的最大 (S_k + S_j) / d(c_k, c_j)
  db_values <- numeric(n_clusters)
  
  for (k in 1:n_clusters) {
    ratios <- (cluster_distances[k] + cluster_distances) / center_distances[k, ]
    db_values[k] <- max(ratios)
  }
  
  dbi_score <- mean(db_values)
  
  # 打印最终结果
  cat("\nDavies-Bouldin Index Calculation Complete:\n")
  cat("Davies-Bouldin Index:", dbi_score, "\n")

  return(dbi_score)
}



silhouette_score <- function(clustered.data) {
  # 打印开始信息
  cat("Calculating Silhouette Score...\n")
  
  X <- clustered.data[, -ncol(clustered.data), drop = FALSE]
  labels <- clustered.data[, ncol(clustered.data)]
  
  
  n_samples <- nrow(X)
  unique_labels <- unique(labels)
  n_clusters <- length(unique_labels)
  
  # 边界情况处理
  if (n_clusters == 1) {
    cat("Warning: All samples belong to a single cluster. Returning 0.\n")
    return(0)
  }
  
  # 计算距离矩阵（欧氏距离）
  cat("Calculating distance matrix...\n")
  dist_matrix <- dist(X)
  
  # 将标签映射到连续整数（从1开始）
  label_to_index <- match(labels, unique_labels)
  cat("Cluster labels mapped to continuous indices.\n")
  
  # 计算每个样本的轮廓系数
  silhouette_values <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    current_cluster <- label_to_index[i]
    
    # 计算a(i)：样本i到同簇其他样本的平均距离
    same_cluster <- which(label_to_index == current_cluster)
    a_i <- mean(as.matrix(dist_matrix)[i, same_cluster[same_cluster != i]])
    
    # 计算b(i)：样本i到其他最近簇的平均距离
    other_clusters <- setdiff(unique(label_to_index), current_cluster)
    b_i_values <- sapply(other_clusters, function(k) {
      mean(as.matrix(dist_matrix)[i, label_to_index == k])
    })
    b_i <- min(b_i_values)
    
    # 计算样本i的轮廓系数
    silhouette_values[i] <- (b_i - a_i) / max(a_i, b_i)
  }
  
  overall_silhouette <- mean(silhouette_values)
  
  # 打印最终结果
  cat("\nSilhouette Score Calculation Complete:\n")
  cat("Silhouette Score:", overall_silhouette, "\n")

  
  return(overall_silhouette)
}


Sub_evaluation_results <- function(object ,
                                   seed = 123) {
  
  set.seed(seed)
  
  if (inherits(object, "Subtyping")) {
    clustered.data <- object@clustered.data
  } else if (is.list(object)) {
    clustered.data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame.")
  }
  set.seed(seed)
  
  if (ncol(clustered.data) < 2) {
    stop("clustered.data must contain features and one cluster label column.")
  }
  
  clustered.data[, ncol(clustered.data)] <-
    as.numeric(as.character(clustered.data[, ncol(clustered.data)]))
  
  labels <- clustered.data[, ncol(clustered.data)]
  n_clusters <- length(unique(labels))
  n_samples  <- nrow(clustered.data)
  
  if (n_clusters < 2) {
    warning("Only one cluster detected. All evaluation metrics will be set to 0.")
  }
  
  dbi_score <- davies_bouldin(clustered.data)
  chi_score <- calinski_harabasz(clustered.data)
  sil_score <- silhouette_score(clustered.data)
  
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
