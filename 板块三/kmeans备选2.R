kmeans_with_optimal_k <- function(data,
                                  palette_name = "Zissou1",
                                  save_plots = TRUE,
                                  save_dir = here('Subtyping', "kmeans_result"),
                                  plot_width = 5,
                                  plot_height = 5,
                                  base_size = 14,
                                  nstart = 50,
                                  seed = 123,
                                  k.max = 15) {
  # 参数验证
  set.seed(seed)
  if (!is.data.frame(data) || !all(sapply(data, is.numeric))) {
    stop("Input data must be a numeric data frame.")
  }
  if (!palette_name %in% names(wesanderson::wes_palettes)) {
    stop(paste("Invalid palette name. Choose from:", 
               paste(names(wesanderson::wes_palettes), collapse = ", ")))
  }
  
  # 创建目录（一次性）
  if (save_plots && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    cat("Created directory:", save_dir, "\n")
  }
  
  # 1. 确定最佳K值
  silhouette <- fviz_nbclust(data, kmeans, method = "silhouette", k.max = k.max) +
    labs(title = "Silhouette Method for Optimal K") +
    scale_color_manual(values = wes_palette(palette_name)) +
    theme_minimal(base_size = base_size)
  
  optimal_k <- as.numeric(silhouette$data[which.max(silhouette$data$y), "clusters"])
  cat("Optimal K value is:", optimal_k, "\n")
  
  # 2. 执行聚类
  kmeans_result <- kmeans(data, centers = optimal_k, nstart = nstart)
  clustered_data <- cbind(data, group = factor(kmeans_result$cluster))
  cluster_labels<-kmeans_result$cluster
  # 3. 可视化
  cluster_plot <- fviz_cluster(kmeans_result, data = data, geom = "point") + 
    ggtitle(paste("K-means Clustering (K =", optimal_k, ")")) + 
    scale_color_manual(values = wes_palette(palette_name)) +
    theme_bw(base_size = base_size) +
    theme(axis.title = element_blank())
  
  # 打印结果
  print(silhouette+cluster_plot)
  
  # 保存图像
  if (save_plots) {
    plots <- list(
      silhouette = list(plot = silhouette, filename = "silhouette_plot.pdf"),
      clusters = list(plot = cluster_plot, filename = "cluster_plot.pdf")
    )
    
    for (p in plots) {
      ggsave(
        filename = file.path(save_dir, p$filename),
        plot = p$plot,
        width = ifelse(grepl("silhouette", p$filename), plot_width, plot_width * 1.2),
        height = plot_height,
        device = "pdf"
      )
      cat("Saved", p$filename, "to", save_dir, "\n")
    }
  }
  
  return(list(
    data = clustered_data,
    cluster_labels =cluster_labels,
    model = kmeans_result,
    optimal_k = optimal_k,
    plots = list(silhouette = silhouette, clusters = cluster_plot)
  ))
}


Sub_kmeans_with_optimal_k <- function(object,
                                      use_scaled_data = TRUE,  
                                      palette_name = "Zissou1",
                                      save_plots = TRUE,
                                      save_dir = here('Subtyping', "cluster_results", "kmeans_result"),
                                      plot_width = 5,
                                      plot_height = 5,
                                      base_size = 14,
                                      nstart = 50,
                                      seed = 123,
                                      k.max = 15) {
  
  # 参数验证
  if (!palette_name %in% names(wesanderson::wes_palettes)) {
    stop(paste("Invalid palette name. Choose from:", 
               paste(names(wesanderson::wes_palettes), collapse = ", ")))
  }
  
  # 数据提取
  if (inherits(object, "Subtyping")) {
    if (use_scaled_data) {
      data <- slot(object, "scale.data")
    } else {
      data <- slot(object, "clean.data")
    }
    raw_data <- slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- object
    raw_data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame")
  }
  
  # 数据验证
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
  
  if (ncol(numeric_data) == 0) {
    stop("No numeric columns found in the data. Please provide numeric data.")
  }
  
  cat("Starting K-means clustering analysis...\n")
  
  # 执行K-means聚类
  kmeans_result <- kmeans_with_optimal_k(
    numeric_data,
    palette_name = palette_name,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    nstart = nstart,
    seed = seed,
    k.max = k.max
  )
  
  # 数据整合
  clustered_data <- kmeans_result$data
  cluster_labels <- kmeans_result$cluster_labels
  
  # 保留原始行名
  raw_data_with_rownames <- raw_data %>% 
    tibble::rownames_to_column("sample_id")
  
  clustered_data_with_rownames <- clustered_data %>% 
    tibble::rownames_to_column("sample_id") %>% 
    dplyr::select(sample_id, group)
  
  combined_data <- raw_data_with_rownames %>% 
    dplyr::left_join(clustered_data_with_rownames, by = "sample_id") %>% 
    tibble::column_to_rownames("sample_id") 
  
  # 更新Subtyping对象或返回结果
  if (inherits(object, "Subtyping")) {
    object@cluster.results[["kmeans.result"]] <- kmeans_result
    object@Optimal.cluster <- kmeans_result$optimal_k
    object@clustered.data <- combined_data
    object@cluster.results[["cluster.labels"]] <- cluster_labels
    
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'cluster.results' slot updated with kmeans.result\n")
    cat("- 'Optimal.cluster' slot updated with optimal K value\n")
    cat("- 'clustered.data' slot updated with clustering results\n")
    cat("- 'cluster.results' slot updated with cluster labels\n")
    
    return(object)
  } else {
    return(list(
      clustered_data = combined_data,
      cluster_labels = cluster_labels,
      kmeans_result = kmeans_result
    ))
  }
}
