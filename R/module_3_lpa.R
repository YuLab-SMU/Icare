#' LPA with Optimal K
#'
#' @param data Data.
#' @param max_clusters Max clusters.
#' @param model_names Model names.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param seed Seed.
#' @param verbose Verbose.
#' @importFrom mclust Mclust plot.mclustBIC
#' @importFrom grDevices recordPlot replayPlot
#' @export
lpa_with_optimal_k <- function(data,
                               max_clusters = 5,
                               model_names = c("EII", "VII", "EEI", "VEI", "EVI", "VVI"),
                               save_plots = TRUE,
                               save_dir = here('Subtyping',"lpa_result"),
                               plot_width = 7,
                               plot_height = 5,
                               base_size = 14,
                               seed = 123,
                               verbose = TRUE) {
  
  
  if (!all(sapply(data, is.numeric))) {
    non_numeric <- names(which(!sapply(data, is.numeric)))
    stop("All columns must be numeric. Non-numeric columns detected: ", 
         paste(non_numeric, collapse = ", "))
  }
  
  if (nrow(data) < 2) {
    stop("Data must contain at least 2 observations")
  }
  
  if (save_plots && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    if (verbose) cat("Created output directory:", save_dir, "\n")
  }
  
  set.seed(seed)
  if (verbose) cat("\nRunning LPA with model selection...\n")
  lpa_result <- mclust::Mclust(
    data,
    G = 1:max_clusters,
    modelNames = model_names,
    verbose = FALSE
  )
  
  optimal_k <- lpa_result$G
  optimal_model <- lpa_result$modelName
  cluster_labels <- lpa_result$classification
  
  if (verbose) {
    cat("Optimal number of clusters (K):", optimal_k, "\n")
    cat("Selected covariance model:", optimal_model, "\n")
  }
  
  bic_plot <- mclust::plot.mclustBIC(lpa_result$BIC, ylab = "BIC", 
                                     main = "LPA Model Selection by BIC")
  bic_plot <- recordPlot()  
  
  if (save_plots) {
    pdf(file.path(save_dir, "lpa_bic_plot.pdf"), 
        width = plot_width, height = plot_height)
    replayPlot(bic_plot)
    dev.off()
    
    if (verbose) cat("BIC plot saved to:", file.path(save_dir, "lpa_bic_plot.pdf"), "\n")
  }
  
  clustered_data <- data.frame(data, cluster = factor(cluster_labels))
  
  if (save_plots) {
    write.csv(clustered_data, 
              file.path(save_dir, "lpa_cluster_assignments.csv"),
              row.names = FALSE)
    if (verbose) cat("Cluster assignments saved to:", 
                     file.path(save_dir, "lpa_cluster_assignments.csv"), "\n")
  }
  
  return(list(
    clustered_data = clustered_data,
    cluster_labels = cluster_labels,
    optimal_model = optimal_model,
    optimal_k = optimal_k,
    lpa_object = lpa_result,
    plot_bic = bic_plot
  ))
}


#' Sub LPA with Optimal K
#'
#' @param object Subtyping object.
#' @param use_scaled_data Logical.
#' @param max_clusters Max clusters.
#' @param model_names Model names.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param seed Seed.
#' @param verbose Verbose.
#' @export
Sub_lpa_with_optimal_k <- function(object,
                                   use_scaled_data = TRUE,
                                   max_clusters = 5,
                                   model_names = c("EII", "VII", "EEI", "VEI", "EVI", "VVI"),
                                   save_plots = TRUE,
                                   save_dir = here('Subtyping',"cluster_results","lpa_result"),
                                   plot_width = 7,
                                   plot_height = 5,
                                   base_size = 14,
                                   seed = 123,
                                   verbose = TRUE) {
  
  if (inherits(object, "Subtyping")) {
    data <- if (use_scaled_data) {
      slot(object, "scale.data")
    } else {
      slot(object, "clean.data")
    }
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
  if (ncol(numeric_data) == 0) {
    stop("No numeric columns available for clustering")
  }
  
  if (verbose) cat("Starting LPA clustering analysis...\n")
  lpa_result <- lpa_with_optimal_k(
    numeric_data,
    max_clusters = max_clusters,
    model_names = model_names,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    seed = seed,
    verbose = verbose
  )
  
  clustered_data <- lpa_result$clustered_data
  colnames(clustered_data)[colnames(clustered_data) == "cluster"] <- "group"
  
  if (inherits(object, "Subtyping")) {
    object@cluster.results[["lpa.result"]] <- lpa_result
    object@Optimal.cluster <- lpa_result$optimal_k
    object@clustered.data <- clustered_data
    
    if (verbose) {
      cat("\nUpdating 'Subtyping' object...\n")
      cat("The 'Subtyping' object has been updated with:\n")
      cat("- 'cluster.results' slot (lpa.result)\n")
      cat("- 'Optimal.cluster' slot:", lpa_result$optimal_k, "\n")
      cat("- 'clustered.data' slot with group assignments\n")
      cat("Selected covariance model:", lpa_result$optimal_model, "\n")
    }
    
    return(object)
  } else {
    return(lpa_result)
  }
}
