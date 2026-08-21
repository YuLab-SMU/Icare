#' Compare clustering results from info.data
#'
#' @param object Subtyping object.
#' @param methods Character vector of clustering column names in info.data (default: "cluster_kmeans", "cluster_lpa", "cluster_nmf").
#' @param output Either "matrix" (ARI matrix) or "data.frame" (sample x method grouping) or "list" (both).
#' @return Depends on output.
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' demo_df <- data.frame(
#'   id    = paste0("S", 1:60),
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2),
#'   feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#' )
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' sub_obj <- Sub_normalize_process(sub_obj, normalize_method = "min_max")
#'
#' sub_obj <- Sub_kmeans_with_optimal_k(sub_obj, k.max = 3, save_plots = FALSE, seed = 1)
#' sub_obj <- Sub_lpa_with_optimal_k(sub_obj, max_clusters = 3, verbose = FALSE,
#'                                   save_plots = FALSE, seed = 1)
#'
#' sub_obj@info.data$cluster_kmeans <- paste0("S", sub_obj@info.data$cluster_kmeans)
#' sub_obj@info.data$cluster_lpa    <- paste0("S", sub_obj@info.data$cluster_lpa)
#'
#' ari_mat <- compare_clusterings(sub_obj, methods = c("cluster_kmeans", "cluster_lpa"),
#'                                output = "matrix")
#' ari_mat
#' }
compare_clusterings <- function(object, 
                                methods = c("cluster_kmeans", "cluster_lpa", "cluster_nmf"),
                                output = "matrix") {
  if (is.null(object@info.data) || nrow(object@info.data) == 0) {
    stop("info.data is empty. Please run clustering methods first.")
  }
  .require_pkgs("mclust")
  present <- methods[methods %in% colnames(object@info.data)]
  if (length(present) == 0) stop("None of the specified clustering columns found in info.data.")
  
  group_df <- object@info.data[, present, drop = FALSE]
  group_df <- na.omit(group_df) 
  if (nrow(group_df) == 0) stop("No samples have complete clustering results.")
  
  if (output == "data.frame") {
    return(group_df)
  }
  
  n_methods <- length(present)
  if (n_methods < 2) {
    stop("Need at least 2 clustering methods to compare (found ", n_methods, 
         "). Provide more columns via 'methods', or use output = 'data.frame' instead.")
  }
  ari_mat <- matrix(1, nrow = n_methods, ncol = n_methods)
  rownames(ari_mat) <- colnames(ari_mat) <- present
  for (i in 1:(n_methods-1)) {
    for (j in (i+1):n_methods) {
      ari <- mclust::adjustedRandIndex(group_df[[i]], group_df[[j]])
      ari_mat[i, j] <- ari_mat[j, i] <- ari
    }
  }
  
  if (output == "matrix") return(ari_mat)
  if (output == "list") return(list(ari_matrix = ari_mat, group_table = group_df))
}

#' Plot Comparison Heatmap of Clustering Results
#'
#' Generates a heatmap of Adjusted Rand Index (ARI) values comparing multiple
#' clustering methods. The heatmap displays pairwise similarity between
#' clustering solutions, with numerical values overlaid on each cell.
#'
#' @param object A \code{Subtyping} S4 object containing clustering results
#'   in the \code{info.data} slot.
#' @param methods Character vector of column names in \code{object@info.data}
#'   that contain clustering assignments. Default is
#'   \code{c("cluster_kmeans", "cluster_lpa", "cluster_nmf")}.
#' @param save_dir Character string specifying the directory to save the plot.
#'   If \code{NULL}, the plot is not saved. Default is \code{NULL}.
#' @param width Numeric. Width of the saved PDF in inches. Default is \code{5}.
#' @param height Numeric. Height of the saved PDF in inches. Default is \code{4.5}.
#' @param base_size Numeric. Base font size for the heatmap. Default is \code{13}.
#' @param ... Additional arguments passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return Invisibly returns the \code{pheatmap} object. The heatmap is
#'   drawn on the current graphics device and optionally saved to PDF.
#'
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom pheatmap pheatmap
#' @export
#'
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
#' sub_obj <- Sub_kmeans_with_optimal_k(sub_obj, k.max = 3, save_plots = FALSE, seed = 1)
#' sub_obj <- Sub_lpa_with_optimal_k(sub_obj, max_clusters = 3, verbose = FALSE,
#'                                   save_plots = FALSE, seed = 1)
#' sub_obj@info.data$cluster_kmeans <- paste0("S", sub_obj@info.data$cluster_kmeans)
#' sub_obj@info.data$cluster_lpa    <- paste0("S", sub_obj@info.data$cluster_lpa)
#'
#' # save_dir omitted / NULL so the example does not write to disk.
#' plot_clustering_comparison(
#'   sub_obj, methods = c("cluster_kmeans", "cluster_lpa")
#' )
#' }
plot_clustering_comparison <- function(object,
                                       methods   = c("cluster_kmeans", "cluster_lpa", "cluster_nmf"),
                                       save_dir  = NULL,
                                       width     = 5,
                                       height    = 4.5,
                                       base_size = 13,
                                       ...) {
  requireNamespace('mclust',quietly = T)
  requireNamespace('NMF',quietly = T)
  ari_mat <- compare_clusterings(object, methods, output = "matrix")
  n_breaks <- 101
  breaks   <- seq(0, 1, length.out = n_breaks)
  colors   <- colorRampPalette(c("#F7F7F7", "#4393C3", "#08306B"))(n_breaks - 1)
  num_mat  <- matrix(sprintf("%.3f", ari_mat),
                     nrow = nrow(ari_mat),
                     dimnames = dimnames(ari_mat))
  diag(num_mat) <- "1.000"
  clean_names <- gsub("^cluster_", "", rownames(ari_mat))
  clean_names <- toupper(clean_names)
  rownames(ari_mat) <- colnames(ari_mat) <- clean_names
  rownames(num_mat) <- colnames(num_mat) <- clean_names

  p <- pheatmap::pheatmap(
    ari_mat,
    color            = colors,
    breaks           = breaks,
    display_numbers  = num_mat,
    number_color     = "black",
    fontsize_number  = base_size - 1,
    cluster_rows     = FALSE,
    cluster_cols     = FALSE,
    border_color     = "white",
    fontsize         = base_size,
    main             = "Adjusted Rand Index between clustering methods",
    ...
  )

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    pdf_path <- file.path(save_dir, "clustering_comparison_ARI.pdf")
    grDevices::pdf(pdf_path, width = width, height = height)
    print(p)
    grDevices::dev.off()
    cat("ARI heatmap saved to:", pdf_path, "\n")
  }

  invisible(p)
}
