# =============================================================================
# K-means Subtyping Module
# =============================================================================
# Required: factoextra, NbClust, cluster, ggplot2, wesanderson, tibble, dplyr

# =============================================================================
# SECTION 1: Core clustering helpers
# =============================================================================

#' K-means Clustering (single K)
#'
#' @param data   Numeric data frame or matrix.
#' @param k      Number of clusters.
#' @param nstart Number of random starts.
#' @return Integer vector of cluster assignments.
#' @export
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(50), y = rnorm(50))
#'   clusters <- kmeans_clustering(df, k = 3)
#'   table(clusters)
#' }
kmeans_clustering <- function(data, k, nstart = 50) {
  kmeans(data, centers = k, nstart = nstart)$cluster
}

# =============================================================================
# SECTION 2: Main function - optimal K selection + full visualisation
# =============================================================================

#' K-means with Optimal K
#'
#' Determines the optimal number of clusters using the silhouette index
#' (via \code{NbClust}) and produces four publication-ready plots:
#' \enumerate{
#'   \item Silhouette-width curve (elbow-style)
#'   \item Cluster scatter (PCA projection, \code{factoextra::fviz_cluster})
#'   \item Silhouette bar chart per sample
#'   \item WSS elbow curve
#' }
#'
#' @param data         Numeric data frame.
#' @param palette_name wesanderson palette name for cluster colours.
#' @param save_plots   Whether to save plots to disk.
#' @param save_dir     Output directory.
#' @param plot_width   Plot width (inches).
#' @param plot_height  Plot height (inches).
#' @param base_size    Base font size.
#' @param nstart       Number of K-means random starts.
#' @param k.max        Maximum K to search.
#' @param seed         Random seed.
#' @importFrom factoextra fviz_nbclust fviz_cluster fviz_silhouette
#' @importFrom NbClust NbClust
#' @importFrom cluster silhouette
#' @importFrom wesanderson wes_palette
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline labs
#'   scale_x_continuous ggsave theme element_blank element_text
#' @export
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(100), y = rnorm(100))
#'   result <- kmeans_with_optimal_k(df, palette_name = "Zissou1", save_plots = FALSE)
#'   print(result$optimal_k)
#' }
kmeans_with_optimal_k <- function(data,
                                  palette_name = "Zissou1",
                                  save_plots   = TRUE,
                                  save_dir     = file.path(get_output_dir("Subtyping", "cluster_results"),
                                                             "kmeans_result"),
                                  plot_width   = 6,
                                  plot_height  = 5,
                                  base_size    = 14,
                                  nstart       = 50,
                                  k.max        = 15,
                                  seed         = 123) {

  set.seed(seed)

  if (!is.data.frame(data) || !all(sapply(data, is.numeric)))
    stop("Input data must be a numeric data frame.")

  if (save_plots && !dir.exists(save_dir))
    dir.create(save_dir, recursive = TRUE)

  pub_theme <- .pub_theme(base_size)

  # ---- 1. Determine optimal K via NbClust (silhouette) ----
  cat("Determining optimal K via NbClust silhouette index...\n")
  nc        <- NbClust::NbClust(data, min.nc = 2, max.nc = k.max,
                                method = "kmeans", index = "silhouette")
  optimal_k <- as.integer(nc$Best.nc[1])
  cat("Optimal K:", optimal_k, "\n")

  # ---- 2. Silhouette-width curve ----
  main_col <- tryCatch(
    wesanderson::wes_palette(palette_name, 1, type = "continuous"),
    error = function(e) "#3A7DC9"
  )
  sil_curve_plt <- factoextra::fviz_nbclust(
      data, kmeans, method = "silhouette", k.max = k.max
    ) +
    ggplot2::geom_vline(xintercept = optimal_k, linetype = "dashed",
                        colour = "#E84040", linewidth = 0.8) +
    ggplot2::labs(title    = "Optimal K \u2013 Silhouette Width",
                  subtitle = paste0("Recommended K = ", optimal_k),
                  x        = "Number of Clusters (K)",
                  y        = "Average Silhouette Width") +
    pub_theme

  if (save_plots) {
    ggplot2::ggsave(file.path(save_dir, "kmeans_silhouette_curve.pdf"),
                    sil_curve_plt, width = plot_width, height = plot_height, dpi = 300)
    ggplot2::ggsave(file.path(save_dir, "kmeans_silhouette_curve.png"),
                    sil_curve_plt, width = plot_width, height = plot_height, dpi = 150)
    cat("- Silhouette curve saved.\n")
  }

  # ---- 3. WSS elbow curve ----
  wss_plt <- factoextra::fviz_nbclust(
      data, kmeans, method = "wss", k.max = k.max
    ) +
    ggplot2::geom_vline(xintercept = optimal_k, linetype = "dashed",
                        colour = "#E84040", linewidth = 0.8) +
    ggplot2::labs(title    = "Within-cluster Sum of Squares (Elbow)",
                  subtitle = paste0("Recommended K = ", optimal_k),
                  x        = "Number of Clusters (K)",
                  y        = "Total Within-cluster SS") +
    pub_theme

  if (save_plots) {
    ggplot2::ggsave(file.path(save_dir, "kmeans_wss_elbow.pdf"),
                    wss_plt, width = plot_width, height = plot_height, dpi = 300)
    ggplot2::ggsave(file.path(save_dir, "kmeans_wss_elbow.png"),
                    wss_plt, width = plot_width, height = plot_height, dpi = 150)
    cat("- WSS elbow curve saved.\n")
  }

  # ---- 4. Final K-means fit ----
  km             <- kmeans(data, centers = optimal_k, nstart = nstart)
  cluster_labels <- km$cluster
  clustered_data <- cbind(data, group = factor(cluster_labels))

  pal <- tryCatch(
    wesanderson::wes_palette(palette_name, optimal_k, type = "continuous"),
    error = function(e) scales::hue_pal()(optimal_k)
  )

  # ---- 5. Cluster scatter (PCA projection) ----
  cluster_plt <- factoextra::fviz_cluster(
      km, data = data, geom = "point", palette = pal,
      ggtheme = pub_theme, main = paste0("K-means Clustering  (K = ", optimal_k, ")")
    ) +
    pub_theme

  if (save_plots) {
    ggplot2::ggsave(file.path(save_dir, "kmeans_cluster_scatter.pdf"),
                    cluster_plt, width = plot_width * 1.2, height = plot_height, dpi = 300)
    ggplot2::ggsave(file.path(save_dir, "kmeans_cluster_scatter.png"),
                    cluster_plt, width = plot_width * 1.2, height = plot_height, dpi = 150)
    cat("- Cluster scatter saved.\n")
  }

  # ---- 6. Silhouette bar ----
  sil_obj     <- cluster::silhouette(km$cluster, stats::dist(data))
  avg_sil     <- mean(sil_obj[, 3])
  sil_bar_plt <- factoextra::fviz_silhouette(sil_obj, palette = pal,
                                             ggtheme = pub_theme) +
    ggplot2::labs(title    = paste0("Silhouette Plot  (K = ", optimal_k, ")"),
                  subtitle = paste0("Average silhouette width = ", round(avg_sil, 3)),
                  x        = "Samples",
                  y        = "Silhouette Width") +
    pub_theme +
    ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  if (save_plots) {
    ggplot2::ggsave(file.path(save_dir, "kmeans_silhouette_bar.pdf"),
                    sil_bar_plt, width = plot_width, height = plot_height, dpi = 300)
    ggplot2::ggsave(file.path(save_dir, "kmeans_silhouette_bar.png"),
                    sil_bar_plt, width = plot_width, height = plot_height, dpi = 150)
    cat("- Silhouette bar saved.\n")
  }

  return(list(
    clustered_data = clustered_data,
    cluster_labels = cluster_labels,
    optimal_k      = optimal_k,
    model          = km,
    plots = list(
      silhouette_curve = sil_curve_plt,
      wss_elbow        = wss_plt,
      cluster_scatter  = cluster_plt,
      silhouette_bar   = sil_bar_plt
    )
  ))
}

# =============================================================================
# SECTION 3: Subtyping-object wrapper
# =============================================================================

#' Sub K-means with Optimal K
#'
#' Wraps \code{kmeans_with_optimal_k} for use with a \code{Subtyping} object.
#'
#' @param object        Subtyping object (or numeric data frame).
#' @param use_scaled_data Use the \code{scale.data} slot when \code{TRUE}.
#' @param palette_name  wesanderson palette name for cluster colours.
#' @param save_plots    Whether to save plots.
#' @param save_dir      Output directory.
#' @param plot_width    Plot width (inches).
#' @param plot_height   Plot height (inches).
#' @param base_size     Base font size.
#' @param nstart        Number of K-means random starts.
#' @param k.max         Maximum K to search.
#' @param seed          Random seed.
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
#' sub_obj <- Sub_kmeans_with_optimal_k(
#'   sub_obj, use_scaled_data = TRUE, k.max = 4,
#'   save_plots = FALSE, seed = 1
#' )
#' table(sub_obj@info.data$cluster_kmeans)
#' }
Sub_kmeans_with_optimal_k <- function(object,
                                      use_scaled_data = TRUE,
                                      palette_name    = "Zissou1",
                                      save_plots      = TRUE,
                                      save_dir        = file.path(get_output_dir("Subtyping", "cluster_results"),
                                                                   "kmeans_result"),
                                      plot_width      = 6,
                                      plot_height     = 5,
                                      base_size       = 14,
                                      nstart          = 50,
                                      k.max           = 15,
                                      seed            = 123) {

  if (inherits(object, "Subtyping")) {
    data     <- if (use_scaled_data) slot(object, "scale.data") else slot(object, "clean.data")
    raw_data <- slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- raw_data <- object
  } else {
    stop("Input must be a 'Subtyping' object or a data frame.")
  }

  if (is.null(data) || nrow(data) == 0) stop("No valid data found in the input.")

  numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
  if (ncol(numeric_data) == 0) stop("No numeric columns found.")

  cat("Starting K-means clustering analysis...\n")

  kmeans_result <- kmeans_with_optimal_k(
    numeric_data,
    palette_name = palette_name,
    save_plots   = save_plots,
    save_dir     = save_dir,
    plot_width   = plot_width,
    plot_height  = plot_height,
    base_size    = base_size,
    nstart       = nstart,
    k.max        = k.max,
    seed         = seed
  )

  # Merge group column back onto raw_data (preserves row names)
  raw_rn  <- tibble::rownames_to_column(as.data.frame(raw_data),  "sample_id")
  grp_rn  <- dplyr::select(
    tibble::rownames_to_column(as.data.frame(kmeans_result$clustered_data), "sample_id"),
    sample_id, group
  )
  combined_data <- tibble::column_to_rownames(
    dplyr::left_join(raw_rn, grp_rn, by = "sample_id"),
    "sample_id"
  )

  # Write cluster_kmeans into info.data
  cluster_df <- data.frame(
    sample  = rownames(kmeans_result$clustered_data),
    cluster = kmeans_result$cluster_labels,
    stringsAsFactors = FALSE
  )

  if (inherits(object, "Subtyping")) {
    if (nrow(object@info.data) == 0)
      object@info.data <- data.frame(row.names = rownames(object@clean.data))
    object@info.data <- object@info.data[rownames(object@clean.data), , drop = FALSE]
    idx <- match(rownames(object@info.data), cluster_df$sample)
    object@info.data$cluster_kmeans <- cluster_df$cluster[idx]

    object@cluster.results[["kmeans.result"]] <- kmeans_result
    object@Optimal.cluster <- kmeans_result$optimal_k
    object@clustered.data  <- combined_data

    cat("Updating 'Subtyping' object...\n")
    cat("- cluster.results[kmeans.result] updated\n")
    cat("- Optimal.cluster:", kmeans_result$optimal_k, "\n")
    cat("- clustered.data updated\n")
    cat("- info.data$cluster_kmeans written\n")
    return(object)
  } else {
    return(kmeans_result)
  }
}
