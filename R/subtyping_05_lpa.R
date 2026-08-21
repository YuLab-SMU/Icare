# =============================================================================
# LPA (Latent Profile Analysis / Gaussian Mixture Modelling) Subtyping Module
# =============================================================================
# Required: mclust, cluster, ggplot2, reshape2, wesanderson, scales, here
#
# Plots follow the mclust book:  https://mclust-org.github.io/mclust-book/
# Figures produced:
#   1. BIC model-selection plot  (multi-model, multi-K)
#   2. PCA scatter with cluster ellipses
#   3. Variable profile plot     (mean per cluster, per variable)
#   4. Classification probability heatmap
#   5. Silhouette bar chart

# =============================================================================
# SECTION 0: Shared helpers  (also sourced by kmeans / nmf modules)
# =============================================================================

#' @keywords internal
.pub_theme <- function(base_size = 14) {
  ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, size = base_size + 2),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = base_size - 1,
                                               color = "grey40"),
      axis.title       = ggplot2::element_text(size = base_size),
      axis.text        = ggplot2::element_text(size = base_size - 2),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(size = base_size - 1),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.key.size  = grid::unit(0.8, "lines"),
      legend.margin    = ggplot2::margin(t = 5)
    )
}

#' @keywords internal
.get_palette <- function(palette_name, n) {
  tryCatch(
    wesanderson::wes_palette(palette_name, n, type = "continuous"),
    error = function(e) {
      warning("wesanderson palette '", palette_name, "' failed - using hue_pal fallback.")
      scales::hue_pal()(n)
    }
  )
}

# =============================================================================
# SECTION 1: Core LPA function
# =============================================================================

#' LPA with Optimal K (core)
#'
#' Fits a Gaussian Mixture Model via \code{mclust::Mclust} over a range of K
#' and covariance-structure models, selects the BIC-optimal solution, and
#' produces five diagnostic plots following the mclust book.
#'
#' @param data             Numeric data frame (samples x variables).
#' @param max_clusters     Maximum number of mixture components (G) to evaluate.
#' @param model_names      mclust covariance model names to evaluate.
#' @param save_plots       Whether to save plots to \code{save_dir}.
#' @param save_dir         Output directory.
#' @param plot_width       Base plot width (inches).
#' @param plot_height      Base plot height (inches).
#' @param base_size        Base font size.
#' @param seed             Random seed.
#' @param verbose          Print progress messages.
#' @param show_silhouette  Whether to draw silhouette plot.
#' @param show_profile     Whether to draw variable-profile plot.
#' @param show_class_prob  Whether to draw classification-probability heatmap.
#' @param color_palette    wesanderson palette for cluster fill colours.
#' @importFrom cluster silhouette
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_bar geom_tile geom_errorbar labs scale_fill_manual scale_colour_manual ggsave
#' @importFrom reshape2 melt
#' @examples
#' \dontrun{
#'   df <- data.frame(x = rnorm(100), y = rnorm(100), z = rnorm(100))
#'   result <- lpa_with_optimal_k(df, max_clusters = 4, save_plots = FALSE)
#'   print(result$optimal_k)
#' }
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' require('mclust')
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
#' sub_obj <- Sub_lpa_with_optimal_k(
#'   sub_obj, use_scaled_data = TRUE, max_clusters = 3,
#'   verbose = FALSE, save_plots = FALSE, seed = 1
#' )
#' table(sub_obj@info.data$cluster_lpa)
#' }
lpa_with_optimal_k <- function(data,
                               max_clusters    = 5,
                               model_names     = c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
                                                   "EEE", "EEV", "VEV", "EVV",
                                                   "VVE", "VEE", "EVE", "VVV"),
                               save_plots      = TRUE,
                               save_dir        = NULL,
                               plot_width      = 7,
                               plot_height     = 5,
                               base_size       = 14,
                               seed            = 123,
                               verbose         = TRUE,
                               show_silhouette = TRUE,
                               show_profile    = TRUE,
                               show_class_prob = TRUE,
                               color_palette   = "Darjeeling1") {

  # ---- Validate ----
  if (!all(sapply(data, is.numeric))) {
    bad <- names(which(!sapply(data, is.numeric)))
    stop("All columns must be numeric. Non-numeric: ", paste(bad, collapse = ", "))
  }
  if (nrow(data) < 2) stop("Data must have at least 2 observations.")
  .require_pkgs("mclust")
  if (save_plots && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    if (verbose) cat("Created output directory:", save_dir, "\n")
  }

  pub_theme <- .pub_theme(base_size)
  get_cols  <- function(n) .get_palette(color_palette, n)

  # ---- Fit Mclust ----
  set.seed(seed)
  if (verbose) cat("Running LPA with model selection...\n")
  lpa_result <- mclust::Mclust(data,
                               G          = 1:max_clusters,
                               modelNames = model_names,
                               initialization = list(subset = NULL), 
                               verbose    = FALSE)

  optimal_k      <- lpa_result$G
  optimal_model  <- lpa_result$modelName
  cluster_labels <- lpa_result$classification

  if (verbose) {
    cat("Optimal number of clusters (K):", optimal_k, "\n")
    cat("Selected covariance model     :", optimal_model, "\n")
  }

  cols <- get_cols(optimal_k)

  # ====================================================================
  # Plot 1: BIC model-selection plot
  #   mclustBIC dims: rows = G values, cols = model names
  # ====================================================================
  bic_plot <- NULL
  bic_obj  <- lpa_result$BIC

  if (!is.null(bic_obj)) {
    bic_mat <- as.matrix(bic_obj)   # [G x models]

    if (is.matrix(bic_mat) && !all(is.na(bic_mat))) {
      g_vals    <- as.integer(rownames(bic_mat))
      model_vec <- colnames(bic_mat)

      DF <- data.frame(
        G     = rep(g_vals,    times = length(model_vec)),
        Model = rep(model_vec, each  = length(g_vals)),
        BIC   = as.vector(bic_mat),
        stringsAsFactors = FALSE
      )
      DF <- DF[!is.na(DF$BIC), ]

      canonical <- tryCatch(mclust::mclust.options("emModelNames"),
                            error = function(e) unique(DF$Model))
      DF$Model <- factor(DF$Model,
                         levels = canonical[canonical %in% unique(DF$Model)])

      n_models   <- length(levels(DF$Model))
      mod_cols   <- tryCatch(mclust::mclust.options("bicPlotColors")[seq_len(n_models)],
                             error = function(e) scales::hue_pal()(n_models))
      mod_shapes <- tryCatch(mclust::mclust.options("bicPlotSymbols")[seq_len(n_models)],
                             error = function(e) rep(16L, n_models))

      bic_plot <- ggplot2::ggplot(DF, ggplot2::aes(x = G, y = BIC,
                                                   colour = Model,
                                                   shape  = Model)) +
        ggplot2::geom_line(linewidth = 1, alpha = 0.85) +
        ggplot2::geom_point(size = 2.8, alpha = 0.95) +
        ggplot2::scale_colour_manual(values = mod_cols,   name = "Model") +
        ggplot2::scale_shape_manual( values = mod_shapes, name = "Model") +
        ggplot2::scale_x_continuous(breaks = seq_len(max_clusters)) +
        ggplot2::labs(
          title    = "Model Selection by BIC (all models)",
          subtitle = paste0("Optimal: K = ", optimal_k, ", Model = ", optimal_model),
          x        = "Number of Components (K)",
          y        = "Bayesian Information Criterion (BIC)"
        ) +
        pub_theme +
        ggplot2::guides(
          colour = ggplot2::guide_legend(ncol = 4, override.aes = list(size = 3)),
          shape  = ggplot2::guide_legend(ncol = 4, override.aes = list(size = 3))
        )

      if (save_plots) {
        ggplot2::ggsave(file.path(save_dir, "lpa_bic_plot.pdf"), bic_plot,
                        width = plot_width + 2, height = plot_height + 1, dpi = 300)
        ggplot2::ggsave(file.path(save_dir, "lpa_bic_plot.png"), bic_plot,
                        width = plot_width + 2, height = plot_height + 1, dpi = 150)
        if (verbose) cat("- BIC plot saved.\n")
      }

      # ---- Extra: best-BIC-per-K summary plot --------------------------------
      # For each K, take the model with the highest BIC. This makes it visually
      # unambiguous why a particular K was chosen.
      best_per_k <- do.call(rbind, lapply(split(DF, DF$G), function(d) {
        d[which.max(d$BIC), , drop = FALSE]
      }))
      best_per_k$Selected <- best_per_k$G == optimal_k

      bic_summary_plot <- ggplot2::ggplot(
          best_per_k,
          ggplot2::aes(x = G, y = BIC)
        ) +
        ggplot2::geom_line(linewidth = 1, colour = "grey50") +
        ggplot2::geom_point(
          ggplot2::aes(colour = Selected, size = Selected, shape = Selected)
        ) +
        ggplot2::scale_colour_manual(
          values = c("FALSE" = "grey40", "TRUE" = "#D6604D"),
          labels = c("FALSE" = "Other K", "TRUE" = paste0("Selected K = ", optimal_k)),
          name   = ""
        ) +
        ggplot2::scale_size_manual(
          values = c("FALSE" = 3, "TRUE" = 5), guide = "none"
        ) +
        ggplot2::scale_shape_manual(
          values = c("FALSE" = 16, "TRUE" = 18), guide = "none"
        ) +
        ggplot2::geom_vline(
          xintercept = optimal_k, linetype = "dashed",
          colour = "#D6604D", linewidth = 0.8
        ) +
        ggplot2::geom_text(
          data = best_per_k[best_per_k$Selected, ],
          ggplot2::aes(
            label = paste0("K=", G, "\n(", Model, ")"),
            y     = BIC
          ),
          vjust = -0.9, hjust = 0.5,
          size = base_size / 4.5, colour = "#D6604D"
        ) +
        ggplot2::scale_x_continuous(breaks = seq_len(max_clusters)) +
        ggplot2::labs(
          title    = "Optimal K Selection \u2013 Best BIC per K",
          subtitle = paste0(
            "Each point = best model at that K  |  ",
            "Selected: K = ", optimal_k, ", Model = ", optimal_model
          ),
          x = "Number of Components (K)",
          y = "Best BIC (across all covariance models)"
        ) +
        pub_theme +
        ggplot2::theme(legend.position = "top")

      if (save_plots) {
        ggplot2::ggsave(file.path(save_dir, "lpa_bic_summary_plot.pdf"),
                        bic_summary_plot,
                        width = plot_width + 1, height = plot_height, dpi = 300)
        ggplot2::ggsave(file.path(save_dir, "lpa_bic_summary_plot.png"),
                        bic_summary_plot,
                        width = plot_width + 1, height = plot_height, dpi = 150)
        if (verbose) cat("- BIC summary (best-per-K) plot saved.\n")
      }
      print(bic_summary_plot)
    }
  } else {
    warning("BIC object is NULL; BIC plot skipped.")
  }

  # ====================================================================
  # Plot 2: PCA scatter with cluster ellipses
  # ====================================================================
  pca_plot_obj <- NULL
  if (ncol(data) >= 2) {
    pca_res  <- prcomp(data, scale. = TRUE)
    var_expl <- summary(pca_res)$importance[2, 1:2] * 100
    pca_df   <- data.frame(
      PC1     = pca_res$x[, 1],
      PC2     = pca_res$x[, 2],
      Cluster = factor(cluster_labels)
    )

    pca_plot_obj <- ggplot2::ggplot(pca_df,
                                    ggplot2::aes(x = PC1, y = PC2,
                                                 colour = Cluster,
                                                 fill   = Cluster)) +
      ggplot2::geom_point(alpha = 0.75, size = 2.5, shape = 21, stroke = 0.5) +
      ggplot2::stat_ellipse(geom = "polygon", alpha = 0.10, level = 0.95,
                            type = "norm", linewidth = 0.8) +
      ggplot2::scale_colour_manual(values = cols, name = "Cluster") +
      ggplot2::scale_fill_manual(  values = cols, name = "Cluster") +
      ggplot2::labs(
        title    = "LPA Clustering \u2013 PCA Projection",
        subtitle = paste0("K = ", optimal_k, ", Model = ", optimal_model),
        x        = paste0("PC1 (", round(var_expl[1], 1), "% variance)"),
        y        = paste0("PC2 (", round(var_expl[2], 1), "% variance)")
      ) +
      pub_theme

    if (save_plots) {
      ggplot2::ggsave(file.path(save_dir, "lpa_pca_plot.pdf"), pca_plot_obj,
                      width = plot_width, height = plot_height, dpi = 300)
      ggplot2::ggsave(file.path(save_dir, "lpa_pca_plot.png"), pca_plot_obj,
                      width = plot_width, height = plot_height, dpi = 150)
      if (verbose) cat("- PCA scatter saved.\n")
    }
  }

  # ====================================================================
  # Plot 3: Variable profile plot  (mclust book 5.3)
  #   Mean +/- SE of each variable per cluster
  # ====================================================================
  profile_plot <- NULL
  if (show_profile) {
    # Compute cluster means and SE
    profile_long <- do.call(rbind, lapply(seq_len(optimal_k), function(k) {
      idx <- which(cluster_labels == k)
      sub <- data[idx, , drop = FALSE]
      data.frame(
        Cluster  = factor(k),
        Variable = colnames(sub),
        Mean     = colMeans(sub),
        SE       = apply(sub, 2, function(x) sd(x) / sqrt(length(x))),
        stringsAsFactors = FALSE
      )
    }))

    profile_plot <- ggplot2::ggplot(
        profile_long,
        ggplot2::aes(x = Variable, y = Mean, colour = Cluster, group = Cluster)
      ) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SE, ymax = Mean + SE),
                             width = 0.2, linewidth = 0.6, alpha = 0.7) +
      ggplot2::scale_colour_manual(values = cols, name = "Cluster") +
      ggplot2::labs(
        title    = "Variable Profile by Cluster",
        subtitle = paste0("K = ", optimal_k, " \u2013 Mean \u00b1 SE per variable"),
        x        = "Variable",
        y        = "Mean value"
      ) +
      pub_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                         size = base_size - 3))

    if (save_plots) {
      ggplot2::ggsave(file.path(save_dir, "lpa_profile_plot.pdf"), profile_plot,
                      width = max(plot_width, ncol(data) * 0.55 + 2),
                      height = plot_height, dpi = 300)
      ggplot2::ggsave(file.path(save_dir, "lpa_profile_plot.png"), profile_plot,
                      width = max(plot_width, ncol(data) * 0.55 + 2),
                      height = plot_height, dpi = 150)
      if (verbose) cat("- Variable profile plot saved.\n")
    }
  }

  # ====================================================================
  # Plot 4: Classification-probability heatmap  (mclust book 5.2)
  #   Rows = clusters, Cols = samples (sorted by assigned cluster)
  # ====================================================================
  class_prob_plot <- NULL
  if (show_class_prob && !is.null(lpa_result$z)) {
    z_mat  <- lpa_result$z           # n x K posterior probabilities
    order_idx <- order(cluster_labels)
    z_ord  <- z_mat[order_idx, , drop = FALSE]
    colnames(z_ord) <- paste0("C", seq_len(ncol(z_ord)))

    z_df <- reshape2::melt(
      data.frame(Sample = seq_len(nrow(z_ord)), z_ord),
      id.vars = "Sample", variable.name = "Component", value.name = "Probability"
    )

    class_prob_plot <- ggplot2::ggplot(z_df,
                                       ggplot2::aes(x = Sample, y = Component,
                                                    fill = Probability)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradientn(
        colours = c("white", "#2166AC", "#D6604D"),
        limits  = c(0, 1),
        name    = "Posterior\nProbability"
      ) +
      ggplot2::labs(
        title    = "Classification Posterior Probabilities",
        subtitle = paste0("K = ", optimal_k, ", samples sorted by assigned cluster"),
        x        = "Samples (sorted by cluster)",
        y        = "Component"
      ) +
      pub_theme +
      ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())

    if (save_plots) {
      ggplot2::ggsave(file.path(save_dir, "lpa_class_prob_heatmap.pdf"), class_prob_plot,
                      width = plot_width + 1, height = plot_height - 1, dpi = 300)
      ggplot2::ggsave(file.path(save_dir, "lpa_class_prob_heatmap.png"), class_prob_plot,
                      width = plot_width + 1, height = plot_height - 1, dpi = 150)
      if (verbose) cat("- Classification probability heatmap saved.\n")
    }
  }

  # ====================================================================
  # Plot 5: Silhouette bar chart
  # ====================================================================
  silhouette_plot <- NULL
  if (show_silhouette && optimal_k > 1) {
    sil    <- cluster::silhouette(cluster_labels, stats::dist(data))
    sil_df <- as.data.frame(sil[, 1:3])
    colnames(sil_df) <- c("cluster", "neighbor", "sil_width")
    sil_df$order <- rank(-sil_df$sil_width + sil_df$cluster * 1e6,
                         ties.method = "first")
    avg_sil <- mean(sil_df$sil_width)

    silhouette_plot <- ggplot2::ggplot(sil_df,
                                       ggplot2::aes(x = order, y = sil_width,
                                                    fill = factor(cluster))) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::geom_hline(yintercept = avg_sil, linetype = "dashed",
                          colour = "#E84040", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = cols, name = "Cluster") +
      ggplot2::labs(
        title    = paste0("Silhouette Plot  (K = ", optimal_k, ")"),
        subtitle = paste0("Average silhouette width = ", round(avg_sil, 3)),
        x        = "Samples (ordered by cluster & silhouette width)",
        y        = "Silhouette Width"
      ) +
      ggplot2::ylim(c(min(-0.2, min(sil_df$sil_width) - 0.05), 1)) +
      pub_theme +
      ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())

    if (save_plots) {
      ggplot2::ggsave(file.path(save_dir, "lpa_silhouette_plot.pdf"), silhouette_plot,
                      width = plot_width, height = plot_height, dpi = 300)
      ggplot2::ggsave(file.path(save_dir, "lpa_silhouette_plot.png"), silhouette_plot,
                      width = plot_width, height = plot_height, dpi = 150)
      if (verbose) cat("- Silhouette plot saved.\n")
    }
  } else if (show_silhouette && optimal_k == 1) {
    warning("Silhouette plot skipped: only one cluster found.")
  }

  # ---- Save cluster assignments ----
  clustered_data <- data.frame(data, cluster = factor(cluster_labels))
  if (save_plots) {
    write.csv(clustered_data,
              file.path(save_dir, "lpa_cluster_assignments.csv"),
              row.names = TRUE)
    if (verbose) cat("- Cluster assignments saved.\n")
  }

  # ---- Print plots ----
  if (!is.null(bic_plot))         print(bic_plot)
  if (!is.null(pca_plot_obj))     print(pca_plot_obj)
  if (!is.null(profile_plot))     print(profile_plot)
  if (!is.null(class_prob_plot))  print(class_prob_plot)
  if (!is.null(silhouette_plot))  print(silhouette_plot)

  return(list(
    clustered_data = clustered_data,
    cluster_labels = cluster_labels,
    optimal_model  = optimal_model,
    optimal_k      = optimal_k,
    lpa_object     = lpa_result,
    plots = list(
      bic         = bic_plot,
      bic_summary = if (exists("bic_summary_plot")) bic_summary_plot else NULL,
      pca         = pca_plot_obj,
      profile     = profile_plot,
      class_prob  = class_prob_plot,
      silhouette  = silhouette_plot
    )
  ))
}

# =============================================================================
# SECTION 2: Subtyping-object wrapper
# =============================================================================

#' Sub LPA with Optimal K
#'
#' Wraps \code{lpa_with_optimal_k} for use with a \code{Subtyping} S4 object.
#'
#' @param object          Subtyping object (or numeric data frame).
#' @param use_scaled_data Use the \code{scale.data} slot when \code{TRUE}.
#' @param max_clusters    Maximum number of components.
#' @param model_names     mclust covariance model names.
#' @param save_plots      Whether to save plots.
#' @param save_dir        Output directory.
#' @param plot_width      Plot width (inches).
#' @param plot_height     Plot height (inches).
#' @param base_size       Base font size.
#' @param seed            Random seed.
#' @param verbose         Print progress messages.
#' @param show_silhouette Include silhouette plot.
#' @param show_profile    Include variable profile plot.
#' @param show_class_prob Include classification probability heatmap.
#' @param color_palette   wesanderson palette for cluster colours.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object with scale.data
#'   obj <- Sub_lpa_with_optimal_k(obj, max_clusters = 5, save_plots = FALSE)
#'   print(obj@Optimal.cluster)
#' }
#' @export
Sub_lpa_with_optimal_k <- function(object,
                                   use_scaled_data = TRUE,
                                   max_clusters    = 5,
                                   model_names     = c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
                                                       "EEE", "EEV", "VEV", "EVV",
                                                       "VVE", "VEE", "EVE", "VVV"),
                                   save_plots      = TRUE,
                                   save_dir        = NULL,
                                   plot_width      = 7,
                                   plot_height     = 5,
                                   base_size       = 14,
                                   seed            = 123,
                                   verbose         = TRUE,
                                   show_silhouette = TRUE,
                                   show_profile    = TRUE,
                                   show_class_prob = TRUE,
                                   color_palette   = "Darjeeling1") {

  if (inherits(object, "Subtyping")) {
    data <- if (use_scaled_data) slot(object, "scale.data") else slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be a 'Subtyping' object or a data frame.")
  }

  if (is.null(data) || nrow(data) == 0) stop("No valid data found in the input.")

  numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
  if (ncol(numeric_data) == 0) stop("No numeric columns available for clustering.")

  if (verbose) cat("Starting LPA clustering analysis...\n")

  lpa_result <- lpa_with_optimal_k(
    numeric_data,
    max_clusters    = max_clusters,
    model_names     = model_names,
    save_plots      = save_plots,
    save_dir        = save_dir,
    plot_width      = plot_width,
    plot_height     = plot_height,
    base_size       = base_size,
    seed            = seed,
    verbose         = verbose,
    show_silhouette = show_silhouette,
    show_profile    = show_profile,
    show_class_prob = show_class_prob,
    color_palette   = color_palette
  )

  # Rename 'cluster' -> 'group' in clustered_data
  clustered_data <- lpa_result$clustered_data
  colnames(clustered_data)[colnames(clustered_data) == "cluster"] <- "group"

  if (inherits(object, "Subtyping")) {
    # Write cluster_lpa into info.data
    if (nrow(object@info.data) == 0)
      object@info.data <- data.frame(row.names = rownames(object@clean.data))
    object@info.data <- object@info.data[rownames(object@clean.data), , drop = FALSE]
    idx <- match(rownames(object@info.data),
                 rownames(lpa_result$clustered_data))
    object@info.data$cluster_lpa <- lpa_result$cluster_labels[idx]

    object@cluster.results[["lpa.result"]] <- lpa_result
    object@Optimal.cluster <- lpa_result$optimal_k
    object@clustered.data  <- clustered_data

    if (verbose) {
      cat("\nUpdating 'Subtyping' object...\n")
      cat("- cluster.results[lpa.result] updated\n")
      cat("- Optimal.cluster:", lpa_result$optimal_k, "\n")
      cat("- clustered.data updated\n")
      cat("- info.data$cluster_lpa written\n")
      cat("- Selected model:", lpa_result$optimal_model, "\n")
    }
    return(object)
  } else {
    return(lpa_result)
  }
}
