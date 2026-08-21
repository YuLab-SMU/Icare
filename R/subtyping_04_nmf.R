# =============================================================================
# NMF Subtyping Module
# =============================================================================
# Required: NMF, nnls, here, ggplot2, gridExtra, pheatmap, cluster,
#           RColorBrewer, wesanderson, reshape2

#' Extract cluster assignments from a fitted NMF model
#'
#' Internal helper: each sample assigned to the component with highest
#' coefficient weight (max row of H matrix).
#'
#' @param fit NMFfit object.
#' @return Integer vector of cluster assignments.
#' @keywords internal
#' @noRd
.nmf_get_groups <- function(fit) {
  if (missing(fit)) stop("'fit' is required.")
  if (!inherits(fit, c("NMFfit", "NMFfitX1"))) {
    stop("'fit' must be an NMF object.")
  }
  H <- NMF::coef(fit)  # components × samples
  if (nrow(H) == 0 || ncol(H) == 0) {
    stop("Coefficient matrix is empty.")
  }
  as.integer(apply(H, 2, which.max))
}

# =============================================================================
# SECTION 1: Rank Evaluation
# =============================================================================

#' Generate NMF rank-evaluation plots
#'
#' Produces four panels: Cophenetic Correlation, Reconstruction Residual,
#' Dispersion, and RSS - saved as PDF + PNG.
#'
#' @param estimate  NMF.rank object (from \code{NMF::nmf(..., rank = k1:k2)}).
#' @param save_dir  Output directory.
#' @param width     PDF width (inches).
#' @param height    PDF height (inches).
#' @param base_size Base font size.
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs scale_x_continuous ggsave
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'estimate' is an NMF.rank object from NMF::nmf()
#'   generate_nmf_rank_plots(estimate, save_dir = "./nmf_results")
#' }
generate_nmf_rank_plots <- function(estimate,
                                    save_dir  = NULL,
                                    width     = 10,
                                    height    = 6,
                                    base_size = 14) {

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  pub_theme <- .pub_theme(base_size)

  # measures is a data.frame with columns: rank, cophenetic, rss, dispersion, residuals, ...
  ms   <- summary(estimate)
  if (is.null(ms) || is.null(nrow(ms)) || nrow(ms) == 0) {
    message("DEBUG estimate class: ", class(estimate))
    message("DEBUG ms class: ", class(ms))
    stop("No measures found in NMF estimate.")
  }
  ranks <- ms$rank

  make_panel <- function(y_col, y_label, title, colour) {
    df <- data.frame(rank = ranks, y = ms[[y_col]])
    ggplot2::ggplot(df, ggplot2::aes(x = rank, y = y)) +
      ggplot2::geom_line(colour = colour, linewidth = 1.2) +
      ggplot2::geom_point(size = 3.5, colour = colour, fill = "white",
                          shape = 21, stroke = 1.5) +
      ggplot2::scale_x_continuous(breaks = ranks) +
      ggplot2::labs(title = title, x = "Rank (k)", y = y_label) +
      pub_theme
  }

  p1 <- make_panel("cophenetic", "Cophenetic Correlation",
                   "Cophenetic Correlation", "#3A7DC9")
  p2 <- make_panel("residuals",  "Residual",
                   "Reconstruction Residual", "#C0392B")
  p3 <- make_panel("dispersion", "Dispersion",
                   "Dispersion", "#27AE60")
  p4 <- make_panel("rss",        "RSS",
                   "RSS", "#8E44AD")

  combined <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

  ggplot2::ggsave(file.path(save_dir, "nmf_rank_metrics.pdf"),
                  combined, width = width, height = height)
  ggplot2::ggsave(file.path(save_dir, "nmf_rank_metrics.png"),
                  combined, width = width, height = height, dpi = 150)
  cat("- NMF rank metrics saved to:", save_dir, "\n")
  grid::grid.newpage(); grid::grid.draw(combined)
  return(invisible(combined))
}

#' Sub NMF Estimate (Subtyping wrapper)
#'
#' Runs multi-rank NMF decomposition, stores the NMF.rank estimate,
#' and saves rank-evaluation plots.
#'
#' @param object     Subtyping object (or numeric matrix).
#' @param rank_range Integer vector of ranks to evaluate.
#' @param seed       Random seed.
#' @param nmf_shift  nmf_shift.
#' @param nrun       Number of NMF runs per rank.
#' @param method     NMF algorithm (default \code{"brunet"}).
#' @param save_dir   Output directory.
#' @param base_size  Base font size for plots.
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' require('NMF')
#'set.seed(1)
#'demo_df <- data.frame(
#'  group = rep(c(0, 1), each = 30),
#'  feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'  feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'  feat3 = rnorm(60, 10, 2),
#'  feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#')
#'rownames(demo_df) <- paste0("S", 1:60)
#'
#'stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                            group_col = "group", na.action = "allow")
#'sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#'sub_obj <- Sub_normalize_process(sub_obj, normalize_method = "min_max")
#'
#'sub_obj <- Sub_nmf_estimate(sub_obj, rank_range = 2:3, nrun = 2,
#'                            method = "brunet", seed = 1)
#' sub_obj <- Sub_nmf_best_rank(sub_obj, nrun = 2, method = "brunet")
#' sub_obj <- Sub_nmf_assign_subtypes(sub_obj)
#' table(sub_obj@info.data$cluster_nmf)
#' }
Sub_nmf_estimate <- function(object,
                             rank_range = 2:4,
                             seed       = 8891,
                             nmf_shift  = 0,
                             nrun       = 10,
                             method     = "brunet",
                             save_dir   = file.path(get_output_dir("m3", "cluster_results"), "nmf_results"),
                             base_size  = 14) {
  .require_pkgs("NMF")
  data <- .extract_data(object, scaled = TRUE)

  if (any(data < 0)) {
    warning("Negative values detected. NMF requires non-negative input. ",
            "Currently applying shift = ", nmf_shift, " (add to all values). ",
            "Consider using a larger shift if your data contains meaningful negatives.")
    data <- data + nmf_shift
  }

  cat("Starting NMF rank estimation (ranks:", paste(rank_range, collapse = "\u2013"), ")...\n")
  # NMF requires features x samples; .extract_data returns samples x features -> transpose
  estimate <- NMF::nmf(t(as.matrix(data)),
                       rank   = rank_range,
                       method = method,
                       nrun   = nrun,
                       seed   = seed)
  cat("NMF decomposition complete.\n")

  if (length(rank_range) > 1) {
    generate_nmf_rank_plots(estimate, save_dir = save_dir, base_size = base_size)
  }

  if (inherits(object, "Subtyping")) {
    object@cluster.results[["nmf.result"]] <- list(estimate = estimate)
    if (length(rank_range) == 1) {
      object@cluster.results[["nmf.result"]][["best_estimate"]] <- estimate
      object@Optimal.cluster <- rank_range
      cat("- cluster.results[nmf.result$best_estimate] also updated since only 1 rank was tested\n")
    }
    cat("Updating 'Subtyping' object...\n")
    cat("- cluster.results[nmf.result$estimate] updated\n")
    return(object)
  }
  return(estimate)
}

# =============================================================================
# SECTION 2: Best Rank Analysis
# =============================================================================

#' Consensus heatmap for a fitted NMF object
#'
#' @param fit      Fitted NMF object (single rank).
#' @param save_dir Output directory.
#' @param width    PDF width.
#' @param height   PDF height.
#' @param k        Rank (used for filename/title).
#' @importFrom pheatmap pheatmap
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'fit' is a fitted NMF object
#'   nmf_consensus_heatmap(fit, save_dir = ".", k = 3)
#' }
nmf_consensus_heatmap <- function(fit,
                                  save_dir = NULL,
                                  width    = 8,
                                  height   = 7,
                                  k        = NULL) {

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  cons   <- NMF::consensus(fit)
  groups <- tryCatch(.nmf_get_groups(fit), error = function(e) NULL)

  ann <- if (!is.null(groups)) {
    a <- data.frame(Subtype = paste0("S", groups))
    rownames(a) <- colnames(cons)
    a
  } else NULL

  title    <- if (!is.null(k)) paste0("Consensus Matrix  (K = ", k, ")") else "Consensus Matrix"
  filename <- if (!is.null(k)) paste0("nmf_consensus_K", k, ".pdf") else "nmf_consensus.pdf"

  pdf(file.path(save_dir, filename), width = width, height = height)
  ph <- pheatmap::pheatmap(
    cons,
    annotation_col = ann,
    annotation_row = ann,
    color          = colorRampPalette(c("white", "#154360"))(100),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    main           = title,
    fontsize       = 8,
    border_color   = NA
  )
  invisible(dev.off())
  cat("- Consensus map saved:", filename, "\n")
}

#' Basis map (W matrix heatmap) for a fitted NMF object
#'
#' Shows the feature (gene/marker) loadings per subtype component.
#' Inspired by EcoTyper's metagene heatmap approach.
#'
#' @param fit           Fitted NMF object.
#' @param save_dir      Output directory.
#' @param top_n         Number of top features per component to highlight.
#' @param palette_name  wesanderson palette for component colours.
#' @param width         PDF width.
#' @param height        PDF height.
#' @param base_size     Base font size.
#' @param k             Rank.
#' @importFrom pheatmap pheatmap
#' @importFrom wesanderson wes_palette
#' @examples
#' \dontrun{
#'   # Assuming 'fit' is a fitted NMF object with rank 3
#'   nmf_basis_heatmap(fit, save_dir = "./nmf_results", k = 3, top_n = 20)
#' }
#' @export
nmf_basis_heatmap <- function(fit,
                              save_dir     = NULL,
                              top_n        = 20,
                              palette_name = "Zissou1",
                              width        = 10,
                              height       = 8,
                              base_size    = 12,
                              k            = NULL) {

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  W      <- NMF::basis(fit)          # features x components
  n_comp <- ncol(W)
  colnames(W) <- paste0("S", seq_len(n_comp))

  # Select top_n features per component by highest loading
  top_idx_list <- lapply(seq_len(n_comp), function(j) {
    order(W[, j], decreasing = TRUE)[seq_len(min(top_n, nrow(W)))]
  })
  # Keep assignment: each feature -> dominant component (for row gap lines)
  feature_comp <- unlist(lapply(seq_len(n_comp), function(j) rep(j, length(top_idx_list[[j]]))))
  top_idx      <- unlist(top_idx_list)
  # Remove duplicates, keep first occurrence (dominant assignment wins)
  dup_mask <- duplicated(top_idx)
  top_idx      <- top_idx[!dup_mask]
  feature_comp <- feature_comp[!dup_mask]

  W_top <- W[top_idx, , drop = FALSE]

  # Row-scale: z-score across components for each feature
  W_scaled <- t(scale(t(W_top)))
  W_scaled[is.nan(W_scaled) | is.infinite(W_scaled)] <- 0

  # Row annotation: which component each feature belongs to
  comp_cols <- tryCatch(
    wesanderson::wes_palette(palette_name, n_comp, type = "continuous"),
    error = function(e) scales::hue_pal()(n_comp)
  )
  comp_lvls <- paste0("S", seq_len(n_comp))
  ann_row   <- data.frame(Component = factor(paste0("S", feature_comp), levels = comp_lvls),
                          row.names = rownames(W_top))
  ann_col_df <- data.frame(Component = factor(comp_lvls, levels = comp_lvls),
                           row.names = comp_lvls)
  col_list <- list(Component = setNames(comp_cols, comp_lvls))

  # Row gaps between components
  gaps_row <- cumsum(table(factor(feature_comp, levels = seq_len(n_comp))))[-n_comp]

  title    <- if (!is.null(k)) paste0("NMF Basis Matrix  (K = ", k, ")") else "NMF Basis Matrix"
  filename <- if (!is.null(k)) paste0("nmf_basis_K", k, ".pdf") else "nmf_basis.pdf"

  # Dynamic height
  fig_height <- max(height, length(top_idx) * 0.22 + 2)

  pdf(file.path(save_dir, filename), width = width, height = fig_height)
  ph <- pheatmap::pheatmap(
    W_scaled,
    annotation_row    = ann_row,
    annotation_col    = ann_col_df,
    annotation_colors = col_list,
    color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
    cluster_rows      = FALSE,   # already ordered by component
    cluster_cols      = FALSE,
    gaps_row          = gaps_row,
    show_rownames     = nrow(W_scaled) <= 80,
    fontsize          = base_size - 2,
    fontsize_row      = max(6, base_size - 4),
    main              = title,
    border_color      = NA
  )
  invisible(dev.off())
  cat("- Basis heatmap saved:", filename, "\n")
}

#' Coefficient map (H matrix heatmap) for a fitted NMF object
#'
#' Shows the sample membership per subtype component (EcoTyper-style).
#'
#' @param fit           Fitted NMF object.
#' @param save_dir      Output directory.
#' @param palette_name  wesanderson palette.
#' @param width         PDF width.
#' @param height        PDF height.
#' @param k             Rank.
#' @importFrom pheatmap pheatmap
#' @examples
#' \dontrun{
#'   # Assuming 'fit' is a fitted NMF object with rank 3
#'   nmf_coef_heatmap(fit, save_dir = "./nmf_results", k = 3)
#' }
#' @export
nmf_coef_heatmap <- function(fit,
                             save_dir     = NULL,
                             palette_name = "Zissou1",
                             width        = 10,
                             height       = 5,
                             k            = NULL) {

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  H      <- NMF::coef(fit)            # components x samples
  n_comp <- nrow(H)
  rownames(H) <- paste0("S", seq_len(n_comp))

  # Column-normalize: each sample sums to 1 (membership proportions)
  H_norm <- sweep(H, 2, colSums(H) + 1e-10, "/")

  # Predicted subtype per sample
  groups <- tryCatch(as.integer(.nmf_get_groups(fit)), error = function(e) NULL)

  if (!is.null(groups)) {
    # Sort samples by subtype -> within subtype by dominant component value
    ord      <- order(groups, -H_norm[cbind(groups, seq_along(groups))])
    H_norm   <- H_norm[, ord, drop = FALSE]
    groups   <- groups[ord]
    # gaps between subtypes
    gaps     <- which(diff(groups) != 0)
  } else {
    gaps <- NULL
  }

  # Annotation bar (top of heatmap)
  comp_cols <- tryCatch(
    wesanderson::wes_palette(palette_name, n_comp, type = "continuous"),
    error = function(e) scales::hue_pal()(n_comp)
  )
  subtype_lvls <- paste0("S", seq_len(n_comp))
  ann_col <- if (!is.null(groups)) {
    data.frame(Subtype = factor(paste0("S", groups), levels = subtype_lvls),
               row.names = colnames(H_norm))
  } else NULL
  ann_colors <- list(Subtype = setNames(comp_cols, subtype_lvls))

  title    <- if (!is.null(k)) paste0("NMF Coefficient Matrix  (K = ", k, ")") else "NMF Coefficient Matrix"
  filename <- if (!is.null(k)) paste0("nmf_coef_K", k, ".pdf") else "nmf_coef.pdf"

  # Dynamic height: give more room when many components
  fig_height <- max(height, n_comp * 0.8 + 2)

  pdf(file.path(save_dir, filename), width = width, height = fig_height)
  ph <- pheatmap::pheatmap(
    H_norm,
    annotation_col    = ann_col,
    annotation_colors = ann_colors,
    color             = colorRampPalette(c("#F7FBFF", "#2171B5", "#08306B"))(100),
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,      # already sorted by subtype
    gaps_col          = gaps,       # vertical dividers between subtypes
    show_colnames     = FALSE,
    fontsize          = 11,
    fontsize_row      = 12,
    main              = title,
    border_color      = NA,
    legend_breaks     = c(0, 0.25, 0.5, 0.75, 1),
    legend_labels     = c("0", "0.25", "0.5", "0.75", "1")
  )
  invisible(dev.off())
  cat("- Coefficient heatmap saved:", filename, "\n")
}

#' Silhouette plot for NMF clustering
#'
#' @param fit       Fitted NMF object.
#' @param save_dir  Output directory.
#' @param palette_name wesanderson palette.
#' @param plot_width  Plot width.
#' @param plot_height Plot height.
#' @param base_size   Base font size.
#' @param k           Rank.
#' @importFrom cluster silhouette
#' @importFrom ggplot2 ggplot aes geom_bar geom_hline geom_text labs ylim ggsave
#' @examples
#' \dontrun{
#'   # Assuming 'fit' is a fitted NMF object with rank 3
#'   sil_info <- nmf_silhouette_plot(fit, save_dir = "./nmf_results", k = 3)
#'   print(sil_info$avg_width)
#' }
#' @export
nmf_silhouette_plot <- function(fit,
                                save_dir     = NULL,
                                palette_name = "Darjeeling1",
                                plot_width   = 7,
                                plot_height  = 5,
                                base_size    = 14,
                                k            = NULL) {

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  pub_theme <- .pub_theme(base_size)
  cons      <- NMF::consensus(fit)
  groups    <- as.numeric(.nmf_get_groups(fit))
  n_groups  <- length(unique(groups))

  sil <- cluster::silhouette(groups, stats::as.dist(1 - cons))
  sil_df <- as.data.frame(sil[, 1:3])
  colnames(sil_df) <- c("cluster", "neighbor", "sil_width")
  sil_df$order <- rank(-sil_df$sil_width + sil_df$cluster * 1e6, ties.method = "first")
  avg_sil <- mean(sil_df$sil_width)

  cols <- tryCatch(
    wesanderson::wes_palette(palette_name, n_groups, type = "continuous"),
    error = function(e) scales::hue_pal()(n_groups)
  )

  title    <- if (!is.null(k)) paste0("Silhouette Plot  (K = ", k, ")") else "Silhouette Plot"
  filename <- if (!is.null(k)) paste0("nmf_silhouette_K", k, ".pdf") else "nmf_silhouette.pdf"

  p <- ggplot2::ggplot(sil_df, ggplot2::aes(x = order, y = sil_width,
                                            fill = factor(cluster))) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::geom_hline(yintercept = avg_sil, linetype = "dashed",
                        colour = "#E84040", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = cols, name = "Subtype") +
    ggplot2::labs(
      title    = title,
      subtitle = paste0("Average silhouette width = ", round(avg_sil, 3)),
      x        = "Samples (ordered by subtype & silhouette width)",
      y        = "Silhouette Width"
    ) +
    ggplot2::ylim(c(min(-0.2, min(sil_df$sil_width) - 0.05), 1)) +
    pub_theme +
    ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  ggplot2::ggsave(file.path(save_dir, filename),        p, width = plot_width, height = plot_height, dpi = 300)
  ggplot2::ggsave(file.path(save_dir, sub(".pdf",".png",filename)), p, width = plot_width, height = plot_height, dpi = 150)
  cat("- Silhouette plot saved:", filename, "\n")

  return(list(plot = p, avg_width = avg_sil, silhouette = sil))
}

#' NMF Best Rank Analysis
#'
#' Identifies optimal rank by cophenetic correlation, re-runs NMF at that
#' rank, and produces consensus heatmap, basis heatmap, coefficient heatmap,
#' and silhouette plot.
#'
#' @param object        Subtyping object (or NMF.rank estimate list).
#' @param nrun          Number of NMF runs.
#' @param seed          Random seed.
#' @param nmf_shift     nmf_shift.
#' @param method        NMF algorithm.
#' @param palette_name  wesanderson palette.
#' @param save_dir      Output directory.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object with nmf.result already computed
#'   obj <- Sub_nmf_best_rank(obj, nrun = 20, palette_name = "Zissou1")
#' }
#' @export
Sub_nmf_best_rank <- function(object,
                              nrun         = 10,
                              seed         = 8891,
                              nmf_shift   =  0,
                              method       = "brunet",
                              palette_name = "Zissou1",
                              save_dir     = file.path(get_output_dir("Subtyping", "cluster_results"), "nmf_results")) {

  if (inherits(object, "Subtyping")) {
    data     <- .extract_data(object, scaled = TRUE)
    estimate <- slot(object, "cluster.results")[["nmf.result"]][["estimate"]]
  } else if (is.list(object)) {
    data     <- object$data
    estimate <- object$estimate
  } else {
    stop("Input must be a 'Subtyping' object or a named list.")
  }

  if (is.null(estimate)) stop("Run 'Sub_nmf_estimate' first.")
  if (is.null(data) || nrow(data) == 0) stop("No valid data found.")

  if (any(data < 0)) {
    warning("Negative values detected. NMF requires non-negative input. ",
            "Currently applying shift = ", nmf_shift, " (add to all values). ",
            "Consider using a larger shift if your data contains meaningful negatives.")
    data <- data + nmf_shift
  }

  ms         <- summary(estimate)
  best_idx   <- which.max(ms$cophenetic)
  best_k     <- ms$rank[best_idx]
  cat("Best rank (max cophenetic):", best_k, "\n")

  cat("Re-running NMF at K =", best_k, "for final model...\n")
  # NMF requires features x samples -> transpose
  fit <- NMF::nmf(t(as.matrix(data)), rank = best_k,
                  method = method, nrun = nrun, seed = seed)

  nmf_consensus_heatmap(fit, save_dir = save_dir, k = best_k)
  nmf_basis_heatmap(fit,     save_dir = save_dir, palette_name = palette_name, k = best_k)
  nmf_coef_heatmap(fit,      save_dir = save_dir, palette_name = palette_name, k = best_k)
  sil_info <- nmf_silhouette_plot(fit, save_dir = save_dir,
                                  palette_name = palette_name, k = best_k)

  if (inherits(object, "Subtyping")) {
    object@cluster.results[["nmf.result"]][["best_estimate"]] <- fit
    object@cluster.results[["nmf.result"]][["best_rank"]]     <- best_k
    object@cluster.results[["nmf.result"]][["silhouette"]]    <- sil_info
    object@Optimal.cluster <- best_k
    cat("Updating 'Subtyping' object...\n")
    cat("- best_estimate, best_rank, silhouette updated\n")
    cat("- Optimal.cluster:", best_k, "\n")
    return(object)
  }
  return(list(best_estimate = fit, best_rank = best_k, silhouette = sil_info))
}

# =============================================================================
# SECTION 3: Group Assignment
# =============================================================================

#' Assign samples to NMF subtypes
#'
#' Uses the H (coefficient) matrix: each sample is assigned to the
#' component with highest coefficient weight (same logic as EcoTyper).
#'
#' @param object Subtyping object.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object with best NMF estimate
#'   obj <- Sub_nmf_assign_subtypes(obj)
#'   table(obj@info.data$cluster_nmf)
#' }
#' @export
Sub_nmf_assign_subtypes <- function(object) {

  if (!inherits(object, "Subtyping"))
    stop("Input must be a 'Subtyping' object.")

  fit <- slot(object, "cluster.results")[["nmf.result"]][["best_estimate"]]
  if (is.null(fit)) stop("Run 'Sub_nmf_best_rank' first.")

  H              <- NMF::coef(fit)               # components x samples
  group_vec      <- apply(H, 2, which.max)
  group_fac      <- factor(paste0("S", group_vec))

  cat("Subtype assignment (H-matrix max coefficient):\n")
  print(table(group_fac))

  if (nrow(object@info.data) == 0)
    object@info.data <- data.frame(row.names = rownames(object@clean.data))
  object@info.data <- object@info.data[rownames(object@clean.data), , drop = FALSE]
  object@info.data$cluster_nmf <- group_fac[match(rownames(object@info.data),
                                                   names(group_vec))]

  cdata            <- object@clean.data
  cdata$group      <- group_fac[match(rownames(cdata), names(group_vec))]
  object@clustered.data <- cdata

  cat("Updating 'Subtyping' object...\n")
  cat("- clustered.data updated\n")
  cat("- info.data$cluster_nmf written\n")
  return(object)
}

# =============================================================================
# SECTION 4: Model Training & Prediction (EcoTyper-inspired)
# =============================================================================

#' Train final NMF model and save all outputs
#'
#' Trains the best-rank NMF with more runs and saves the model, consensus
#' matrix, subtype assignments, and silhouette summary.
#'
#' @param object     Subtyping object.
#' @param best_k     Optimal rank.
#' @param nrun       Number of NMF runs for final training.
#' @param method     NMF algorithm.
#' @param model_name File stem for the saved model (.rds).
#' @param save_dir   Output directory.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object with best_k already determined
#'   obj <- Sub_nmf_train_model(obj, best_k = 3, nrun = 50, save_dir = "./nmf_model")
#' }
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' require('NMF')
#' demo_df <- data.frame(
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2),
#'   feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#'  )
#' rownames(demo_df) <- paste0("S", 1:60)  
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' split <- SplitSubtypingObject(sub_obj, p = 0.7, stratify_by = "group")
#' sub_train <- Sub_normalize_process(split$train, normalize_method = "min_max")
#' norm_params <- Sub_extract_norm_params(sub_train, verbose = FALSE)
#' sub_test <- Sub_apply_norm_params(split$test, norm_params = norm_params, verbose = FALSE)
#'
#' sub_train <- Sub_nmf_estimate(sub_train, rank_range = 2:3, nrun = 2, seed = 1)
#' sub_train <- Sub_nmf_best_rank(sub_train, nrun = 2)
#' sub_train <- Sub_nmf_assign_subtypes(sub_train)
#' sub_train <- Sub_nmf_train_model(
#'   sub_train, best_k = sub_train@Optimal.cluster, nrun = 2,
#'   model_name = "demo_nmf_model",save_dir   ="./"
#' )
#'
#' # Predict subtypes on the held-out validation split.
#' sub_test <- Sub_predict_subtypes(sub_test, train_object = sub_train,
#'                                  method = "nmf", verbose = FALSE)
#' table(sub_test@info.data$cluster_nmf)
#' }
Sub_nmf_train_model <- function(object,
                                best_k,
                                nrun       = 50,
                                method     = "brunet",
                                model_name = "nmf_model",
                                save_dir   = NULL) {

  if (!inherits(object, "Subtyping")) stop("Input must be a 'Subtyping' object.")

  train_data <- t(as.matrix(.extract_data(object, scaled = TRUE)))   # features x samples
  if (any(train_data < 0)) train_data[train_data < 0] <- 0

  cat("Training final NMF model (K =", best_k, ", nrun =", nrun, ")...\n")
  set.seed(8891)
  fit <- NMF::nmf(train_data, rank = best_k, method = method, nrun = nrun)

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # ---- Save model ----
  model_path <- file.path(save_dir, paste0(model_name, ".rds"))
  saveRDS(list(
    fit           = fit,
    feature_names = rownames(NMF::basis(fit)),
    H_matrix      = NMF::coef(fit),
    best_k        = best_k
  ), model_path)
  cat("- Model saved to:", model_path, "\n")

  # ---- Save consensus matrix ----
  cons      <- NMF::consensus(fit)
  cons_path <- file.path(save_dir, "consensus_matrix.csv")
  write.csv(cons, cons_path)
  cat("- Consensus matrix saved to:", cons_path, "\n")

  # ---- Save subtype assignments ----
  H         <- NMF::coef(fit)
  groups    <- apply(H, 2, which.max)
  group_df  <- data.frame(
    Sample  = colnames(H),
    Subtype = paste0("S", groups),
    row.names = NULL
  )
  assign_path <- file.path(save_dir, "training_subtypes.csv")
  write.csv(group_df, assign_path, row.names = FALSE)
  cat("- Subtype assignments saved to:", assign_path, "\n")

  # ---- Silhouette ----
  sil_info   <- nmf_silhouette_plot(fit, save_dir = save_dir,
                                    k = best_k)
  sil_path   <- file.path(save_dir, "silhouette_summary.txt")
  sink(sil_path)
  cat("Average silhouette width:", round(sil_info$avg_width, 4), "\n\n")
  print(sil_info$silhouette)
  sink()
  cat("- Silhouette summary saved to:", sil_path, "\n")

  # ---- Update object ----
  object@Optimal.cluster <- best_k
  object@cluster.results[["nmf.result"]][["best_estimate"]]            <- fit
  object@cluster.results[["nmf.result"]][["silhouette"]]               <- sil_info
  object@cluster.results[["nmf.result"]][["saved_model_path"]]         <- model_path
  object@cluster.results[["nmf.result"]][["consensus_matrix_path"]]    <- cons_path
  object@cluster.results[["nmf.result"]][["subtype_assignments_path"]] <- assign_path

  # Assign subtypes to clustered.data
  group_fac            <- factor(paste0("S", groups))
  cdata                <- object@clean.data
  cdata$group          <- group_fac[match(rownames(cdata), names(groups))]
  object@clustered.data <- cdata

  cat("All outputs saved. Returning updated Subtyping object.\n")
  return(object)
}

#' Predict subtypes for new samples using a saved NMF model
#'
#' Projects new data onto the saved H matrix via NNLS (same approach
#' as EcoTyper's recovery step).
#'
#' @param model_path Path to .rds model saved by \code{Sub_nmf_train_model}.
#' @param new_data   New data matrix: samples x features (rows = samples).
#' @return Data frame with predicted subtypes and membership probabilities.
#' @examples
#' \dontrun{
#'   # Assuming 'new_data' is a data frame with same features as training data
#'   pred <- Sub_nmf_predict("./nmf_model/nmf_model.rds", new_data)
#'   head(pred)
#' }
#' @export
Sub_nmf_predict <- function(model_path, new_data) {

  if (!file.exists(model_path)) stop("Model file not found: ", model_path)

  mdl      <- readRDS(model_path)
  new_data <- as.matrix(new_data)

  feature_names  <- mdl$feature_names
  common_features <- intersect(feature_names, colnames(new_data))

  if (length(common_features) < length(feature_names))
    warning(length(feature_names) - length(common_features),
            " model features missing in new_data - filled with 0.")

  input_mat <- matrix(0,
                      nrow = nrow(new_data),
                      ncol = length(feature_names),
                      dimnames = list(rownames(new_data), feature_names))
  input_mat[, common_features] <- new_data[, common_features]

  H     <- mdl$H_matrix                                        # components x features
  W_new <- t(apply(input_mat, 1, function(x) nnls::nnls(t(H), x)$x))
  prob  <- t(apply(W_new, 1, function(x) {
    s <- sum(x)
    if (s == 0) x else x / s
  }))
  colnames(prob) <- paste0("Prob_S", seq_len(ncol(prob)))

  data.frame(
    SampleID          = rownames(new_data),
    Predicted_Subtype = paste0("S", apply(W_new, 1, which.max)),
    Confidence        = apply(prob, 1, max),
    prob,
    stringsAsFactors  = FALSE
  )
}

#' Visualise NMF prediction results
#'
#' Produces a boxplot of prediction confidence and a membership
#' probability heatmap.
#'
#' @param pred_df   Data frame from \code{Sub_nmf_predict}.
#' @param save_dir  Output directory.
#' @param palette_name wesanderson palette.
#' @param plot_width  Plot width.
#' @param plot_height Plot height.
#' @param base_size   Base font size.
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange
#' @examples
#' \dontrun{
#'   # Assuming 'pred' is from Sub_nmf_predict()
#'   Sub_nmf_plot_prediction(pred, save_dir = "./nmf_results")
#' }
#' @export
Sub_nmf_plot_prediction <- function(pred_df,
                                    save_dir     = NULL,
                                    palette_name = "Darjeeling1",
                                    plot_width   = 12,
                                    plot_height  = 5,
                                    base_size    = 14) {

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  pub_theme <- .pub_theme(base_size)

  n_sub  <- length(unique(pred_df$Predicted_Subtype))
  cols   <- tryCatch(
    wesanderson::wes_palette(palette_name, n_sub, type = "continuous"),
    error = function(e) scales::hue_pal()(n_sub)
  )

  p1 <- ggplot2::ggplot(pred_df,
                        ggplot2::aes(x = Predicted_Subtype, y = Confidence,
                                     fill = Predicted_Subtype)) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.size = 1.5) +
    ggplot2::geom_jitter(width = 0.18, alpha = 0.55, size = 1.5) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::labs(title = "Prediction Confidence by Subtype",
                  x = "Predicted Subtype", y = "Confidence Score") +
    pub_theme

  prob_cols <- grep("^Prob_S", colnames(pred_df), value = TRUE)
  prob_mat  <- as.matrix(pred_df[, prob_cols, drop = FALSE])
  rownames(prob_mat) <- pred_df$SampleID

  ann_row  <- data.frame(Subtype = pred_df$Predicted_Subtype,
                         row.names = pred_df$SampleID)

  p2 <- pheatmap::pheatmap(
    t(prob_mat),
    annotation_col  = ann_row,
    main            = "Membership Probability Heatmap",
    color           = colorRampPalette(c("white", "#C0392B"))(100),
    show_colnames   = nrow(pred_df) <= 80,
    display_numbers = nrow(pred_df) <= 30,
    border_color    = NA,
    silent          = TRUE
  )

  combined <- gridExtra::grid.arrange(p1, p2$gtable, ncol = 2)

  ggplot2::ggsave(file.path(save_dir, "nmf_prediction_results.pdf"),
                  combined, width = plot_width, height = plot_height)
  cat("- Prediction plot saved to:", save_dir, "\n")

  return(invisible(combined))
}

# =============================================================================
# SECTION 5: Convenience (consensusmap / NMF native plots)
# =============================================================================

#' Plot NMF rank evaluation using NMF native graphics
#'
#' @param object   Subtyping object.
#' @param save_dir Output directory.
#' @param width    PDF width.
#' @param height   PDF height.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a Subtyping object with nmf.result
#'   Sub_nmf_plot_eval(obj, save_dir = "./nmf_results")
#' }
#' @export
Sub_nmf_plot_eval <- function(object,
                              save_dir = NULL,
                              width    = 15,
                              height   = 12) {

  estimate <- object@cluster.results[["nmf.result"]][["estimate"]]
  if (is.null(estimate)) stop("Run 'Sub_nmf_estimate' first.")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  pdf_path <- file.path(save_dir, "NMF_rank_metrics_native.pdf")
  grDevices::pdf(pdf_path, width = width, height = height)
  on.exit({
    tryCatch(grDevices::dev.off(), error = function(e) NULL)
  }, add = TRUE)
  tryCatch(print(plot(estimate)),
           error = function(e) message("NMF rank plot error: ", e$message))
  tryCatch(grDevices::dev.off(), error = function(e) NULL)
  on.exit(NULL)
  cat("- Native rank evaluation plot saved to:", pdf_path, "\n")
  # Also print to screen
  tryCatch(print(plot(estimate)),
           error = function(e) message("NMF rank plot screen error: ", e$message))
}

#' Plot NMF Consensus Maps for All Ranks
#'
#' Generates consensus heatmaps for each rank evaluated in the NMF estimate
#' using \code{NMF::consensusmap}. All maps are saved as a multi-page PDF file
#' and also printed to the screen.
#'
#' @param object   A \code{Subtyping} S4 object containing NMF estimation results
#'   in \code{object@cluster.results$nmf.result$estimate}.
#' @param save_dir Output directory where the PDF will be saved. If \code{NULL},
#'   the current working directory is used. The function creates the directory
#'   if it does not exist.
#' @param width    Width of the PDF page in inches. Default: \code{12}.
#' @param height   Height of the PDF page in inches. Default: \code{10}.
#'
#' @return Invisibly returns \code{NULL}. The function is called for its
#'   side effects: generating and saving consensus maps.
#'
#' @importFrom grDevices pdf dev.off
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'sub_obj' has nmf.result$estimate
#' Sub_nmf_plot_consensus_all(sub_obj, save_dir = "./nmf_output", width = 10, height = 8)
#' }
Sub_nmf_plot_consensus_all <- function(object,
                                       save_dir = NULL,
                                       width    = 12,
                                       height   = 10) {

  estimate <- object@cluster.results[["nmf.result"]][["estimate"]]
  if (is.null(estimate)) stop("Run 'Sub_nmf_estimate' first.")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # NMF::consensusmap manages its own graphics device internally (grid-based).
  # Wrapping it in pdf() causes device conflicts and produces a corrupt file.
  # Solution: redirect each rank's map to its own pdf via grDevices::recordPlot
  # OR use the NMF-recommended approach of letting it draw to screen and
  # capturing with pdf() opened BEFORE the call but with on.exit safety.

  pdf_path <- file.path(save_dir, "NMF_consensus_maps_all.pdf")

  # Open pdf device; ensure it is always closed even on error
  grDevices::pdf(pdf_path, width = width, height = height)
  on.exit({
    tryCatch(grDevices::dev.off(), error = function(e) NULL)
  }, add = TRUE)

  tryCatch(
    NMF::consensusmap(estimate),
    error   = function(e) message("Consensus map error: ", e$message),
    warning = function(w) {
      message("Consensus map warning: ", w$message)
      invokeRestart("muffleWarning")
    }
  )

  # Close PDF explicitly before on.exit (on.exit is a safety net only)
  tryCatch(grDevices::dev.off(), error = function(e) NULL)
  on.exit(NULL)   # cancel safety net now that we closed manually

  cat("- All consensus maps saved to:", pdf_path, "\n")

  # Print to screen (separate device - no outer pdf open)
  tryCatch(NMF::consensusmap(estimate),
           error = function(e) message("Consensus map screen print error: ", e$message))
}

# =============================================================================
# SECTION 6: Internal helpers (shared across modules)
# =============================================================================

#' @keywords internal
.extract_data <- function(object, scaled = TRUE) {
  if (inherits(object, "Subtyping")) {
    d <- if (scaled && nrow(slot(object, "scale.data")) > 0)
      slot(object, "scale.data")
    else
      slot(object, "clean.data")
    # IMPORTANT: Sub_normalize_process()/normalize_data() only overwrites the
    # columns it identifies as numeric IN PLACE -- it never drops the
    # remaining columns from the data frame it returns (this is intentional
    # for Module 1's descriptive/reporting use case, where you typically
    # want to keep ID/categorical columns around). If `group_col` wasn't
    # passed to Sub_normalize_process() and doesn't match the real outcome
    # column name, or if clean.data has any other non-numeric column (ID,
    # site, free text, etc.), that column silently survives into scale.data.
    #
    # as.matrix() on a data.frame with even ONE non-numeric column coerces
    # the ENTIRE matrix to character -- not just that column -- and
    # NMF::nmf() (or k-means/LPA) then fails deep inside with a cryptic,
    # unrelated-looking error such as "'x' must be numeric". Catch this
    # here, at the point of conversion, with a clear and specific message
    # instead.
    non_numeric <- names(d)[!vapply(d, is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
      stop(
        "Non-numeric column(s) found in ", if (scaled) "scale.data" else "clean.data",
        ": ", paste(non_numeric, collapse = ", "), ". Clustering methods ",
        "(NMF/K-means/LPA) require a purely numeric matrix. This usually ",
        "means the real outcome/ID column name doesn't match the ",
        "`group_col` argument you passed (or didn't pass) to ",
        "Sub_normalize_process() -- it only excludes a column from ",
        "normalization by exact name match, and never removes non-numeric ",
        "columns from its output. Either pass the correct `group_col`, or ",
        "drop these columns manually before clustering, e.g.:\n",
        "  object@clean.data <- object@clean.data[, setdiff(colnames(object@clean.data), c(",
        paste(sprintf('\"%s\"', non_numeric), collapse = ", "), ")), drop = FALSE]",
        call. = FALSE
      )
    }
    as.matrix(d)
  } else {
    as.matrix(object)
  }
}
