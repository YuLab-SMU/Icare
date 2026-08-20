## ============================================================
##  viz_feature_selection.R
##  Publication-quality Feature Selection Visualization Functions
##  
##  Organized by method:
##     1. RFE    -- Profile plot, variable importance
##     2. GA     -- Evolution trace, fitness history
##     3. SA     -- Temperature schedule, acceptance trace
##     4. SBF    -- Score distribution, selection frequency
##     5. Builtin -- Variable importance comparison
##     6. Integrated -- UpSet, Venn, consensus heatmap, stability
##
##  Every function follows the viz_functions.R framework:
##    - accepts S4 object OR raw result
##    - returns ggplot2 / grob object
##    - saves when save_plot = TRUE
##    - uses .pub_theme + wesanderson palettes
## ============================================================

## -- Load dependencies (if not already loaded) --------------------------------
if (!exists(".pub_theme")) {
  .pub_theme <- function(base_size = 13) {
    if (requireNamespace("ggprism", quietly = TRUE)) {
      ggprism::theme_prism(base_size = base_size)
    } else {
      ggplot2::theme_bw(base_size = base_size)
    } +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                              size = base_size + 1),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
        axis.title    = ggplot2::element_text(face = "bold"),
        legend.title  = ggplot2::element_text(face = "bold"),
        strip.text    = ggplot2::element_text(face = "bold")
      )
  }
}

if (!exists(".get_palette")) {
  .get_palette <- function(palette_name, n) {
    tryCatch({
      if (requireNamespace("wesanderson", quietly = TRUE)) {
        as.character(wesanderson::wes_palette(
          n = n, 
          name = palette_name,
          type = if (n > 5) "continuous" else "discrete"
        ))
      } else {
        RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)]
      }
    }, error = function(e) {
      RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)]
    })
  }
}

if (!exists(".save_plot")) {
  .save_plot <- function(p, dir, filename, width, height, format = "pdf") {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    path <- file.path(dir, paste0(tools::file_path_sans_ext(filename), ".", format))
    ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
    cat("Plot saved to:", path, "\n")
    invisible(path)
  }
}

if (!exists(".get_viz_output_dir")) {
  .get_viz_output_dir <- function(type = "Model") {
    if (exists("get_output_root", mode = "function")) {
      root <- get_output_root()
      if (!is.null(root)) {
        viz_dir <- file.path(root, "Figures", type, "FeatureSelection")
      } else {
        viz_dir <- here::here("Figures", type, "FeatureSelection")
      }
    } else {
      viz_dir <- here::here("Figures", type, "FeatureSelection")
    }
    if (!dir.exists(viz_dir)) {
      dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(viz_dir)
  }
}


## =============================================================================
##   1  RFE -- Recursive Feature Elimination
## =============================================================================

#' RFE Profile Plot (Enhanced)
#'
#' Plots cross-validated performance across different subset sizes with
#' optimal point highlighted.
#'
#' @param rfe_result An \code{rfe} object from caret, or the result list
#'   from \code{FeatureSelectRFE}.
#' @param metric Metric to plot: "Accuracy", "Kappa", "ROC", etc.
#' @param palette_name Wesanderson palette name. Default "Darjeeling1".
#' @param show_optimal Highlight the optimal subset size. Default TRUE.
#' @param base_size Base font size. Default 13.
#' @param save_plot Save plot to file? Default FALSE.
#' @param save_dir Output directory. Auto-detected if NULL.
#' @param width,height Plot dimensions (inches). Defaults 8 x 6.
#' @param format File format: "pdf", "png", "svg". Default "pdf"
#' @param y_limits y_limits.
#' @return A ggplot object.
#' @export
PlotRFE <- function(rfe_result,
                    metric       = "Accuracy",
                    palette_name = "Darjeeling1",
                    show_optimal = TRUE,
                    base_size    = 13,
                    save_plot    = FALSE,
                    save_dir     = NULL,
                    width        = 8,
                    height       = 6,
                    y_limits = c(0, 1),
                    format       = "pdf") {
  
  cat("Generating RFE profile plot...\n")
  if (is.null(save_dir) && save_plot) {
    if (exists(".get_viz_output_dir")) {
      save_dir <- .get_viz_output_dir("Model")
    } else {
      save_dir <- file.path(".", "Figures", "Model", "FeatureSelection")
    }
  }
  
  if (is.list(rfe_result) && "result" %in% names(rfe_result)) {
    rfe_obj <- rfe_result$result
  } else if (inherits(rfe_result, "rfe")) {
    rfe_obj <- rfe_result
  } else {
    stop("Input must be an 'rfe' object or result from FeatureSelectRFE().")
  }
  
  perf_df <- rfe_obj$results
  if (!metric %in% colnames(perf_df)) {
    stop("Metric '", metric, "' not found in results. Available: ",
         paste(colnames(perf_df), collapse = ", "))
  }
  
  opt_size <- rfe_obj$optsize
  
  metric_sd <- paste0(metric, "SD")
  has_sd <- metric_sd %in% colnames(perf_df)
  
  perf_df$Variables   <- perf_df$Variables
  perf_df$Performance <- perf_df[[metric]]
  
  cols <- if (exists(".get_palette")) {
    .get_palette(palette_name, 3)
  } else {
    c("#1b9e77", "#d95f02", "#7570b3")
  }
  
  p <- ggplot2::ggplot(perf_df, ggplot2::aes(x = Variables, y = Performance)) +
    ggplot2::geom_line(colour = cols[1], linewidth = 1) +
    ggplot2::geom_point(colour = cols[2], size = 2.5, alpha = 0.7)
  
  # ---- FIX: use explicit data column ----
  if (has_sd) {
    perf_df$SD <- perf_df[[metric_sd]]
    p <- p + ggplot2::geom_ribbon(
      data = perf_df,
      mapping = ggplot2::aes(x = Variables,
                             ymin = Performance - SD,
                             ymax = Performance + SD),
      alpha = 0.2,
      fill = cols[1],
      inherit.aes = FALSE
    )
  }
  
  if (show_optimal && nrow(perf_df) > 0) {
    opt_row <- perf_df[perf_df$Variables == opt_size, ]
    if (nrow(opt_row) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = opt_row,
          mapping = ggplot2::aes(x = Variables, y = Performance),
          colour = cols[3], size = 5, shape = 21, fill = cols[3], stroke = 1.5,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_vline(
          xintercept = opt_size,
          linetype = "dashed",
          colour = "grey40",
          linewidth = 0.6
        ) +
        ggplot2::annotate(
          "text",
          x = opt_size,
          y = max(perf_df$Performance) * 0.95,
          label = sprintf("Optimal: %d variables", opt_size),
          colour = cols[3],
          fontface = "bold",
          hjust = ifelse(opt_size > median(perf_df$Variables), 1.1, -0.1)
        )
    }
  }
  
  p <- p +
    ggplot2::labs(
      title = "RFE Profile - Recursive Feature Elimination",
      subtitle = sprintf("Selected %d features out of %d tested",
                         opt_size, max(perf_df$Variables)),
      x = "Number of Variables",
      y = metric
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    .pub_theme(base_size)
  
  if (!is.null(y_limits)) {
    if (is.numeric(y_limits) && length(y_limits) == 2) {
      p <- p + ggplot2::scale_y_continuous(limits = y_limits)
    } else {
      warning("y_limits must be a numeric vector of length 2; using adaptive limits.")
    }
  }
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, paste0("RFE_profile.", format))
    ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
    cat("Plot saved to:", path, "\n")
  }
  
  return(p)
}


#' RFE Variable Importance Plot
#'
#' Bar plot of variable importance from the final RFE model.
#'
#' @param rfe_result An \code{rfe} object or result list.
#' @param top_n Number of top variables to display. Default 20.
#' @param palette_name Palette. Default "Moonrise2".
#' @param base_size Font size. Default 13.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 7 x 8.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotRFEImportance <- function(rfe_result,
                              top_n        = 20,
                              palette_name = "Moonrise2",
                              base_size    = 13,
                              save_plot    = FALSE,
                              save_dir     = NULL,
                              width        = 7,
                              height       = 8,
                              format       = "pdf") {
  
  cat("Generating RFE importance plot...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (is.list(rfe_result) && "result" %in% names(rfe_result)) {
    rfe_obj <- rfe_result$result
  } else {
    rfe_obj <- rfe_result
  }
  
  # Extract variable importance
  if (is.null(rfe_obj$variables)) {
    stop("No variable importance information in RFE result.")
  }
  
  # Get importance scores for optimal variables
  opt_vars <- rfe_obj$optVariables
  var_imp <- rfe_obj$variables
  
  # Aggregate importance across resamples
  imp_df <- var_imp %>%
    dplyr::filter(var %in% opt_vars) %>%
    dplyr::group_by(var) %>%
    dplyr::summarise(
      Importance = mean(Overall, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(Importance)) %>%
    head(top_n)
  
  imp_df$var <- factor(imp_df$var, levels = rev(imp_df$var))
  
  cols <- .get_palette(palette_name, 2)
  
  p <- ggplot2::ggplot(imp_df, ggplot2::aes(x = Importance, y = var)) +
    ggplot2::geom_col(fill = cols[1], colour = cols[2], linewidth = 0.3) +
    ggplot2::labs(
      title = "RFE Variable Importance",
      subtitle = sprintf("Top %d selected features", nrow(imp_df)),
      x = "Mean Importance",
      y = NULL
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "bold", colour = "black")
    )
  
  if (save_plot) {
    .save_plot(p, save_dir, "RFE_importance", width, height, format)
  }
  
  return(p)
}


## =============================================================================
##   2  GA -- Genetic Algorithm
## =============================================================================

#' GA Evolution Trace Plot
#'
#' Plots fitness evolution across generations showing best, mean, and
#' worst fitness.
#'
#' @param ga_result A \code{gafs} object or result list from
#'   \code{FeatureSelectGA}.
#' @param metric Fitness metric. Default "Fitness".
#' @param palette_name Palette. Default "Royal2".
#' @param base_size Font size. Default 13.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 8 x 6.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotGA <- function(ga_result,
                   metric       = "Accuracy",
                   palette_name = "Royal2",
                   base_size    = 13,
                   save_plot    = FALSE,
                   save_dir     = NULL,
                   width        = 8,
                   height       = 6,
                   format       = "pdf") {
  
  cat("Generating GA evolution trace...\n")
  if (is.null(save_dir) && save_plot) {
    if (exists(".get_viz_output_dir")) {
      save_dir <- .get_viz_output_dir("Model")
    } else {
      save_dir <- file.path(".", "Figures", "Model", "FeatureSelection")
    }
  }
  
  if (is.list(ga_result) && "result" %in% names(ga_result)) {
    ga_obj <- ga_result$result
  } else if (inherits(ga_result, "gafs")) {
    ga_obj <- ga_result
  } else {
    stop("Input must be a 'gafs' object or result from FeatureSelectGA().")
  }
  
  ext_perf <- ga_obj$external
  if (is.null(ext_perf)) stop("No external performance data in GA result.")
  
  metric_col <- grep(metric, colnames(ext_perf), value = TRUE)[1]
  if (is.na(metric_col)) {
    stop("Metric '", metric, "' not found. Available: ",
         paste(colnames(ext_perf), collapse = ", "))
  }
  
  ext_perf$Generation   <- seq_len(nrow(ext_perf))
  ext_perf$Performance  <- ext_perf[[metric_col]]
  
  # ---- FIX: align feature counts with external rows ----
  n_vars_all <- sapply(ga_obj$ga$final, length)
  n_iters    <- length(n_vars_all)
  n_rows     <- nrow(ext_perf)
  
  if (n_iters == n_rows) {
    ext_perf$NumFeatures <- n_vars_all
  } else if (n_iters > n_rows) {
    # pick evenly spaced indices
    idx <- round(seq(1, n_iters, length.out = n_rows))
    ext_perf$NumFeatures <- n_vars_all[idx]
  } else {
    # pad with last value
    ext_perf$NumFeatures <- c(n_vars_all, rep(n_vars_all[n_iters], n_rows - n_iters))
  }
  
  opt_iter <- which.max(ext_perf$Performance)
  opt_perf <- ext_perf$Performance[opt_iter]
  opt_feat <- ext_perf$NumFeatures[opt_iter]
  
  cols <- if (exists(".get_palette")) {
    .get_palette(palette_name, 3)
  } else {
    c("#1b9e77", "#d95f02", "#7570b3")
  }
  
  p <- ggplot2::ggplot(ext_perf, ggplot2::aes(x = Generation, y = Performance)) +
    ggplot2::geom_line(colour = cols[1], linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(size = NumFeatures),
                        colour = cols[2], alpha = 0.6) +
    ggplot2::scale_size_continuous(name = "# Features", range = c(2, 6)) +
    ggplot2::labs(
      title = "GA Evolution - Genetic Algorithm Feature Selection",
      subtitle = sprintf("Converged to %d features after %d generations",
                         length(ga_obj$optVariables), n_rows),
      x = "Generation",
      y = metric_col
    ) +
    .pub_theme(base_size)+
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.title    = ggplot2::element_text(face = "bold")
    )
  
  # Add optimal point if visible
  if (opt_iter <= n_rows) {
    opt_row <- ext_perf[opt_iter, ]
    p <- p +
      ggplot2::geom_point(
        data = opt_row,
        mapping = ggplot2::aes(x = Generation, y = Performance),
        colour = cols[3], size = 5, shape = 21, fill = cols[3], stroke = 1.5,
        inherit.aes = FALSE
      )
  }
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, paste0("GA_evolution.", format))
    ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
    cat("Plot saved to:", path, "\n")
  }
  
  return(p)
}


#' GA Feature Selection Frequency
#'
#' Heatmap showing how often each feature was selected across generations.
#'
#' @param ga_result A \code{gafs} object or result list.
#' @param top_n Show top N features. Default 30.
#' @param palette_name Palette. Default "Zissou1".
#' @param base_size Font size. Default 11.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 8 x 10.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotGAFrequency <- function(ga_result,
                            top_n        = 30,
                            palette_name = "Zissou1",
                            base_size    = 11,
                            save_plot    = FALSE,
                            save_dir     = NULL,
                            width        = 8,
                            height       = 10,
                            format       = "pdf") {
  
  cat("Generating GA selection frequency heatmap...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (is.list(ga_result) && "result" %in% names(ga_result)) {
    ga_obj <- ga_result$result
  } else {
    ga_obj <- ga_result
  }
  
  # Extract feature sets across generations
  all_sets <- ga_obj$ga$final
  all_features <- unique(unlist(all_sets))
  
  # Build binary matrix: generation x feature
  gen_mat <- matrix(0, nrow = length(all_sets), ncol = length(all_features))
  colnames(gen_mat) <- all_features
  rownames(gen_mat) <- paste0("Gen", seq_along(all_sets))
  
  for (i in seq_along(all_sets)) {
    gen_mat[i, all_sets[[i]]] <- 1
  }
  
  # Calculate selection frequency
  freq <- colSums(gen_mat) / nrow(gen_mat)
  freq_df <- data.frame(
    Feature = names(freq),
    Frequency = freq,
    row.names = NULL
  ) %>%
    dplyr::arrange(desc(Frequency)) %>%
    head(top_n)
  
  # Subset matrix to top features
  gen_mat_top <- gen_mat[, freq_df$Feature, drop = FALSE]
  
  # Convert to long format
  heatmap_df <- reshape2::melt(gen_mat_top, 
                               varnames = c("Generation", "Feature"),
                               value.name = "Selected")
  heatmap_df$Selected <- factor(heatmap_df$Selected, levels = c(0, 1))
  
  cols <- .get_palette(palette_name, 2)
  
  p <- ggplot2::ggplot(heatmap_df, 
                       ggplot2::aes(x = Generation, y = Feature, fill = Selected)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      values = c("0" = "grey90", "1" = cols[1]),
      labels = c("0" = "Not Selected", "1" = "Selected")
    ) +
    ggplot2::labs(
      title = "GA Feature Selection Across Generations",
      subtitle = sprintf("Top %d features by selection frequency", nrow(freq_df)),
      x = "Generation",
      y = NULL,
      fill = NULL
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold", colour = "black"),
      legend.position = "bottom"
    )
  
  if (save_plot) {
    .save_plot(p, save_dir, "GA_frequency", width, height, format)
  }
  
  return(p)
}


## =============================================================================
##   3  SA -- Simulated Annealing
## =============================================================================

#' SA Acceptance Trace Plot
#'
#' Shows fitness evolution and acceptance probability across SA iterations.
#'
#' @param sa_result A \code{safs} object or result list from
#'   \code{FeatureSelectSA}.
#' @param metric Fitness metric. Default "Fitness".
#' @param palette_name Palette. Default "Cavalcanti1".
#' @param base_size Font size. Default 13.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 8 x 6.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotSA <- function(sa_result,
                   metric       = "Fitness",
                   palette_name = "Cavalcanti1",
                   base_size    = 13,
                   save_plot    = FALSE,
                   save_dir     = NULL,
                   width        = 8,
                   height       = 6,
                   format       = "pdf") {
  
  cat("Generating SA acceptance trace...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (is.list(sa_result) && "result" %in% names(sa_result)) {
    sa_obj <- sa_result$result
  } else if (inherits(sa_result, "safs")) {
    sa_obj <- sa_result
  } else {
    stop("Input must be a 'safs' object or result from FeatureSelectSA().")
  }
  
  # Extract external performance
  ext_perf <- sa_obj$external
  
  if (is.null(ext_perf)) {
    stop("No external performance data in SA result.")
  }
  
  metric_col <- grep(metric, colnames(ext_perf), value = TRUE)[1]
  if (is.na(metric_col)) {
    stop("Metric '", metric, "' not found in external performance.")
  }
  
  ext_perf$Iteration <- seq_len(nrow(ext_perf))
  ext_perf$Performance <- ext_perf[[metric_col]]
  
  # Get number of features
  n_vars <- sapply(sa_obj$sa$currentSet, length)
  ext_perf$NumFeatures <- n_vars
  
  opt_iter <- which.max(ext_perf$Performance)
  
  cols <- .get_palette(palette_name, 3)
  
  p <- ggplot2::ggplot(ext_perf, 
                       ggplot2::aes(x = Iteration, y = Performance)) +
    ggplot2::geom_line(colour = cols[1], linewidth = 0.8, alpha = 0.7) +
    ggplot2::geom_point(
      ggplot2::aes(size = NumFeatures),
      colour = cols[2],
      alpha = 0.5
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = FALSE,
      colour = cols[3],
      linewidth = 1.2,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      data = ext_perf[opt_iter, ],
      colour = "red",
      size = 5,
      shape = 21,
      fill = "red",
      stroke = 1.5
    ) +
    ggplot2::scale_size_continuous(
      name = "# Features",
      range = c(1.5, 5)
    ) +
    ggplot2::labs(
      title = "SA Trace - Simulated Annealing Feature Selection",
      subtitle = sprintf("Optimized to %d features after %d iterations",
                        length(sa_obj$optVariables), nrow(ext_perf)),
      x = "Iteration",
      y = metric_col
    ) +
    .pub_theme(base_size)
  
  if (save_plot) {
    .save_plot(p, save_dir, "SA_trace", width, height, format)
  }
  
  return(p)
}


## =============================================================================
##   4  SBF -- Selection By Filtering
## =============================================================================

#' SBF Score Distribution Plot
#'
#' Violin + box plot showing score distribution for selected vs removed features.
#'
#' @param sbf_result An \code{sbf} object or result list from
#'   \code{FeatureSelectSBF}.
#' @param palette_name Palette. Default "Darjeeling2".
#' @param base_size Font size. Default 13.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 7 x 6.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotSBF <- function(sbf_result,
                    palette_name = "Darjeeling2",
                    base_size    = 13,
                    save_plot    = FALSE,
                    save_dir     = NULL,
                    width        = 7,
                    height       = 6,
                    format       = "pdf") {
  
  cat("Generating SBF score distribution...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (is.list(sbf_result) && "result" %in% names(sbf_result)) {
    sbf_obj <- sbf_result$result
  } else if (inherits(sbf_result, "sbf")) {
    sbf_obj <- sbf_result
  } else {
    stop("Input must be an 'sbf' object or result from FeatureSelectSBF().")
  }
  
  # Extract variable importance/scores
  var_scores <- sbf_obj$variables
  
  if (is.null(var_scores)) {
    stop("No variable scores in SBF result.")
  }
  
  # Get selected variables
  selected_vars <- sbf_obj$optVariables
  
  # Aggregate scores across resamples
  score_df <- var_scores %>%
    dplyr::group_by(var) %>%
    dplyr::summarise(
      Score = mean(Overall, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Status = ifelse(var %in% selected_vars, "Selected", "Removed")
    )
  
  cols <- .get_palette(palette_name, 2)
  
  p <- ggplot2::ggplot(score_df, 
                       ggplot2::aes(x = Status, y = Score, fill = Status)) +
    ggplot2::geom_violin(alpha = 0.5, colour = NA) +
    ggplot2::geom_boxplot(
      width = 0.2,
      outlier.shape = 21,
      outlier.size = 2,
      colour = "grey30"
    ) +
    ggplot2::geom_jitter(
      width = 0.1,
      alpha = 0.3,
      size = 1.5,
      colour = "grey20"
    ) +
    ggplot2::scale_fill_manual(values = c("Selected" = cols[1], "Removed" = cols[2])) +
    ggplot2::labs(
      title = "SBF Score Distribution - Selection By Filtering",
      subtitle = sprintf("%d selected, %d removed",
                        sum(score_df$Status == "Selected"),
                        sum(score_df$Status == "Removed")),
      x = NULL,
      y = "Filter Score",
      fill = "Status"
    ) +
    ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      size = 5
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = "none")
  
  if (save_plot) {
    .save_plot(p, save_dir, "SBF_distribution", width, height, format)
  }
  
  return(p)
}


## =============================================================================
##   5  Built-in Variable Importance
## =============================================================================

#' Built-in Variable Importance Comparison
#'
#' Heatmap showing which features were selected by each model.
#'
#' @param builtin_result Result list from \code{FeatureSelectBuiltin}.
#' @param palette_name Palette. Default "Royal1".
#' @param base_size Font size. Default 11.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 8 x 10.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotBuiltinImportance <- function(builtin_result,
                                  palette_name = "Royal1",
                                  base_size    = 11,
                                  save_plot    = FALSE,
                                  save_dir     = NULL,
                                  width        = 8,
                                  height       = 10,
                                  format       = "pdf") {
  
  cat("Generating built-in importance heatmap...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (!is.list(builtin_result) || !"importance_table" %in% names(builtin_result)) {
    stop("Input must be result from FeatureSelectBuiltin().")
  }
  
  imp_table <- builtin_result$importance_table
  
  # Convert to long format for heatmap
  model_cols <- setdiff(colnames(imp_table), c("Feature", "Selected"))
  
  long_df <- imp_table %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(model_cols),
      names_to = "Model",
      values_to = "Selected_Status"
    ) %>%
    dplyr::mutate(
      Value = ifelse(Selected_Status == "Yes", 1, 0)
    )
  
  # Order features by total selection count
  feature_order <- imp_table %>%
    dplyr::arrange(desc(rowSums(.[, model_cols] == "Yes"))) %>%
    dplyr::pull(Feature)
  
  long_df$Feature <- factor(long_df$Feature, levels = rev(feature_order))
  
  cols <- .get_palette(palette_name, 2)
  
  p <- ggplot2::ggplot(long_df, 
                       ggplot2::aes(x = Model, y = Feature, fill = factor(Value))) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = c("0" = "grey90", "1" = cols[1]),
      labels = c("0" = "Not Selected", "1" = "Selected"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Built-in Variable Importance Across Models",
      subtitle = sprintf("%d features, %d models compared",
                        length(unique(long_df$Feature)),
                        length(unique(long_df$Model))),
      x = NULL,
      y = NULL
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold", colour = "black"),
      legend.position = "bottom"
    )
  
  if (save_plot) {
    .save_plot(p, save_dir, "builtin_importance", width, height, format)
  }
  
  return(p)
}


## =============================================================================
##   6  INTEGRATED -- Multi-method Comparison & Consensus
## =============================================================================

#' UpSet Plot for Feature Selection Methods
#'
#' Shows intersections of features selected by different methods.
#'
#' @param multi_result Result list from \code{FeatureSelectMulti} containing
#'   \code{results} and \code{selected_features}.
#' @param palette_name Palette. Default "GrandBudapest1".
#' @param base_size Font size. Default 13.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 10 x 7.
#' @param format Format. Default "pdf".
#' @return NULL (displays UpSetR plot directly).
#' @export
PlotFeatureUpSet <- function(multi_result,
                             palette_name = "GrandBudapest1",
                             base_size    = 13,
                             save_plot    = FALSE,
                             save_dir     = NULL,
                             width        = 10,
                             height       = 7,
                             format       = "pdf") {
  
  cat("Generating UpSet plot...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Package 'UpSetR' required. Install with: install.packages('UpSetR')")
  }
  
  if (!is.list(multi_result) || !"results" %in% names(multi_result)) {
    stop("Input must be result from FeatureSelectMulti().")
  }
  
  # Extract feature lists
  feature_lists <- lapply(multi_result$results, function(x) x$opt_vars)
  
  # Create binary matrix
  upset_mat <- UpSetR::fromList(feature_lists)
  
  cols <- .get_palette(palette_name, length(feature_lists))
  
  # Create upset plot
  p <- UpSetR::upset(
    upset_mat,
    sets = names(feature_lists),
    order.by = "freq",
    decreasing = TRUE,
    text.scale = c(1.3, 1.2, 1.1, 1, 1.5, 1.2),
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Features Selected",
    point.size = 3.5,
    line.size = 1.2,
    mb.ratio = c(0.6, 0.4),
    sets.bar.color = cols,
    main.bar.color = "grey30",
    matrix.color = "grey30"
  )
  
  if (save_plot) {
    pdf(file.path(save_dir, paste0("feature_upset.", format)),
        width = width, height = height)
    print(p)
    dev.off()
    cat("Plot saved to:", file.path(save_dir, paste0("feature_upset.", format)), "\n")
  } else {
    print(p)
  }
  
  invisible(p)
}


#' Feature Selection Consensus Heatmap
#'
#' Heatmap showing which features were selected by which methods, ordered
#' by consensus score.
#'
#' @param multi_result Result from \code{FeatureSelectMulti}.
#' @param palette_name Palette. Default "Moonrise3".
#' @param show_counts Annotate cells with selection counts. Default FALSE.
#' @param base_size Font size. Default 11.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 8 x 12.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotFeatureConsensus <- function(multi_result,
                                 palette_name = "Moonrise3",
                                 show_counts  = FALSE,
                                 base_size    = 11,
                                 save_plot    = FALSE,
                                 save_dir     = NULL,
                                 width        = 8,
                                 height       = 12,
                                 format       = "pdf") {
  
  cat("Generating consensus heatmap...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (!is.list(multi_result) || !"results" %in% names(multi_result)) {
    stop("Input must be result from FeatureSelectMulti().")
  }
  
  # Extract feature lists
  feature_lists <- lapply(multi_result$results, function(x) x$opt_vars)
  all_features <- unique(unlist(feature_lists))
  
  # Build binary matrix
  consensus_mat <- matrix(0, nrow = length(all_features), ncol = length(feature_lists))
  rownames(consensus_mat) <- all_features
  colnames(consensus_mat) <- names(feature_lists)
  
  for (i in seq_along(feature_lists)) {
    consensus_mat[feature_lists[[i]], i] <- 1
  }
  
  # Calculate consensus score (proportion of methods selecting each feature)
  consensus_score <- rowSums(consensus_mat) / ncol(consensus_mat)
  consensus_mat <- consensus_mat[order(-consensus_score), , drop = FALSE]
  
  # Convert to long format
  heatmap_df <- reshape2::melt(
    consensus_mat,
    varnames = c("Feature", "Method"),
    value.name = "Selected"
  )
  
  # Add consensus score
  heatmap_df$Consensus <- consensus_score[as.character(heatmap_df$Feature)]
  heatmap_df$Selected_Factor <- factor(heatmap_df$Selected, levels = c(0, 1))
  
  # Order features by consensus
  feature_order <- rownames(consensus_mat)
  heatmap_df$Feature <- factor(heatmap_df$Feature, levels = rev(feature_order))
  
  cols <- .get_palette(palette_name, 2)
  
  p <- ggplot2::ggplot(heatmap_df, 
                       ggplot2::aes(x = Method, y = Feature, fill = Selected_Factor)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3)
  
  if (show_counts) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = Selected),
      colour = "white",
      fontface = "bold",
      size = 3
    )
  }
  
  p <- p +
    ggplot2::scale_fill_manual(
      values = c("0" = "grey90", "1" = cols[1]),
      labels = c("0" = "Not Selected", "1" = "Selected"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Feature Selection Consensus Across Methods",
      subtitle = sprintf("%d features, %d methods", 
                        nrow(consensus_mat), ncol(consensus_mat)),
      x = NULL,
      y = NULL
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold", colour = "black", size = 8),
      legend.position = "bottom"
    )
  
  if (save_plot) {
    .save_plot(p, save_dir, "feature_consensus", width, height, format)
  }
  
  return(p)
}


#' Feature Stability Barplot
#'
#' Bar chart showing how many methods selected each feature, sorted by
#' stability (consensus).
#'
#' @param multi_result Result from \code{FeatureSelectMulti}.
#' @param top_n Show top N features. NULL = all. Default NULL.
#' @param palette_name Palette. Default "Darjeeling1".
#' @param base_size Font size. Default 13.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 8 x 10.
#' @param format Format. Default "pdf".
#' @return A ggplot object.
#' @export
PlotFeatureStability <- function(multi_result,
                                 top_n        = NULL,
                                 palette_name = "BottleRocket2",
                                 base_size    = 13,
                                 save_plot    = FALSE,
                                 save_dir     = NULL,
                                 width        = 8,
                                 height       = 10,
                                 format       = "pdf") {
  
  cat("Generating feature stability plot...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (!is.list(multi_result) || !"results" %in% names(multi_result)) {
    stop("Input must be result from FeatureSelectMulti().")
  }
  
  # Extract feature lists
  feature_lists <- lapply(multi_result$results, function(x) x$opt_vars)
  all_features <- unique(unlist(feature_lists))
  n_methods <- length(feature_lists)
  
  # Count selections
  selection_counts <- sapply(all_features, function(feat) {
    sum(sapply(feature_lists, function(fl) feat %in% fl))
  })
  
  stability_df <- data.frame(
    Feature = names(selection_counts),
    Count = as.numeric(selection_counts),
    Proportion = as.numeric(selection_counts) / n_methods,
    row.names = NULL
  ) %>%
    dplyr::arrange(desc(Count), Feature)
  
  if (!is.null(top_n)) {
    stability_df <- head(stability_df, top_n)
  }
  
  stability_df$Feature <- factor(stability_df$Feature, 
                                 levels = rev(stability_df$Feature))
  
  # Color gradient based on stability
  cols <- .get_palette(palette_name, 5)
  
  p <- ggplot2::ggplot(stability_df, 
                       ggplot2::aes(x = Count, y = Feature, fill = Proportion)) +
    ggplot2::geom_col(colour = "grey30", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(
      low = cols[1], 
      high = cols[5],
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::percent,
      name = "Stability"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      hjust = -0.3,
      colour = "black",
      fontface = "bold",
      size = 3.5
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15)),
      breaks = seq(0, n_methods, 1)
    ) +
    ggplot2::labs(
      title = "Feature Selection Stability",
      subtitle = sprintf("Number of methods selecting each feature (out of %d)", 
                        n_methods),
      x = "Number of Methods",
      y = NULL
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "bold", colour = "black"),
      panel.grid.major.x = ggplot2::element_line(colour = "grey90", linewidth = 0.3)
    )
  
  if (save_plot) {
    .save_plot(p, save_dir, "feature_stability", width, height, format)
  }
  
  return(p)
}


#' Feature Selection Methods Comparison Summary
#'
#' Composite plot showing: (1) number of features selected, 
#' (2) performance metrics, (3) overlap statistics.
#'
#' @param multi_result Result from \code{FeatureSelectMulti}.
#' @param metric Performance metric to compare. Default "Accuracy".
#' @param palette_name Palette. Default "Royal2".
#' @param base_size Font size. Default 12.
#' @param save_plot Save? Default FALSE.
#' @param save_dir Output dir.
#' @param width,height Dimensions. Defaults 12 x 8.
#' @param format Format. Default "pdf".
#' @return A patchwork composite plot.
#' @export
PlotFeatureComparison <- function(multi_result,
                                  metric       = "Accuracy",
                                  palette_name = "Royal2",
                                  base_size    = 12,
                                  save_plot    = FALSE,
                                  save_dir     = NULL,
                                  width        = 12,
                                  height       = 8,
                                  format       = "pdf") {
  
  cat("Generating method comparison plot...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required. Install with: install.packages('patchwork')")
  }
  
  if (!is.list(multi_result) || !"results" %in% names(multi_result)) {
    stop("Input must be result from FeatureSelectMulti().")
  }
  
  results <- multi_result$results
  feature_lists <- lapply(results, function(x) x$opt_vars)
  
  # Panel A: Number of features
  n_features <- data.frame(
    Method = names(feature_lists),
    Count = sapply(feature_lists, length),
    row.names = NULL
  )
  n_features$Method <- factor(n_features$Method, levels = n_features$Method)
  
  cols <- .get_palette(palette_name, length(feature_lists))
  
  p_count <- ggplot2::ggplot(n_features, 
                             ggplot2::aes(x = Method, y = Count, fill = Method)) +
    ggplot2::geom_col(colour = "grey30", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      vjust = -0.5,
      fontface = "bold",
      size = 4
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      title = "Number of Features Selected",
      x = NULL,
      y = "Count"
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  # Panel B: Performance (if available)
  perf_data <- lapply(names(results), function(method) {
    res <- results[[method]]$result
    # Try to extract best performance
    if ("results" %in% names(res)) {
      perf <- max(res$results[[metric]], na.rm = TRUE)
    } else if ("external" %in% names(res)) {
      metric_col <- grep(metric, colnames(res$external), value = TRUE)[1]
      if (!is.na(metric_col)) {
        perf <- max(res$external[[metric_col]], na.rm = TRUE)
      } else {
        perf <- NA
      }
    } else {
      perf <- NA
    }
    data.frame(Method = method, Performance = perf)
  })
  perf_df <- dplyr::bind_rows(perf_data)
  perf_df <- perf_df[!is.na(perf_df$Performance), ]
  
  if (nrow(perf_df) > 0) {
    perf_df$Method <- factor(perf_df$Method, levels = names(feature_lists))
    
    p_perf <- ggplot2::ggplot(perf_df, 
                              ggplot2::aes(x = Method, y = Performance, fill = Method)) +
      ggplot2::geom_col(colour = "grey30", linewidth = 0.3) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.3f", Performance)),
        vjust = -0.5,
        fontface = "bold",
        size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = cols[perf_df$Method]) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
      ggplot2::labs(
        title = paste("Best", metric),
        x = NULL,
        y = metric
      ) +
      .pub_theme(base_size) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  } else {
    p_perf <- ggplot2::ggplot() + 
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = "Performance data unavailable",
                       size = 6) +
      ggplot2::theme_void()
  }
  
  # Panel C: Overlap matrix
  n_methods <- length(feature_lists)
  overlap_mat <- matrix(0, nrow = n_methods, ncol = n_methods)
  rownames(overlap_mat) <- colnames(overlap_mat) <- names(feature_lists)
  
  for (i in seq_len(n_methods)) {
    for (j in seq_len(n_methods)) {
      overlap_mat[i, j] <- length(intersect(feature_lists[[i]], feature_lists[[j]]))
    }
  }
  
  overlap_df <- reshape2::melt(overlap_mat, varnames = c("Method1", "Method2"),
                               value.name = "Overlap")
  
  p_overlap <- ggplot2::ggplot(overlap_df, 
                               ggplot2::aes(x = Method1, y = Method2, fill = Overlap)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = Overlap),
      colour = "black",
      fontface = "bold",
      size = 3.5
    ) +
    ggplot2::scale_fill_gradient(
      low = "white", 
      high = cols[1],
      name = "# Shared"
    ) +
    ggplot2::labs(
      title = "Feature Overlap Matrix",
      x = NULL,
      y = NULL
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # Combine with patchwork
  composite <- (p_count | p_perf) / p_overlap +
    patchwork::plot_annotation(
      title = "Feature Selection Methods Comparison",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = base_size + 3
      ))
    )
  
  if (save_plot) {
    .save_plot(composite, save_dir, "method_comparison", width, height, format)
  }
  
  return(composite)
}

