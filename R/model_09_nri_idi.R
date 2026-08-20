# =============================================================================
# model_nri_idi.R
# NRI & IDI Analysis Module -- Independent from clinical_analysis.R
# =============================================================================
# Adapted from: Leong, L.T. (2021) "Area Under the Curve and Beyond"
# https://medium.com/data-science/area-under-the-curve-and-beyond-f87a8ec6937b
#
# Supports: Train_Model, caret train, ensemble, or raw probabilities
# Visual style: ggprism + wesanderson palettes

# 0. Package check -----------------------------------------------------------
#' @keywords internal
.check_nri_pkgs <- function() {
  required <- c("ggplot2", "dplyr", "tidyr", "wesanderson", "ggprism",
                "pROC", "scales", "ggrepel", "cowplot", "viridis")
  missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Please install them.")
  }
  invisible(TRUE)
}

# -- Internal helpers -------------------------------------------------------
.pub_theme <- function(base_size = 13) {
  ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            size = base_size + 2),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.title    = ggplot2::element_text(face = "bold"),
      legend.title  = ggplot2::element_text(face = "bold")
    )
}

.get_palette <- function(palette_name, n) {
  tryCatch(
    as.character(wesanderson::wes_palette(n = n, name = palette_name,
                                          type = if (n > 5) "continuous" else "discrete")),
    error = function(e) RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)]
  )
}

.save_plot <- function(p, dir, filename, width, height) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, paste0(tools::file_path_sans_ext(filename), ".pdf"))
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
  cat("Plot saved to:", path, "\n")
  invisible(path)
}

# =============================================================================
# Helper: Unified Probability and Outcome Extractor (Column-Safe)
# =============================================================================
#' @keywords internal
.extract_probs_and_truth <- function(model_obj,
                                     newdata        = NULL,
                                     model_name     = NULL,
                                     positive_class = NULL) {
  if (inherits(model_obj, "Train_Model")) {
    if (is.null(newdata)) {
      newdata <- model_obj@filtered.set$testing
      if (is.null(newdata)) newdata <- model_obj@split.data$testing
    }
    gc <- model_obj@group_col
    truth <- factor(newdata[[gc]])
    
    pos_cls <- if (!is.null(positive_class)) positive_class else levels(truth)[2]
    
    if (is.null(model_name)) {
      best <- model_obj@best.model.result$model
      if (is.null(best)) best <- model_obj@train.models[[1]]
      prob_mat <- predict(best, newdata, type = "prob")
    } else if (model_name == "ensemble") {
      ens <- model_obj@best.model.result$ensemble
      if (is.null(ens)) stop("No ensemble model found in model_obj.")
      return(list(
        truth    = truth,
        probs    = ens$predict_fn(newdata),
        positive = pos_cls,
        negative = setdiff(levels(truth), pos_cls)[1]
      ))
    } else {
      prob_mat <- predict(model_obj@train.models[[model_name]], newdata, type = "prob")
    }
    
  } else if (inherits(model_obj, "train")) {
    if (is.null(newdata)) stop("`newdata` argument is required for caret train object.")
    gc <- model_obj@group_col
    truth <- factor(newdata[[gc]])
    pos_cls <- if (!is.null(positive_class)) positive_class else levels(truth)[2]
    prob_mat <- predict(model_obj, newdata, type = "prob")
  } else {
    stop("model_obj must be a Train_Model or caret train object.")
  }
  
  # Robust column name matching including make.names check
  probs <- if (is.matrix(prob_mat) || is.data.frame(prob_mat)) {
    if (ncol(prob_mat) == 1) {
      prob_mat[, 1]
    } else if (!is.null(colnames(prob_mat)) && pos_cls %in% colnames(prob_mat)) {
      prob_mat[, pos_cls]
    } else if (!is.null(colnames(prob_mat)) && make.names(pos_cls) %in% colnames(prob_mat)) {
      prob_mat[, make.names(pos_cls)]
    } else {
      warning("Positive class '", pos_cls, "' not found in prediction columns. Defaulting to column 2.")
      prob_mat[, 2]
    }
  } else {
    prob_mat
  }
  
  neg_cls <- setdiff(levels(truth), pos_cls)[1]
  list(truth = truth, probs = probs, positive = pos_cls, negative = neg_cls)
}

# -- 1. Category-based NRI -------------------------------------------------
#' Category-Based Net Reclassification Index (with Per-Model Thresholds)
#'
#' Calculates NRI events, NRI non-events, and total NRI using possibly
#' different risk category thresholds for the reference and new model.
#'
#' @param truth           Binary outcome (factor or 0/1 vector).
#' @param ref_prob        Reference model predicted probabilities.
#' @param new_prob        New model predicted probabilities.
#' @param risk_thresholds Numeric vector of cut points applied to BOTH models
#'   when \code{ref_thresholds} and \code{new_thresholds} are NULL.
#' @param ref_thresholds  Optional separate thresholds for the reference model.
#' @param new_thresholds  Optional separate thresholds for the new model.
#' @return A list with NRI components and reclassification tables.
#' @export
CalculateCategoryNRI <- function(truth,
                                 ref_prob,
                                 new_prob,
                                 risk_thresholds = c(0.02, 0.1, 0.5, 0.95),
                                 ref_thresholds  = NULL,
                                 new_thresholds  = NULL) {
  .check_nri_pkgs()
  binary <- if (is.factor(truth)) as.numeric(truth == levels(truth)[2]) else truth
  
  if (is.null(ref_thresholds)) ref_thresholds <- risk_thresholds
  if (is.null(new_thresholds)) new_thresholds <- risk_thresholds
  
  n_ref_cats <- length(ref_thresholds) + 1
  n_new_cats <- length(new_thresholds) + 1
  
  .get_cat <- function(p, t, max_cat) {
    cat <- rep(1L, length(p))
    for (i in seq_along(t)) cat[p > t[i]] <- i + 1L
    factor(cat, levels = seq_len(max_cat))
  }
  
  ref_cat <- .get_cat(ref_prob, ref_thresholds, n_ref_cats)
  new_cat <- .get_cat(new_prob, new_thresholds, n_new_cats)
  
  ev_idx  <- which(binary == 1)
  nev_idx <- which(binary == 0)
  
  ref_num <- as.numeric(ref_cat)
  new_num <- as.numeric(new_cat)
  
  ev_up   <- sum(new_num[ev_idx]  > ref_num[ev_idx])
  ev_down <- sum(new_num[ev_idx]  < ref_num[ev_idx])
  nri_ev  <- if (length(ev_idx) > 0) (ev_up - ev_down) / length(ev_idx) else NA_real_
  
  nev_up   <- sum(new_num[nev_idx] > ref_num[nev_idx])
  nev_down <- sum(new_num[nev_idx] < ref_num[nev_idx])
  nri_nev  <- if (length(nev_idx) > 0) (nev_down - nev_up) / length(nev_idx) else NA_real_
  
  nri_tot <- if (!is.na(nri_ev) && !is.na(nri_nev)) nri_ev + nri_nev else NA_real_
  
  list(nri_events     = nri_ev,
       nri_nonevents  = nri_nev,
       nri_total      = nri_tot,
       events_up      = ev_up,
       events_down    = ev_down,
       nonevents_up   = nev_up,
       nonevents_down = nev_down,
       event_table    = table(ref_cat[ev_idx],  new_cat[ev_idx]),
       nonevent_table = table(ref_cat[nev_idx], new_cat[nev_idx]),
       ref_thresholds = ref_thresholds,
       new_thresholds = new_thresholds)
}

# -- 2. Bootstrap AUC & ROC -------------------------------------------------
# =============================================================================
# Helper: Bootstrap AUC Confidence Interval and Averaged ROC Curve
# =============================================================================
BootstrapROC <- function(truth, prob, n_boot = 500) {
  .check_nri_pkgs()
  
  n <- length(truth)
  
  # Define a fixed False Positive Rate (FPR) evaluation grid from 0 to 1
  base_fpr <- seq(0, 1, length.out = 101)
  tprs <- matrix(NA_real_, n_boot, 101)
  aucs <- numeric(n_boot)
  
  # Perform bootstrap resampling
  for (i in seq_len(n_boot)) {
    idx <- sample(n, replace = TRUE)
    
    # Ensure bootstrap sample contains at least two distinct class levels
    if (length(unique(truth[idx])) < 2) next
    
    roc_obj <- tryCatch(
      pROC::roc(truth[idx], prob[idx], direction = "auto", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) next
    
    aucs[i] <- as.numeric(pROC::auc(roc_obj))
    
    # Extract FPR (1 - Specificity) and TPR (Sensitivity)
    fpr <- 1 - roc_obj$specificities
    tpr <- roc_obj$sensitivities
    
    # Sort FPR in strictly ascending order (0 -> 1) to comply with approx() prerequisites
    ord <- order(fpr)
    fpr_sorted <- fpr[ord]
    tpr_sorted <- tpr[ord]
    
    # Remove duplicate FPR values to avoid interpolation ties
    non_dup <- !duplicated(fpr_sorted)
    
    # Interpolate True Positive Rate (TPR) over the fixed base FPR grid
    tprs[i, ] <- approx(x = fpr_sorted[non_dup],
                        y = tpr_sorted[non_dup],
                        xout = base_fpr,
                        rule = 2)$y
  }
  
  # Filter out failed bootstrap iterations
  valid_rows <- rowSums(!is.na(tprs)) > 0
  aucs_valid <- aucs[!is.na(aucs) & aucs > 0]
  tprs_valid <- tprs[valid_rows, , drop = FALSE]
  
  if (nrow(tprs_valid) == 0) {
    stop("Bootstrap evaluation failed across all iterations. Check sample size and class balance.")
  }
  
  list(
    thresholds = base_fpr,  # Added for 100% backward compatibility with PlotIDICurve
    fpr        = base_fpr,
    mean_fpr   = base_fpr,
    mean_tpr   = colMeans(tprs_valid, na.rm = TRUE),
    tpr_lower  = apply(tprs_valid, 2, quantile, 0.025, na.rm = TRUE),
    tpr_upper  = apply(tprs_valid, 2, quantile, 0.975, na.rm = TRUE),
    mean_auc   = mean(aucs_valid, na.rm = TRUE),
    auc_sd     = sd(aucs_valid, na.rm = TRUE)
  )
}

# -- 3. ROC curve comparison -----------------------------------------------
#' ROC Curve Comparison with Bootstrap CI
#'
#' Plots two ROC curves with optional 95% confidence ribbons, styled with
#' wesanderson palette and ggprism theme.
#'
#' @param truth      Binary outcome.
#' @param ref_prob   Reference model probabilities.
#' @param new_prob   New model probabilities.
#' @param labels     Character vector of length 2 (model names).
#' @param show_ci    Logical. Display bootstrap CI ribbons? Default TRUE.
#' @param n_boot     Bootstrap iterations.
#' @param save_plot  Logical.
#' @param save_dir   Output directory.
#' @return A ggplot object.
#' @export
PlotROCCompare <- function(truth,
                           ref_prob,
                           new_prob,
                           labels     = c("Reference", "New"),
                           show_ci    = TRUE,
                           n_boot     = 500,
                           save_plot  = FALSE,
                           save_dir   = NULL) {
  .check_nri_pkgs()
  ref_boot <- BootstrapROC(truth, ref_prob, n_boot)
  new_boot <- BootstrapROC(truth, new_prob, n_boot)
  
  auc_ref <- round(ref_boot$mean_auc, 3)
  auc_new <- round(new_boot$mean_auc, 3)
  ci_ref  <- round(1.96 * ref_boot$auc_sd, 3)
  ci_new  <- round(1.96 * new_boot$auc_sd, 3)
  
  df <- rbind(
    data.frame(fpr  = ref_boot$mean_fpr, tpr  = ref_boot$mean_tpr,
               lower = ref_boot$tpr_lower, upper = ref_boot$tpr_upper,
               Model = sprintf("%s (AUC=%.3f, CI=%.3f-%.3f)", labels[1], auc_ref, auc_ref - ci_ref, auc_ref + ci_ref)),
    data.frame(fpr  = new_boot$mean_fpr, tpr  = new_boot$mean_tpr,
               lower = new_boot$tpr_lower, upper = new_boot$tpr_upper,
               Model = sprintf("%s (AUC=%.3f, CI=%.3f-%.3f)", labels[2], auc_new, auc_new - ci_new, auc_new + ci_new))
  )
  
  cols <- .get_palette("Darjeeling1", 2)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = fpr, y = tpr, color = Model, fill = Model))
  
  if (show_ci) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA)
  }
  
  p <- p +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "ROC Curve Comparison with 95% CI",
                  x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") +
    ggplot2::coord_equal() +
    .pub_theme(13) +
    ggplot2::theme(legend.position = c(0.68, 0.18),
                   legend.background = ggplot2::element_rect(fill = "white", colour = "grey80"))
  
  print(p)
  if (save_plot) .save_plot(p, save_dir, "roc_comparison", 7, 6)
  invisible(p)
}

# -- 4. IDI curve ----------------------------------------------------------

#' Bootstrap Resampling for IDI Curve Interpolation
#'
#' @description
#' Internal helper function that performs bootstrap resampling to compute average
#' Sensitivity and 1-Specificity (False Positive Rate) profiles across a fixed grid
#' of risk cutoff thresholds (0 to 1) for IDI curve plotting.
#'
#' @param truth A numeric vector of binary outcomes (1 for cases/events, 0 for controls/non-events).
#' @param prob A numeric vector of predicted probabilities for the positive class.
#' @param n_boot Integer. Number of bootstrap iterations. Default is 500.
#'
#' @return A named list containing:
#' \item{cutoffs}{Numeric vector of 101 risk cutoff threshold points from 0 to 1.}
#' \item{mean_sens}{Numeric vector of bootstrapped mean Sensitivity values at each cutoff.}
#' \item{mean_spec}{Numeric vector of bootstrapped mean 1-Specificity (FPR) values at each cutoff.}
#'
#' @keywords internal
#' @noRd
BootstrapIDICurveData <- function(truth, prob, n_boot = 500) {
  .check_nri_pkgs()
  n <- length(truth)
  cutoff_grid <- seq(0, 1, length.out = 101)
  
  sens_mat <- matrix(NA_real_, n_boot, 101)
  spec_mat <- matrix(NA_real_, n_boot, 101) # Stores 1 - Specificity (FPR)
  
  for (i in seq_len(n_boot)) {
    idx <- sample(n, replace = TRUE)
    t_boot <- truth[idx]
    p_boot <- prob[idx]
    
    cases <- p_boot[t_boot == 1]
    ctrls <- p_boot[t_boot == 0]
    
    if (length(cases) == 0 || length(ctrls) == 0) next
    
    sens_mat[i, ] <- sapply(cutoff_grid, function(t) mean(cases >= t))
    spec_mat[i, ] <- sapply(cutoff_grid, function(t) mean(ctrls >= t))
  }
  
  valid_rows <- rowSums(!is.na(sens_mat)) > 0
  if (sum(valid_rows) == 0) {
    stop("Bootstrap IDI evaluation failed across all iterations. Check sample size or class balance.")
  }
  
  list(
    cutoffs   = cutoff_grid,
    mean_sens = colMeans(sens_mat[valid_rows, , drop = FALSE], na.rm = TRUE),
    mean_spec = colMeans(spec_mat[valid_rows, , drop = FALSE], na.rm = TRUE)
  )
}

#' Plot Integrated Discrimination Improvement (IDI) Curve
#'
#' @description
#' Calculates and visualizes the IDI, Integrated Sensitivity (IS), and Integrated Specificity (IP).
#'
#' @param truth A vector or factor containing binary outcomes.
#' @param ref_prob Predicted probabilities from the reference model.
#' @param new_prob Predicted probabilities from the new model.
#' @param positive The value in \code{truth} treated as the event. Defaults to the second unique value.
#' @param risk_thresholds Numeric vector for vertical reference lines.
#' @param n_boot Integer. Number of bootstrap iterations for curve smoothing.
#' @param colors Vector of length 2: \code{c(IP_color, IDI_label_fill)}.
#' @param sensitivity_color Color for Sensitivity lines and IS area.
#' @param linewidth Numeric. Width of lines.
#' @param alpha Numeric. Transparency for reference (dashed) lines.
#' @param ribbon_alpha Numeric. Transparency for shaded areas.
#' @param annotation_position List of \code{x, y, ip_y, idi_y} coordinates for labels.
#' @param annotation_size Numeric. Text size for labels.
#' @param title,subtitle,xlab,ylab,caption Plot labels and titles.
#' @param xlim,ylim,xbreaks,ybreaks Axis limits and breaks.
#' @param theme A \code{ggplot2} theme object.
#' @param save_plot Logical. Whether to save the plot.
#' @param save_dir Directory to save the plot.
#' @param save_filename Filename for the saved plot.
#' @param ... Additional arguments passed to \code{ggplot2::ggsave}.
#'
#' @return A list containing the \code{ggplot} object and calculated \code{is}, \code{ip}, and \code{idi}.
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. Generate synthetic data
#' set.seed(123)
#' labels <- factor(sample(c("Health", "Sick"), 200, replace = TRUE))
#' p_ref  <- runif(200, 0, 0.6)
#' p_new  <- pmin(p_ref + rnorm(200, 0.05, 0.03), 1) # Improved model
#' 
#' # 2. Run function (Assuming BootstrapROC is available in your environment)
#' res <- PlotIDICurve(truth = labels, ref_prob = p_ref, new_prob = p_new, positive = "Sick")
#' plot(res$plot)
#' }
PlotIDICurve <- function(truth,
                         ref_prob,
                         new_prob,
                         positive = NULL,
                         risk_thresholds = c(0.02, 0.1, 0.5, 0.95),
                         n_boot = 500,
                         colors = c("#E58601", "#46ACC8"),
                         sensitivity_color = "grey30",
                         linewidth = 1.1,
                         alpha = 0.6,
                         ribbon_alpha = 0.15,
                         annotation_position = list(x = 0.65, y = 0.85, ip_y = 0.15, idi_y = 0.5),
                         annotation_size = 4.5,
                         title = "Integrated Discrimination Improvement (IDI) Curve",
                         subtitle = "Solid: New Model | Dashed: Reference Model",
                         xlab = "Calculated Risk (Probability Threshold)",
                         ylab = "Sensitivity (Black) / 1-Specificity (Orange)",
                         caption = paste("IS: Improvement in mean probability for cases.",
                                         "\nIP: Improvement in mean probability for controls (reduction in 1-Spec)."),
                         xlim = c(0, 1),
                         ylim = c(0, 1),
                         xbreaks = seq(0, 1, 0.2),
                         ybreaks = seq(0, 1, 0.2),
                         theme = NULL,
                         save_plot = FALSE,
                         save_dir = NULL,
                         save_filename = "idi_curve.png",
                         ...) {
  
  if (is.factor(truth)) truth <- as.character(truth)
  u_vals <- unique(truth)
  if (length(u_vals) != 2) stop("truth must be binary.")
  
  if (is.null(positive)) {
    positive <- sort(u_vals, decreasing = TRUE)[1]
  }
  
  truth_num <- as.numeric(truth == positive)
  
  # NA-safe IDI Statistics
  is_val  <- mean(new_prob[truth_num == 1], na.rm = TRUE) - mean(ref_prob[truth_num == 1], na.rm = TRUE)
  ip_val  <- mean(ref_prob[truth_num == 0], na.rm = TRUE) - mean(new_prob[truth_num == 0], na.rm = TRUE)
  idi_val <- is_val + ip_val
  
  # Bootstrapped IDI Curves over Risk Cutoffs t
  ref_idi_data <- BootstrapIDICurveData(truth_num, ref_prob, n_boot)
  new_idi_data <- BootstrapIDICurveData(truth_num, new_prob, n_boot)
  
  df <- data.frame(
    threshold = new_idi_data$cutoffs,
    ref_sens  = ref_idi_data$mean_sens, new_sens = new_idi_data$mean_sens,
    ref_spec  = ref_idi_data$mean_spec, new_spec = new_idi_data$mean_spec
  )
  
  if (is.null(theme)) {
    theme <- ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     plot.caption = ggplot2::element_text(hjust = 0, size = 8, face = "italic"),
                     plot.title = ggplot2::element_text(face = "bold"))
  }
  
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_line(ggplot2::aes(x = threshold, y = ref_sens), color = sensitivity_color, linetype = "dashed", alpha = alpha) +
    ggplot2::geom_line(ggplot2::aes(x = threshold, y = new_sens), color = sensitivity_color, linewidth = linewidth) +
    ggplot2::geom_line(ggplot2::aes(x = threshold, y = ref_spec), color = colors[1], linetype = "dashed", alpha = alpha) +
    ggplot2::geom_line(ggplot2::aes(x = threshold, y = new_spec), color = colors[1], linewidth = linewidth) +
    ggplot2::geom_ribbon(ggplot2::aes(x = threshold, ymin = ref_sens, ymax = new_sens), fill = sensitivity_color, alpha = ribbon_alpha) +
    ggplot2::geom_ribbon(ggplot2::aes(x = threshold, ymin = new_spec, ymax = ref_spec), fill = colors[1], alpha = ribbon_alpha) +
    ggplot2::annotate("text", x = annotation_position$x, y = annotation_position$y, label = sprintf("IS = %.4f", is_val), fontface = "bold", hjust = 0) +
    ggplot2::annotate("text", x = annotation_position$x, y = annotation_position$ip_y, label = sprintf("IP = %.4f", ip_val), fontface = "bold", color = colors[1], hjust = 0) +
    ggplot2::annotate("label", x = annotation_position$x, y = annotation_position$idi_y, label = sprintf("IDI = %.4f", idi_val), color = "white", fill = colors[2], fontface = "bold") +
    ggplot2::scale_x_continuous(limits = xlim, breaks = xbreaks) +
    ggplot2::scale_y_continuous(limits = ylim, breaks = ybreaks) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab, caption = caption) + theme
  
  if (save_plot) {
    dir <- if (is.null(save_dir)) getwd() else save_dir
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    ggplot2::ggsave(file.path(dir, save_filename), p, width = 8, height = 6, ...)
  }
  
  return(list(plot = p, is = is_val, ip = ip_val, idi = idi_val))
}

# -- 5. NRI reclassification heatmap ---------------------------------------
#' NRI Reclassification Heatmap (Per-Model Labels, Auto-Match)
#'
#' Faceted heatmap showing reclassification counts for events and non-events.
#' Axis labels for risk categories can be customized independently. If the
#' provided labels do not match the number of categories, a warning is issued
#' and default labels are used.
#'
#' @param nri_result           Output of \code{CalculateCategoryNRI}.
#' @param ref_category_labels  Character vector of risk category names for
#'   the reference model (y-axis).
#' @param new_category_labels  Character vector for the new model (x-axis).
#' @param save_plot            Logical.
#' @param save_dir             Output directory.
#' @return A ggplot object.
#' @export
PlotNRIHeatmap <- function(nri_result,
                           ref_category_labels = NULL,
                           new_category_labels = NULL,
                           save_plot = FALSE,
                           save_dir = NULL) {
  .check_nri_pkgs()
  ev_df <- as.data.frame(nri_result$event_table)
  colnames(ev_df) <- c("Ref", "New", "Count")
  ev_df$Type <- "Events (Case)"
  
  nev_df <- as.data.frame(nri_result$nonevent_table)
  colnames(nev_df) <- c("Ref", "New", "Count")
  nev_df$Type <- "Non-Events (Control)"
  
  combined <- rbind(ev_df, nev_df)
  
  ref_ncat <- length(levels(combined$Ref))
  new_ncat <- length(levels(combined$New))
  
  # Safe label assignment with length validation
  if (!is.null(ref_category_labels)) {
    if (length(ref_category_labels) != ref_ncat) {
      warning("ref_category_labels length (", length(ref_category_labels),
              ") does not match categories (", ref_ncat, "). Falling back to defaults.")
      ref_category_labels <- paste0("R", 1:ref_ncat)
    }
  } else {
    ref_category_labels <- paste0("R", 1:ref_ncat)
  }
  
  if (!is.null(new_category_labels)) {
    if (length(new_category_labels) != new_ncat) {
      warning("new_category_labels length (", length(new_category_labels),
              ") does not match categories (", new_ncat, "). Falling back to defaults.")
      new_category_labels <- paste0("N", 1:new_ncat)
    }
  } else {
    new_category_labels <- paste0("N", 1:new_ncat)
  }
  
  levels(combined$Ref) <- ref_category_labels
  levels(combined$New) <- new_category_labels
  
  p <- ggplot2::ggplot(combined, ggplot2::aes(x = New, y = Ref, fill = Count)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = Count), colour = "white", fontface = "bold", size = 5) +
    ggplot2::facet_wrap(~ Type, ncol = 2) +
    ggplot2::scale_fill_viridis_c(option = "C", begin = 0.2, end = 0.9) +
    ggplot2::labs(
      title = "NRI Reclassification Heatmap",
      subtitle = sprintf("NRI Events = %.3f | NRI Non-Events = %.3f | Total NRI = %.3f",
                         nri_result$nri_events, nri_result$nri_nonevents, nri_result$nri_total),
      x = "New Model Risk Category", y = "Reference Model Risk Category", fill = "Count"
    ) +
    .pub_theme(13) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 12))
  
  print(p)
  if (save_plot) .save_plot(p, save_dir, "nri_heatmap", 9, 5)
  invisible(p)
}

# -- 6. NRI bar plot -------------------------------------------------------
#' NRI Reclassification Bar Plot
#'
#' Shows the number of patients reclassified correctly (green) or
#' incorrectly (red).
#'
#' @param nri_result Output of \code{CalculateCategoryNRI}.
#' @param save_plot  Logical.
#' @param save_dir   Output directory.
#' @return A ggplot object.
#' @export
PlotNRIBars <- function(nri_result, save_plot = FALSE, save_dir = NULL) {
  .check_nri_pkgs()
  df <- data.frame(
    Category  = c("Events Up", "Events Down", "Non-Events Up", "Non-Events Down"),
    Count     = c(nri_result$events_up, nri_result$events_down,
                  nri_result$nonevents_up, nri_result$nonevents_down),
    Type      = rep(c("Events", "Non-Events"), each = 2),
    Direction = rep(c("Upward", "Downward"), 2)
  )
  df$Fill <- ifelse(df$Direction == "Upward" & df$Type == "Events", "#06D6A0",
                    ifelse(df$Direction == "Downward" & df$Type == "Events", "#EF476F",
                           ifelse(df$Direction == "Upward" & df$Type == "Non-Events", "#EF476F", "#06D6A0")))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Category, y = Count, fill = Fill)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5, fontface = "bold", size = 4.5) +
    ggplot2::scale_fill_identity() +
    ggplot2::facet_wrap(~ Type, scales = "free_x") +
    ggplot2::labs(title = "NRI Reclassification Counts",
                  subtitle = "Green = Correct  |  Red = Incorrect",
                  y = "Number of Patients", x = NULL) +
    .pub_theme(13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
  
  print(p)
  if (save_plot) .save_plot(p, save_dir, "nri_bars", 8, 5)
  invisible(p)
}

# -- 7. Prediction probability distribution ---------------------------------
#' Prediction Probability Distribution by Outcome
#'
#' Violin + boxplot comparing predicted probabilities for reference and new
#' models, stratified by true outcome.  Axis labels are customizable.
#'
#' @param truth          Binary outcome.
#' @param ref_prob       Reference model probabilities.
#' @param new_prob       New model probabilities.
#' @param labels         Model names (length 2).
#' @param outcome_labels Labels for the outcome classes (length 2), e.g.
#'                       c("Benign", "Cancer").
#' @param save_plot      Logical.
#' @param save_dir       Output directory.
#' @return A ggplot object.
#' @export
PlotPredDist <- function(truth,
                         ref_prob,
                         new_prob,
                         labels    = c("Reference", "New"),
                         outcome_labels = c("Benign", "Malignant"),
                         save_plot = FALSE,
                         save_dir  = NULL) {
  .check_nri_pkgs()
  
  u_truth <- if (is.factor(truth)) levels(truth) else sort(unique(truth))
  formatted_truth <- if (length(outcome_labels) == length(u_truth)) {
    factor(truth, levels = u_truth, labels = outcome_labels)
  } else {
    factor(truth)
  }
  
  df <- data.frame(
    outcome   = formatted_truth,
    Reference = ref_prob,
    New       = new_prob
  )
  long_df <- tidyr::pivot_longer(df, -outcome, names_to = "Model", values_to = "Probability")
  cols <- .get_palette("Darjeeling1", 2)
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = outcome, y = Probability, fill = Model)) +
    ggplot2::geom_violin(alpha = 0.6, position = ggplot2::position_dodge(0.8)) +
    ggplot2::geom_boxplot(width = 0.2, position = ggplot2::position_dodge(0.8),
                          outlier.shape = NA, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "Prediction Probability Distribution by Outcome",
                  x = "True Outcome", y = "Predicted Probability") +
    .pub_theme(13) +
    ggplot2::theme(legend.position = "top")
  
  print(p)
  if (save_plot) .save_plot(p, save_dir, "pred_distribution", 8, 6)
  invisible(p)
}

# -- 8. Threshold-specific NRI curve ---------------------------------------
#' NRI as a Function of Risk Threshold
#'
#' Sweeps a single threshold from 0.01 to 0.99 and plots NRI components.
#'
#' @param truth     Binary outcome.
#' @param ref_prob  Reference model probabilities.
#' @param new_prob  New model probabilities.
#' @param save_plot Logical.
#' @param save_dir  Output directory.
#' @return A ggplot object.
#' @export
PlotThresholdNRI <- function(truth, ref_prob, new_prob,
                             save_plot = FALSE, save_dir = NULL) {
  .check_nri_pkgs()
  thresholds <- seq(0.01, 0.99, by = 0.01)
  res <- lapply(thresholds, function(t) {
    nri <- CalculateCategoryNRI(truth, ref_prob, new_prob, c(t))
    c(nri$nri_events, nri$nri_nonevents, nri$nri_total)
  })
  df <- data.frame(
    threshold    = rep(thresholds, 3),
    nri_value    = c(sapply(res, `[`, 1), sapply(res, `[`, 2), sapply(res, `[`, 3)),
    nri_type     = factor(rep(c("NRI Events", "NRI Non-Events", "Total NRI"), each = length(thresholds)),
                          levels = c("NRI Events", "NRI Non-Events", "Total NRI"))
  )
  
  cols <- c("NRI Events" = "#06D6A0", "NRI Non-Events" = "#E63946", "Total NRI" = "#4A4E69")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = nri_value, colour = nri_type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(title = "NRI as a Function of Risk Threshold",
                  x = "Risk Threshold", y = "NRI Value", colour = "Component") +
    .pub_theme(13) +
    ggplot2::theme(legend.position = "top")
  
  print(p)
  if (save_plot) .save_plot(p, save_dir, "threshold_nri", 8, 5)
  invisible(p)
}

# -- 9. Master NRI/IDI Analysis Pipeline -----------------------------------
#' Complete NRI/IDI Analysis Pipeline
#'
#' Runs the full suite of NRI/IDI analyses (ROC, IDI curve, NRI heatmap,
#' NRI bar plot, prediction distribution, threshold-NRI curve) for two
#' models.  Accepts Train_Model objects or raw probability vectors.
#'
#' @param model_obj_ref   Reference model (Train_Model, caret, or NULL).
#' @param model_obj_new   New model (Train_Model, caret, or NULL).
#' @param ref_prob        Raw reference probabilities (optional).
#' @param new_prob        Raw new probabilities (optional).
#' @param truth           Raw binary outcome (optional).
#' @param newdata         Common data frame for prediction.
#' @param positive The value in \code{truth} treated as the event. Defaults to the second unique value.
#' @param risk_thresholds Numeric vector of risk category boundaries.
#' @param ref_thresholds  Optional separate thresholds for the reference model.
#' @param new_thresholds  Optional separate thresholds for the new model.
#' @param labels          Model name labels (length 2).
#' @param outcome_labels  Labels for outcome classes (for violin plot).
#' @param ref_category_labels Labels for reference risk categories (for heatmap).
#' @param new_category_labels Labels for new model risk categories (for heatmap).
#' @param show_ci         Show bootstrap CI in ROC plot? Default TRUE.
#' @param n_boot          Bootstrap iterations.
#' @param save_plots      Save all plots?
#' @param save_dir        Output directory.
#' @return Invisible list with all results.
#' @export
NRI_IDI_Analysis <- function(model_obj_ref = NULL,
                             model_obj_new = NULL,
                             ref_prob       = NULL,
                             new_prob       = NULL,
                             truth          = NULL,
                             newdata        = NULL,
                             positive       = NULL,
                             risk_thresholds = c(0.02, 0.1, 0.5, 0.95),
                             ref_thresholds  = NULL,
                             new_thresholds  = NULL,
                             labels         = c("Reference", "New"),
                             outcome_labels = c("Benign", "Malignant"),
                             ref_category_labels = NULL,
                             new_category_labels = NULL,
                             show_ci        = TRUE,
                             n_boot         = 500,
                             save_plots     = TRUE,
                             save_dir       = "./NRIDI_Results/") {
  .check_nri_pkgs()
  if (is.null(save_dir) && save_plots) save_dir <- "./NRIDI_Results/"
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # Resolve inputs
  if (!is.null(model_obj_ref)) {
    res_ref <- .extract_probs_and_truth(model_obj_ref, newdata)
    ref_p   <- res_ref$probs
    truth   <- res_ref$truth
  } else if (!is.null(ref_prob) && !is.null(truth)) {
    ref_p <- ref_prob
  } else {
    stop("Provide either model_obj_ref or ref_prob + truth.")
  }
  
  if (!is.null(model_obj_new)) {
    res_new <- .extract_probs_and_truth(model_obj_new, newdata)
    new_p   <- res_new$probs
  } else if (!is.null(new_prob)) {
    new_p <- new_prob
  } else {
    stop("Provide either model_obj_new or new_prob.")
  }
  
  cat("--- NRI/IDI Analysis Pipeline ---\n")
  
  # 1. ROC comparison
  cat("[1/6] ROC comparison...\n")
  PlotROCCompare(truth, ref_p, new_p, labels, show_ci, n_boot, save_plots, save_dir)
  
  # 2. IDI curve
  cat("[2/6] IDI curve...\n")
  idi_res <- PlotIDICurve(
    truth = truth,
    ref_prob = ref_p,
    new_prob = new_p,
    positive = NULL,                
    risk_thresholds = risk_thresholds,
    n_boot = n_boot,
    save_plot = save_plots,
    save_dir = save_dir
  )
  cat(sprintf("  IS = %.4f | IP = %.4f | IDI = %.4f\n", idi_res$is, idi_res$ip, idi_res$idi))
  
  # 3. NRI
  cat("[3/6] Calculating NRI...\n")
  nri_res <- CalculateCategoryNRI(truth, ref_p, new_p,
                                  risk_thresholds = risk_thresholds,
                                  ref_thresholds  = ref_thresholds,
                                  new_thresholds  = new_thresholds)
  cat(sprintf("  NRI Events = %.4f | NRI Non-Events = %.4f | Total NRI = %.4f\n",
              nri_res$nri_events, nri_res$nri_nonevents, nri_res$nri_total))
  
  # 4. NRI heatmap
  cat("[4/6] NRI heatmap...\n")
  PlotNRIHeatmap(nri_res, ref_category_labels, new_category_labels, save_plots, save_dir)
  
  # 5. NRI bars
  cat("[5/6] NRI bar plot...\n")
  PlotNRIBars(nri_res, save_plots, save_dir)
  
  # 6. Prediction distribution
  cat("[6/6] Prediction distribution...\n")
  PlotPredDist(truth, ref_p, new_p, labels, outcome_labels, save_plots, save_dir)
  
  cat("--- NRI/IDI analysis complete. Plots saved to", save_dir, "---\n")
  
  invisible(list(nri = nri_res, idi = idi_res))
}
