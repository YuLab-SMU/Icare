# =============================================================================
# model_04_explain.R
# Model Explanation Module -- DALEX / SHAP / Beeswarm
#
# Interfaces:
#   Train_Model S4 object (module_01_object.R)
#   ExtractModel()         (module_03_compare.R)
#
# Functions:
#   CreateExplainer()          -- Create DALEX explainer
#   ExplainVariableImportance()-- Permutation importance (with error bars)
#   ExplainSHAP()              -- Single-observation SHAP waterfall
#   ExplainSHAPBeeswarm()      -- Global SHAP beeswarm *
#   ExplainBreakDown()         -- Break Down decomposition (with interactions)
#   ExplainCeterisParibus()    -- Ceteris Paribus what-if curves
#   ExplainPartialDependence() -- Global PDP / ALE
#   ExplainAll()               -- One-click full pipeline (including beeswarm)
#
# Dependencies: DALEX, iBreakDown, ingredients, ggplot2, shapviz (optional)
# =============================================================================

# -- 0. Package check ---------------------------------------------------------

.check_xai_packages <- function() {
  required <- c("DALEX", "iBreakDown", "ingredients", "ggplot2")
  missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "The following packages are missing. Please install them:\n",
      paste(" -", missing, collapse = "\n"), "\n\n",
      "install.packages(c('DALEX', 'iBreakDown', 'ingredients', 'ggplot2'))"
    )
  }
  invisible(TRUE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


# -- 1. CreateExplainer --------------------------------------------------------

#' Create a DALEX Explainer (Unified Entry for Single & Ensemble Models)
#'
#' Wraps a caret model, a \code{Train_Model} S4 object, or an ensemble
#' produced by \code{TrainEnsemble} into a DALEX explainer.  The argument
#' \code{model} controls which model is explained:
#' \itemize{
#'   \item \code{NULL} (default): the best model stored in the
#'         \code{Train_Model} object is used automatically.
#'   \item \code{"ensemble"}: the ensemble stored in
#'         \code{object@best.model.result$ensemble} is used.  The ensemble
#'         must contain a \code{predict_fn}.  Voting results are
#'         automatically converted to positive-class probabilities.
#'   \item A character string (e.g. \code{"rf"}, \code{"gbm"}): the
#'         corresponding model from \code{object@train.models} is used.
#'   \item A \code{caret::train} object: this model is used directly
#'         (requires \code{data} to be supplied when \code{object} is not
#'         a \code{Train_Model}).
#' }
#'
#' @param object    A \code{Train_Model} S4 object, or a \code{caret::train}
#'   object.
#' @param model     Which model to explain (see Description).  Ignored when
#'   \code{object} is a \code{caret::train} object.
#' @param data      Background data (data.frame).  Defaults to
#'   \code{object@filtered.set$training} for \code{Train_Model} objects.
#'   Required when \code{object} is a bare \code{caret::train} object.
#' @param y         Response vector (numeric 0/1).  If \code{NULL}, it is
#'   extracted from the \code{group_col} column of \code{data}.
#' @param group_col Name of the response column (default \code{"group"}).
#' @param label     Explainer label.  Auto-generated when \code{NULL}.
#' @param verbose   Print DALEX verbose output? Default \code{FALSE}.
#'
#' @return A \code{DALEX::explainer} object.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # ---- 1. Default: automatically use the best model from a Train_Model ----
#' explainer <- CreateExplainer(model_obj)
#'
#' # ---- 2. Explain a specific model by name ----
#' explainer_gbm <- CreateExplainer(model_obj, model = "gbm")
#'
#' # ---- 3. Explain a fine-tuned model (caret train object) ----
#' tuned_rf <- model_obj@best.model.result$fine_tuned_model
#' explainer_tuned <- CreateExplainer(model_obj, model = tuned_rf)
#'
#' # ---- 4. Explain a stacking ensemble ----
#' model_obj <- TrainEnsemble(model_obj, strategy = "stacking",
#'                            meta_method = "glm", top_n = 4)
#' explainer_stack <- CreateExplainer(model_obj, model = "ensemble")
#'
#' # ---- 5. Explain a voting ensemble ----
#' model_obj <- TrainEnsemble(model_obj, strategy = "voting", top_n = 4)
#' explainer_vote <- CreateExplainer(model_obj, model = "ensemble")
#'
#' # ---- 6. Explain a bare caret model (no Train_Model) ----
#' my_rf <- model_obj@train.models$rf
#' explainer_bare <- CreateExplainer(my_rf,
#'                                   data = model_obj@filtered.set$training,
#'                                   group_col = model_obj@group_col)
#'
#' # ---- Use the explainer with any interpretation function ----
#' ExplainVariableImportance(explainer_stack, top_n = 15,
#'                           save_plots = TRUE, save_dir = "./Explain/Ensemble/")
#' ExplainSHAPBeeswarm(explainer_stack, N = 50,
#'                     save_plots = TRUE, save_dir = "./Explain/Ensemble/")
#' }
CreateExplainer <- function(object,
                            model     = NULL,
                            data      = NULL,
                            y         = NULL,
                            group_col = NULL,
                            label     = NULL,
                            verbose   = FALSE) {
  
  .check_xai_packages()
  
  # -- 1. Resolve group_col and data --------------------------------------
  if (inherits(object, "Train_Model")) {
    if (is.null(group_col)) group_col <- object@group_col
    if (is.null(data)) {
      # Use scaled training data if available, else standard training data
      data <- object@split.scale.data$training %||% object@filtered.set$training
      if (is.null(data)) stop("No training data found in object. Provide 'data' manually.")
    }
  } else if (inherits(object, "train")) {
    if (is.null(group_col)) group_col <- "group"
    if (is.null(data)) stop("When passing a caret model, 'data' is required.")
  }
  
  # Ensure the target column is a factor
  if (!is.factor(data[[group_col]])) {
    data[[group_col]] <- as.factor(data[[group_col]])
  }
  
  # -- 2. Determine Positive Class (Handle numeric/X-prefixed levels) -----
  # DALEX typically needs probabilities for the 'positive' class (usually the 2nd level)
  raw_levels <- levels(data[[group_col]])
  positive_class <- raw_levels[2] 
  
  # -- 3. Resolve actual model to explain --------------------------------
  if (inherits(object, "Train_Model")) {
    if (is.null(model)) {
      # Automatically use fine_tuned_model if it exists, else best_model
      best_model <- object@best.model.result$fine_tuned_model %||% 
        object@train.models[[object@best.model.result$model_type]]
      if (is.null(label)) label <- object@best.model.result$model_type %||% "BestModel"
    } else if (is.character(model) && model == "ensemble") {
      ens <- object@best.model.result$ensemble
      if (is.null(ens)) stop("No ensemble found. Run TrainEnsemble first.")
      best_model <- object
      if (is.null(label)) label <- paste0("Ensemble (", ens$strategy, ")")
    } else if (is.character(model)) {
      best_model <- object@train.models[[model]]
      if (is.null(label)) label <- model
    } else if (inherits(model, "train")) {
      best_model <- model
      if (is.null(label)) label <- best_model$method
    }
  } else {
    best_model <- object
    if (is.null(label)) label <- object$method
  }
  
  # -- 4. Build Predict Function (Safe against numeric levels) ----------
  # This function ensures that newdata doesn't contain y, and extracts the correct prob col
  if (is.null(model) || !is.character(model) || model != "ensemble") {
    
    predict_fn <- function(mod, newdata) {
      # Remove target column if present
      X_input <- newdata[, setdiff(colnames(newdata), group_col), drop = FALSE]
      probs <- predict(mod, X_input, type = "prob")
      
      # Caret might have changed "0" to "X0" via make.names
      # We check which column name in 'probs' matches our positive_class
      col_idx <- which(colnames(probs) == positive_class)
      if (length(col_idx) == 0) {
        # Fallback: if names don't match, try make.names version
        col_idx <- which(colnames(probs) == make.names(positive_class))
      }
      
      if (length(col_idx) == 0) {
        # Final fallback: use the second column (standard for binary)
        return(probs[, 2])
      }
      return(probs[, col_idx])
    }
    
  } else {
    # Ensemble Logic
    predict_fn <- function(mod, newdata) {
      # The ensemble's stored predict_fn should already handle internal column dropping
      pred <- object@best.model.result$ensemble$predict_fn(newdata)
      
      if (is.numeric(pred)) return(pred)
      
      # If ensemble returns class labels, convert to numeric probability based on pos_class
      if (is.factor(pred) || is.character(pred)) {
        res <- as.numeric(make.names(as.character(pred)) == make.names(positive_class))
        return(res)
      }
      stop("Unsupported ensemble output type.")
    }
  }
  
  # -- 5. Prepare y (numeric 0/1) and X --------------------------------
  if (is.null(y)) {
    # DALEX works best when y is numeric 0 and 1
    # We ensure y matches the factor underlying levels
    y <- as.numeric(data[[group_col]]) - 1L
  }
  
  X <- data[, setdiff(colnames(data), group_col), drop = FALSE]
  
  # -- 6. Create Explainer ---------------------------------------------
  explainer <- DALEX::explain(
    model            = best_model,
    data             = X,
    y                = y,
    predict_function = predict_fn,
    label            = label %||% "model",
    verbose          = verbose
  )
  
  cat(sprintf("[OK] Explainer created | Label: '%s' | Pos Class: '%s' | %d features\n", 
              explainer$label, positive_class, ncol(X)))
  
  return(explainer)
}

# -- 2. ExplainVariableImportance ----------------------------------------------
#' Variable Importance (Permutation) - Universal for Any Caret Model (Fixed)
#'
#' @param explainer     A DALEX explainer.
#' @param B             Number of permutations (default 10).
#' @param type          Loss type: "difference", "ratio", "raw".
#' @param show_boxplots (Ignored; kept for compatibility)
#' @param top_n         Number of most important features to display.
#' @param filter_zero   Remove features with near-zero importance.
#' @param save_plots    Save plot? Default FALSE.
#' @param save_dir      Output directory.
#' @param plot_width,plot_height  Dimensions.
#' @param base_size     Base font size.
#'
#' @return Invisibly returns the DALEX model_parts object.
#' @export
ExplainVariableImportance <- function(explainer,
                                      B             = 10,
                                      type          = "difference",
                                      show_boxplots = TRUE,
                                      top_n         = 20,
                                      filter_zero   = TRUE,
                                      save_plots    = FALSE,
                                      save_dir      = NULL,
                                      plot_width    = 8,
                                      plot_height   = 5,
                                      base_size     = 13) {
  
  .check_xai_packages()
  
  cat("-- Variable Importance (Permutation) -------------------------------------\n")
  vi <- DALEX::model_parts(explainer, B = B, type = type)
  print(vi)
  
  # Keep only features (remove baseline and full model)
  vi_plot <- vi[!vi$variable %in% c("_baseline_", "_full_model_"), ]
  
  # Optionally filter zero-importance
  if (filter_zero) {
    mean_loss <- aggregate(vi_plot$dropout_loss, by = list(vi_plot$variable), FUN = mean)
    colnames(mean_loss) <- c("variable", "mean_loss")
    zero_vars <- mean_loss$variable[abs(mean_loss$mean_loss) < 1e-9]
    if (length(zero_vars) > 0) {
      cat("Excluding", length(zero_vars), "features with zero importance.\n")
      vi_plot <- vi_plot[!vi_plot$variable %in% zero_vars, ]
    }
  }
  
  # Calculate mean and SE per variable
  importance_summary <- vi_plot %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      mean_loss = mean(dropout_loss, na.rm = TRUE),
      se_loss   = sd(dropout_loss, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(abs(mean_loss))) %>%
    head(top_n)
  
  # Order factor levels
  importance_summary$variable <- factor(
    importance_summary$variable,
    levels = rev(importance_summary$variable)
  )
  
  # Build dot-and-whisker plot
  p <- ggplot2::ggplot(importance_summary, ggplot2::aes(x = mean_loss, y = variable)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", size = 0.3) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = mean_loss - se_loss, xmax = mean_loss + se_loss),
      height = 0.15, color = "black", size = 0.4, alpha = 0.6
    ) +
    ggplot2::geom_point(size = 2.5, color = "black", fill = "black", shape = 21, stroke = 0.3) +
    ggplot2::labs(
      title    = paste("Variable Importance -", explainer$label),
      subtitle = sprintf("Permutation importance (B = %d), mean +/- SE", B),
      x        = "Drop-out loss (1 - AUC)",
      y        = NULL
    ) +
    ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.text.y   = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  
  if (save_plots) {
    if (is.null(save_dir)) {
      if (exists(".get_viz_output_dir")) {
        save_dir <- .get_viz_output_dir("Model")
      } else {
        save_dir <- file.path(".", "Figures", "Model")
      }
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    filename <- file.path(save_dir, "variable_importance.pdf")
    ggplot2::ggsave(filename, plot = p, width = plot_width, height = plot_height,
                    device = "pdf", dpi = 300)
    cat("Plot saved to:", filename, "\n")
  }
  
  invisible(vi)
}


# -- 3. ExplainSHAP ------------------------------------------------------------

#' Single-Observation SHAP Waterfall (Custom Style, Flexible Filtering)
#'
#' Computes Shapley values and draws a publication-ready waterfall plot.
#' By default, variables with exactly zero SHAP contribution are hidden.
#'
#' @param explainer         A DALEX explainer.
#' @param new_observation   Single-row data frame or integer row index.
#' @param B                 Number of random permutations (default 25).
#' @param remove_zero       Logical; remove variables with zero SHAP? Default TRUE.
#' @param save_plots        Save plot? Default FALSE.
#' @param save_dir          Output directory.
#' @param plot_width,plot_height  Dimensions (inches).
#'
#' @return A ggplot object. Invisibly returns the DALEX predict_parts (shap) object.
#' @export
ExplainSHAP <- function(explainer,
                        new_observation,
                        B             = 25,
                        remove_zero   = TRUE,
                        save_plots    = FALSE,
                        save_dir      = "ModelExplain",
                        plot_width    = 8,
                        plot_height   = 5) {
  .check_xai_packages()
  
  obs  <- .resolve_observation(explainer, new_observation)
  pred <- explainer$predict_function(explainer$model, obs)
  
  cat("-- Single-Observation SHAP Waterfall ----------------------------\n")
  cat(sprintf("Predicted probability for this observation: %.4f\n", pred))
  
  shap_obj <- DALEX::predict_parts(
    explainer       = explainer,
    new_observation = obs,
    type            = "shap",
    B               = B
  )
  
  shap_df <- as.data.frame(shap_obj)
  # Remove internal rows
  shap_df <- shap_df[!shap_df$variable %in% c("_intercept_", "_baseline_"), ]
  if (B > 0) {
    shap_df <- shap_df[shap_df$B == 0, ]
  }
  
  # Optionally drop zero-importance features
  if (remove_zero) {
    shap_df <- shap_df[abs(shap_df$contribution) > 1e-12, ]
  }
  
  # Order by absolute contribution
  shap_df <- shap_df[order(abs(shap_df$contribution), decreasing = TRUE), ]
  shap_df$variable <- factor(shap_df$variable, levels = rev(shap_df$variable))
  
  cols <- wesanderson::wes_palette("Royal1", 2, type = "discrete")
  shap_df$fill_col <- ifelse(shap_df$contribution >= 0, cols[1], cols[2])
  
  p <- ggplot2::ggplot(shap_df, ggplot2::aes(x = contribution, y = variable, fill = fill_col)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.4f", contribution),
                   hjust = ifelse(contribution >= 0, -0.2, 1.2)),
      size = 3.5, colour = "black"
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(
      title    = paste0("SHAP Waterfall -- ", explainer$label),
      subtitle = sprintf("Shapley values (B = %d permutations)", B),
      x        = "SHAP contribution",
      y        = NULL
    ) +
    .pub_theme(13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.text.y   = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  if (save_plots) {
    if (is.null(save_dir) || !nzchar(save_dir)) save_dir <- "ModelExplain"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = file.path(save_dir, "shap_individual.pdf"), plot = p, width = plot_width, height = plot_height)
    cat("SHAP individual saved to:", file.path(save_dir, "shap_individual.pdf"), "\n")
  }
  
  invisible(shap_obj)
}


# -- 5. ExplainSHAPBeeswarm * -------------------------------------------------

#' Global SHAP Beeswarm
#'
#' Computes SHAP values for N samples and generates a beeswarm plot.
#'
#' Plot interpretation:
#' \itemize{
#'   \item y-axis: features, sorted by mean |SHAP| (most important at the top)
#'   \item x-axis: SHAP value (positive -> increases predicted probability,
#'         negative -> decreases it)
#'   \item colour: feature value (red = high, blue = low)
#'   \item each point represents one sample
#' }
#'
#' The function prefers the `shapviz` package for standard beeswarm plots;
#' if it is not installed, it falls back to a ggplot2 + geom_jitter
#' implementation (visually equivalent).
#'
#' @param explainer    A DALEX explainer.
#' @param N            Number of samples (default 100, suggest 50-200).
#' @param B            SHAP permutations per sample (default 10).
#' @param seed         Random seed (default 42).
#' @param max_features Maximum number of features to display (default 15).
#' @param save_plots   Save PDF? Default FALSE.
#' @param save_dir     Output directory (default "ModelExplain").
#' @param plot_width,plot_height  Dimensions (inches, default 9 x 6).
#'
#' @return Invisibly returns a named list:
#' \itemize{
#'   \item shap_matrix  -- N x p matrix of SHAP values
#'   \item X_sample     -- corresponding feature matrix
#'   \item shapviz_obj  -- shapviz object (if available, otherwise NULL)
#'   \item plot         -- ggplot2 object
#'   \item importance   -- data.frame of feature importance rankings
#' }
#' @export
ExplainSHAPBeeswarm <- function(explainer,
                                N            = 100,
                                B            = 10,
                                seed         = 42,
                                max_features = 15,
                                save_plots   = FALSE,
                                save_dir     = "ModelExplain",
                                plot_width   = 9,
                                plot_height  = 6) {
  .check_xai_packages()
  
  cat("-- Global SHAP Beeswarm ------------------------------------------\n")
  cat(sprintf("  Samples N = %d | Permutations per sample B = %d | Seed = %d\n", N, B, seed))
  cat("  Computing SHAP values, please wait...\n")
  
  set.seed(seed)
  X_full <- explainer$data
  n_use  <- min(N, nrow(X_full))
  idx    <- sample(nrow(X_full), n_use)
  X_samp <- X_full[idx, , drop = FALSE]
  
  feat_names <- colnames(X_samp)
  shap_mat   <- matrix(NA_real_, nrow = n_use, ncol = length(feat_names),
                       dimnames = list(NULL, feat_names))
  
  for (i in seq_len(n_use)) {
    if (i %% 10 == 0 || i == n_use)
      cat(sprintf("  Progress: %d / %d\r", i, n_use))
    
    sp <- tryCatch(
      DALEX::predict_parts(
        explainer       = explainer,
        new_observation = X_samp[i, , drop = FALSE],
        type            = "shap",
        B               = B
      ),
      error = function(e) NULL
    )
    
    if (!is.null(sp)) {
      avg     <- tapply(sp$contribution, sp$variable_name, mean)
      matched <- intersect(names(avg), feat_names)
      shap_mat[i, matched] <- avg[matched]
    }
  }
  cat("\n  SHAP matrix computed.\n")
  
  # Top features by mean absolute SHAP
  mean_abs  <- colMeans(abs(shap_mat), na.rm = TRUE)
  top_feats <- names(sort(mean_abs, decreasing = TRUE))[
    seq_len(min(max_features, length(mean_abs)))
  ]
  shap_top <- shap_mat[, top_feats, drop = FALSE]
  X_top    <- X_samp[, top_feats,  drop = FALSE]
  
  # Prefer shapviz if available
  sv_obj <- NULL
  p      <- NULL
  
  if (requireNamespace("shapviz", quietly = TRUE)) {
    cat("  Using shapviz for standard beeswarm...\n")
    sv_obj <- tryCatch(shapviz::shapviz(shap_top, X = X_top),
                       error = function(e) { message(e$message); NULL })
    if (!is.null(sv_obj)) {
      p <- tryCatch(
        shapviz::sv_importance(sv_obj,
                               kind         = "beeswarm",
                               show_numbers = FALSE,
                               max_display  = max_features) +
          ggplot2::labs(
            title    = paste0("SHAP Beeswarm -- ", explainer$label),
            subtitle = sprintf("N = %d samples, B = %d permutations | colour: feature value (red=high, blue=low)",
                               n_use, B),
            x        = "SHAP value (impact on predicted probability)",
            y        = NULL
          ) +
          ggprism::theme_prism(base_size = 13) +
          ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold"),
            legend.position = "right"
          ),
        error = function(e) { message(e$message); NULL }
      )
    }
  }
  
  # Fallback: ggplot2 manual beeswarm
  if (is.null(p)) {
    cat("  Using ggplot2 beeswarm (shapviz unavailable or failed)...\n")
    p <- .ggplot_beeswarm(shap_top, X_top,
                          model_label = explainer$label,
                          N = n_use, B = B)
  }
  
  print(p)
  if (save_plots) {
    if (is.null(save_dir) || !nzchar(save_dir)) save_dir <- "ModelExplain"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = file.path(save_dir, "shap_beeswarm.pdf"), plot = p, width = plot_width, height = plot_height)
    cat("SHAP beeswarm saved to:", file.path(save_dir, "shap_beeswarm.pdf"), "\n")
  }
  
  importance_df <- data.frame(
    Rank          = seq_along(top_feats),
    Feature       = top_feats,
    Mean_Abs_SHAP = round(mean_abs[top_feats], 6),
    row.names     = NULL
  )
  cat(sprintf("\n  Top %d features (by |mean SHAP|):\n", nrow(importance_df)))
  print(importance_df, row.names = FALSE)
  
  invisible(list(
    shap_matrix  = shap_mat,
    X_sample     = X_samp,
    shapviz_obj  = sv_obj,
    plot         = p,
    importance   = importance_df
  ))
}


# Internal: ggplot2 beeswarm fallback
.ggplot_beeswarm <- function(shap_mat, X_mat, model_label, N, B) {
  feat_names <- colnames(shap_mat)
  
  shap_long <- do.call(rbind, lapply(feat_names, function(f) {
    data.frame(
      feature    = f,
      shap_value = shap_mat[, f],
      feat_raw   = if (is.numeric(X_mat[[f]])) X_mat[[f]]
      else as.numeric(as.factor(X_mat[[f]])),
      stringsAsFactors = FALSE
    )
  }))
  shap_long <- shap_long[!is.na(shap_long$shap_value), ]
  
  shap_long$feat_scaled <- ave(
    shap_long$feat_raw, shap_long$feature,
    FUN = function(x) {
      r <- range(x, na.rm = TRUE)
      if (diff(r) == 0) return(rep(0.5, length(x)))
      (x - r[1]) / diff(r)
    }
  )
  
  feat_order <- names(sort(colMeans(abs(shap_mat), na.rm = TRUE),
                           decreasing = FALSE))
  shap_long$feature <- factor(shap_long$feature, levels = feat_order)
  
  ggplot2::ggplot(shap_long,
                  ggplot2::aes(x = shap_value, y = feature, color = feat_scaled)) +
    ggplot2::geom_jitter(height = 0.22, alpha = 0.70, size = 1.6) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5,
                        linetype = "dashed", color = "grey45") +
    ggplot2::scale_color_gradientn(
      colours = c("#3a86ff", "#74b9ff", "#f8f9fa", "#fdcb6e", "#d63031"),
      name    = "Feature value",
      breaks  = c(0, 0.5, 1),
      labels  = c("Low", "Mid", "High")
    ) +
    ggplot2::labs(
      title    = paste0("SHAP Beeswarm -- ", model_label),
      subtitle = sprintf("N = %d samples, B = %d permutations | colour: feature value (red=high, blue=low)", N, B),
      x        = "SHAP value (impact on predicted probability)",
      y        = NULL
    ) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey92",
                                                 linetype = "dotted"),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "right"
    )
}


# -- 4. ExplainBreakDown -------------------------------------------------------

#' Break Down Decomposition (Custom Style, Flexible Filtering)
#'
#' Computes variable attributions for a single prediction and draws a
#' horizontal bar chart with Wes Anderson colors.
#'
#' @param explainer         A DALEX explainer.
#' @param new_observation   Single-row data frame or integer row index.
#' @param type              "break_down_interactions" (default) or "break_down".
#' @param order             Optional vector of variable names.
#' @param remove_zero       Logical; remove variables with zero contribution?
#'                          Default TRUE.
#' @param remove_intercept  Logical; remove the intercept (baseline) bar?
#'                          Default TRUE.
#' @param save_plots        Save plot? Default FALSE.
#' @param save_dir          Output directory.
#' @param plot_width,plot_height  Dimensions.
#'
#' @return A ggplot object. Invisibly returns the DALEX predict_parts object.
#' @export
ExplainBreakDown <- function(explainer,
                             new_observation,
                             type             = "break_down_interactions",
                             order            = NULL,
                             remove_zero      = TRUE,
                             remove_intercept = TRUE,
                             save_plots       = FALSE,
                             save_dir         = "ModelExplain",
                             plot_width       = 8,
                             plot_height      = 5) {
  .check_xai_packages()
  
  obs <- .resolve_observation(explainer, new_observation)
  cat("-- Break Down Decomposition -------------------------------------\n")
  
  args <- list(explainer = explainer, new_observation = obs, type = type)
  if (!is.null(order)) args$order <- order
  bd <- do.call(DALEX::predict_parts, args)
  
  bd_df <- as.data.frame(bd)
  
  # Identify variable and contribution columns
  var_col <- intersect(c("variable_name", "variable"), colnames(bd_df))[1]
  if (is.na(var_col)) stop("Could not find variable name column in break_down object.")
  val_col <- intersect(c("contribution"), colnames(bd_df))[1]
  if (is.na(val_col)) stop("Could not find contribution column in break_down object.")
  
  bd_df <- bd_df[, c(var_col, val_col)]
  colnames(bd_df) <- c("variable", "value")
  
  # Label the intercept explicitly
  intercept_idx <- which(bd_df$variable %in% c("intercept", "baseline", "_intercept_", "_baseline_") |
                           is.na(bd_df$variable) | bd_df$variable == "")
  if (length(intercept_idx) > 0) {
    bd_df$variable[intercept_idx] <- "Intercept"
  }
  
  # Filtering
  if (remove_zero) {
    zero_idx <- which(bd_df$value == 0 & bd_df$variable != "Intercept")
    if (length(zero_idx) > 0) bd_df <- bd_df[-zero_idx, ]
  }
  if (remove_intercept) {
    intercept_idx <- which(bd_df$variable == "Intercept")
    if (length(intercept_idx) > 0) bd_df <- bd_df[-intercept_idx, ]
  }
  
  if (nrow(bd_df) == 0) stop("No variables left after filtering.")
  
  # Order for visual clarity
  bd_df <- bd_df[order(abs(bd_df$value), decreasing = TRUE), ]
  bd_df$variable <- factor(bd_df$variable, levels = rev(bd_df$variable))
  
  cols <- wesanderson::wes_palette("Royal1", 2, type = "discrete")
  bd_df$fill_col <- ifelse(bd_df$value >= 0, cols[1], cols[2])
  
  p <- ggplot2::ggplot(bd_df, ggplot2::aes(x = value, y = variable, fill = fill_col)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.4f", value),
                   hjust = ifelse(value >= 0, -0.2, 1.2)),
      size = 3.5, colour = "black"
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(
      title    = paste0("Break Down -- ", explainer$label),
      subtitle = if (type == "break_down_interactions") "With interaction detection" else "Sequential variable attribution",
      x        = "Contribution to prediction",
      y        = NULL
    ) +
    .pub_theme(13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.text.y   = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  
  if (save_plots) {
    if (is.null(save_dir) || !nzchar(save_dir)) save_dir <- "ModelExplain"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = file.path(save_dir, "breakdown.pdf"), plot = p, width = plot_width, height = plot_height)
    cat("breakdown saved to:", file.path(save_dir, "breakdown.pdf"), "\n")
  }
  
  invisible(bd)
}


# -- 5. ExplainCeterisParibus -------------------------------------------------

#' Ceteris Paribus What-If Curves
#'
#' Shows how the model prediction changes when a single variable varies
#' while all other features are held constant.
#'
#' @param explainer            A DALEX explainer.
#' @param new_observation      Single-row data frame or row index.
#' @param variables            Variables to plot (default: all continuous).
#' @param categorical_variables  Categorical variables (displayed as bar charts).
#' @param save_plots           Save plot? Default FALSE.
#' @param save_dir             Output directory.
#' @param plot_width,plot_height  Dimensions.
#'
#' @return Invisibly returns the predict_profile object.
#' @export
ExplainCeterisParibus <- function(explainer,
                                  new_observation,
                                  variables             = NULL,
                                  categorical_variables = NULL,
                                  save_plots            = FALSE,
                                  save_dir              = "ModelExplain",
                                  plot_width            = 8,
                                  plot_height           = 5) {
  .check_xai_packages()
  
  obs <- .resolve_observation(explainer, new_observation)
  cat("-- Ceteris Paribus Curves ---------------------------------------\n")
  
  cp_args <- list(explainer = explainer, new_observation = obs)
  if (!is.null(variables)) cp_args$variables <- variables
  cp <- do.call(DALEX::predict_profile, cp_args)
  
  cont_vars <- if (!is.null(variables)) {
    setdiff(variables, categorical_variables)
  } else {
    names(Filter(is.numeric, explainer$data))
  }
  
  if (length(cont_vars) > 0) {
    p_cont <- plot(cp, variables = cont_vars) +
      ggplot2::labs(
        title    = paste0("Ceteris Paribus -- ", explainer$label),
        subtitle = "Predicted probability when varying a single feature"
      ) +
      ggprism::theme_prism(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    print(p_cont)
    
    if (save_plots) {
      if (is.null(save_dir) || !nzchar(save_dir)) save_dir <- "ModelExplain"
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      ggsave(filename = file.path(save_dir, "cp_continuous.pdf"),
             plot = p_cont, width = plot_width, height = plot_height)
      cat("Continuous CP plot saved to:", file.path(save_dir, "cp_continuous.pdf"), "\n")
    }
  }

  if (!is.null(categorical_variables) && length(categorical_variables) > 0) {
    available_cat <- intersect(categorical_variables, names(cp))
    if (length(available_cat) == 0) {
      warning("None of the specified categorical variables are available in the profile object. Skipping categorical plot.")
    } else {
      p_cat <- plot(cp, variables = available_cat, categorical_type = "bars") +
        ggplot2::labs(
          title    = paste0("Ceteris Paribus (categorical) -- ", explainer$label),
          subtitle = "Predicted probability per category"
        ) +
        ggprism::theme_prism(base_size = 13) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
      print(p_cat)
      
      if (save_plots) {
        if (is.null(save_dir) || !nzchar(save_dir)) save_dir <- "ModelExplain"
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
        ggsave(filename = file.path(save_dir, "cp_categorical.pdf"),
               plot = p_cat, width = plot_width, height = plot_height)
        cat("Categorical CP plot saved to:", file.path(save_dir, "cp_categorical.pdf"), "\n")
      }
    }
  }
  invisible(cp)
}


# -- 8. ExplainPartialDependence -----------------------------------------------

#' Partial Dependence Plots (PDP / ALE)
#'
#' Averages Ceteris Paribus profiles over the background dataset to obtain
#' global partial dependence. Supports grouping, clustering, and Accumulated
#' Local Effects (ALE).
#'
#' @param explainer   A DALEX explainer.
#' @param variables   Variables to plot (default: first 4 continuous).
#' @param groups      Grouping variable (factor column name).
#' @param k           Number of clusters (integer).
#' @param type        Type: `"partial"` (PDP, default), `"conditional"`, or
#'   `"accumulated"` (ALE).
#' @param N           Number of samples for estimation (default 500).
#' @param save_plots  Save PDF? Default FALSE.
#' @param save_dir    Output directory.
#' @param plot_width,plot_height  Dimensions.
#'
#' @return Invisibly returns the model_profile object.
#' @export
ExplainPartialDependence <- function(explainer,
                                     variables   = NULL,
                                     groups      = NULL,
                                     k           = NULL,
                                     type        = "partial",
                                     N           = 500,
                                     save_plots  = FALSE,
                                     save_dir    = "ModelExplain",
                                     plot_width  = 8,
                                     plot_height = 5) {
  .check_xai_packages()
  cat("-- Partial Dependence (PDP) -------------------------------------\n")
  
  mp_args <- list(explainer = explainer, type = type, N = N)
  if (!is.null(variables)) mp_args$variables <- variables
  if (!is.null(groups))    mp_args$groups    <- groups
  if (!is.null(k))         mp_args$k         <- k
  pdp <- do.call(DALEX::model_profile, mp_args)
  
  plot_vars <- variables %||%
    names(Filter(is.numeric, explainer$data))[
      seq_len(min(4L, sum(sapply(explainer$data, is.numeric))))
    ]
  
  p <- plot(pdp, variables = plot_vars) +
    ggplot2::labs(
      title    = paste0("Partial Dependence -- ", explainer$label),
      subtitle = paste0(
        switch(type,
               partial = "PDP (averaged CP profiles)",
               conditional = "Conditional Marginal Profiles",
               accumulated = "ALE (Accumulated Local Effects)"),
        if (!is.null(groups)) paste0(" | group: ", groups) else "",
        if (!is.null(k))      paste0(" | ", k, " clusters")  else ""
      )
    ) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  print(p)
  
  if (save_plots) {
    if (is.null(save_dir) || length(save_dir) == 0 || !nzchar(save_dir)) {
      save_dir <- "."
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    }
    base_name <- paste0("pdp_", type, if (!is.null(groups)) paste0("_by_", groups) else "")
    file_path <- file.path(save_dir, paste0(base_name, ".pdf"))
    ggplot2::ggsave(filename = file_path, plot = p, width = plot_width, height = plot_height)
    cat("Plot saved to:", file_path, "\n")
  }
  invisible(pdp)
}




# -- Internal helper functions ------------------------------------------------

#' @keywords internal
.resolve_observation <- function(explainer, new_observation) {
  if (is.numeric(new_observation) && length(new_observation) == 1L &&
      new_observation == as.integer(new_observation)) {
    idx <- as.integer(new_observation)
    if (idx < 1L || idx > nrow(explainer$data))
      stop(sprintf("Row index %d is out of range (1:%d).", idx, nrow(explainer$data)))
    return(explainer$data[idx, , drop = FALSE])
  }
  if (!is.data.frame(new_observation))
    stop("new_observation must be a single-row data frame or an integer row index.")
  new_observation
}

#' @keywords internal
.save_plot <- function(
    p,
    save_plots,
    save_dir,
    name,
    width,
    height
) {

    if (!isTRUE(save_plots))
        return(invisible(NULL))

    if (is.null(save_dir) ||
        length(save_dir) != 1 ||
        !is.character(save_dir) ||
        is.na(save_dir) ||
        !nzchar(save_dir)) {

        save_dir <- "ModelExplain"
    }

    if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
    }

    path <- file.path(save_dir, paste0(name, ".pdf"))

    ggplot2::ggsave(
        filename = path,
        plot = p,
        width = width,
        height = height,
        dpi = 300
    )

    cat("Plot saved to:", path, "\n")
    invisible(path)
}
