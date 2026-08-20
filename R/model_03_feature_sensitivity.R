# =============================================================================
# Feature Sensitivity Analysis - Backward Elimination by Importance
# =============================================================================
# This module performs backward feature elimination using caret models.
# For each algorithm it records CV AUC (mean +/- SE) while removing the least
# important feature at each step. Elbow points are identified by one of four
# methods:
#   - "difference"         : maximum absolute drop in AUC
#   - "ratio"              : maximum relative drop in AUC
#   - "perf_tolerance"     : min features with AUC >= max_AUC * (1 - tol)
#   - "stability_window"   : first stable window (range <= tolerance) of length window_size
# =============================================================================

# 0. Package check -----------------------------------------------------------
.check_fs_packages <- function() {
  required <- c("caret", "ggplot2", "wesanderson", "ggprism", "dplyr", "tidyr",
                "reshape2", "ggrepel", "doParallel", "foreach", "parallel",
                "RColorBrewer")
  missing  <- required[!sapply(required, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })]
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nPlease install them.", call. = FALSE)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Internal helpers (same as before)
# ------------------------------------------------------------------------------
.extract_xy <- function(object) {
  if (inherits(object, "Train_Model")) {
    cd <- object@clean.df
    gc <- object@group_col
  } else if (inherits(object, "Stat")) {
    cd <- object@clean.data
    gc <- object@group_col
  } else {
    stop("Object must be Train_Model or Stat.", call. = FALSE)
  }
  x <- cd[, setdiff(colnames(cd), gc), drop = FALSE]
  colnames(x) <- make.names(colnames(x), unique = TRUE)
  y <- cd[[gc]]
  list(x = x, y = y, group_col = make.names(gc, unique = TRUE))
}

.validate_model <- function(algo) {
  model_info <- tryCatch(caret::getModelInfo(algo, regex = FALSE), error = function(e) NULL)
  if (is.null(model_info) || length(model_info) == 0) {
    return(list(valid = FALSE, message = sprintf("Model '%s' not found", algo)))
  }
  model_info <- model_info[[1]]
  needed_pkgs <- unique(model_info$library)
  missing_pkg <- needed_pkgs[!sapply(needed_pkgs, function(pkg) requireNamespace(pkg, quietly = TRUE))]
  if (length(missing_pkg) > 0) {
    return(list(valid = FALSE, message = sprintf("Missing packages: %s", paste(missing_pkg, collapse=", "))))
  }
  if (!"Classification" %in% model_info$type) {
    return(list(valid = FALSE, message = "Not a classification model"))
  }
  return(list(valid = TRUE, message = "OK"))
}

.run_algorithm_sensitivity <- function(df, full_features, algo, metric, number, min_features, verbose) {
  current_features <- full_features
  step_results <- data.frame()
  step_features <- list()
  iteration <- 0
  max_iterations <- length(full_features) + 10
  while (length(current_features) >= 1 && iteration < max_iterations) {
    iteration <- iteration + 1
    n_feat <- length(current_features)
    if (verbose) cat(sprintf("  Training with %d features... ", n_feat))
    step_features[[as.character(n_feat)]] <- current_features
    keep_cols <- intersect(c(current_features, "group"), colnames(df))
    if (!"group" %in% keep_cols) stop("Response column 'group' not found")
    sub_df <- df[, keep_cols, drop = FALSE]
    complete_cases <- complete.cases(sub_df)
    if (sum(complete_cases) < number) {
      warning(sprintf("Insufficient complete cases (%d) for %d-fold CV", sum(complete_cases), number), call. = FALSE)
      break
    }
    fit_result <- .train_single_model(sub_df = sub_df, algo = algo, metric = metric, number = number)
    if (is.null(fit_result$fit)) { if (verbose) cat("FAILED\n"); break }
    if (verbose) cat(sprintf("AUC = %.4f\n", fit_result$auc_mean))
    step_results <- rbind(step_results, data.frame(n_features = n_feat,
                                                   AUC_mean = fit_result$auc_mean,
                                                   AUC_se = fit_result$auc_se,
                                                   stringsAsFactors = FALSE))
    if (n_feat == 1) break
    removal_result <- .remove_least_important(fit = fit_result$fit, current_features = current_features, verbose = verbose)
    if (is.null(removal_result)) break
    current_features <- removal_result$remaining_features
  }
  if (iteration >= max_iterations) warning("Maximum iterations reached")
  if (nrow(step_results) > 0) step_results <- step_results[order(-step_results$n_features), ]
  return(list(step_results = step_results, step_features = step_features))
}
.train_single_model <- function(sub_df, algo, metric, number) {
  ctrl <- caret::trainControl(method = "cv", number = number, classProbs = TRUE,
                              summaryFunction = caret::twoClassSummary,
                              savePredictions = TRUE, verboseIter = FALSE)
  
  fit <- tryCatch({
    if (algo == "glm"|algo == "glmStepAIC") {
      caret::train(group ~ ., data = sub_df, method = algo, trControl = ctrl,
                   metric = metric, tuneLength = 1, 
                   family = "binomial")
    } else {
      caret::train(group ~ ., data = sub_df, method = algo, trControl = ctrl,
                   metric = metric, tuneLength = 1)
    }
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(list(fit = NULL, auc_mean = NA, auc_se = NA))
  perf <- .extract_performance(fit, metric, number)
  return(list(fit = fit, auc_mean = perf$auc_mean, auc_se = perf$auc_se))
}

.extract_performance <- function(fit, metric, number) {
  perf_cols <- colnames(fit$results)
  metric_col <- grep(paste0("^", metric), perf_cols, value = TRUE)[1]
  if (is.na(metric_col)) {
    warning(sprintf("Metric '%s' not found", metric), call. = FALSE)
    return(list(auc_mean = NA, auc_se = NA))
  }
  sd_col <- paste0(metric_col, "SD")
  if (nrow(fit$results) == 1) {
    best_row <- fit$results[1, ]
  } else {
    best_tune <- fit$bestTune
    if (!is.null(best_tune) && nrow(best_tune) == 1) {
      idx <- which(apply(fit$results[, names(best_tune), drop = FALSE], 1, function(r) all(r == best_tune)))
      if (length(idx) == 0) idx <- which.max(fit$results[[metric_col]])
    } else {
      idx <- which.max(fit$results[[metric_col]])
    }
    best_row <- fit$results[idx[1], ]
  }
  auc_mean <- best_row[[metric_col]]
  if (!is.null(sd_col) && sd_col %in% colnames(fit$results)) {
    auc_se <- best_row[[sd_col]] / sqrt(number)
  } else if (!is.null(fit$resample) && metric_col %in% colnames(fit$resample)) {
    auc_se <- sd(fit$resample[[metric_col]], na.rm = TRUE) / sqrt(number)
  } else {
    auc_se <- 0
  }
  return(list(auc_mean = auc_mean, auc_se = auc_se))
}

.remove_least_important <- function(fit, current_features, verbose) {
  imp <- tryCatch(caret::varImp(fit, scale = TRUE)$importance, error = function(e) NULL)
  if (is.null(imp) || nrow(imp) == 0) {
    if (verbose) cat("  Unable to compute feature importance\n")
    return(NULL)
  }
  imp_vec <- imp[, 1]
  names(imp_vec) <- make.names(rownames(imp), unique = TRUE)
  common_feats <- intersect(names(imp_vec), current_features)
  if (length(common_feats) == 0) {
    if (verbose) cat("  No matching features in importance scores\n")
    return(NULL)
  }
  imp_vec <- imp_vec[common_feats]
  sorted_imp <- sort(imp_vec, decreasing = TRUE)
  feat_to_remove <- names(sorted_imp)[length(sorted_imp)]
  if (verbose) cat(sprintf("  Removing feature: %s (importance: %.4f)\n", feat_to_remove, sorted_imp[length(sorted_imp)]))
  remaining_features <- setdiff(current_features, feat_to_remove)
  return(list(removed_feature = feat_to_remove, remaining_features = remaining_features))
}

.apply_smoothing <- function(step_results, smooth_span) {
  n <- nrow(step_results)
  if (smooth_span > 0 && n >= 5) {
    step_results$AUC_smooth <- tryCatch(
      predict(loess(AUC_mean ~ n_features, data = step_results, span = smooth_span), step_results$n_features),
      error = function(e) { warning("Smoothing failed; using raw AUC", call. = FALSE); step_results$AUC_mean }
    )
  } else {
    step_results$AUC_smooth <- step_results$AUC_mean
  }
  return(step_results)
}

# ------------------------------------------------------------------------------
# Elbow detection methods
# ------------------------------------------------------------------------------
.elbow_difference <- function(auc_vals, n_features, n) {
  changes <- c(NA, diff(auc_vals))
  max_drop_idx <- which.min(changes)
  elbow_idx <- if (max_drop_idx > 1) max_drop_idx - 1 else 1
  return(n_features[elbow_idx])
}

.elbow_ratio <- function(auc_vals, n_features, n) {
  change_rates <- c(NA, diff(auc_vals) / auc_vals[-n])
  max_drop_idx <- which.min(change_rates)
  elbow_idx <- if (max_drop_idx > 1) max_drop_idx - 1 else 1
  return(n_features[elbow_idx])
}

.elbow_perf_tolerance <- function(auc_vals, n_features, tol) {
  max_auc <- max(auc_vals)
  target <- max_auc * (1 - tol)
  valid <- which(auc_vals >= target)
  if (length(valid) == 0) return(n_features[1])
  optimal_idx <- valid[length(valid)]
  return(n_features[optimal_idx])
}

.elbow_stability_window <- function(auc_vals, n_features, window_size, tolerance) {
  n <- length(auc_vals)
  if (n < window_size) {
    warning("Number of points less than window_size; returning maximum features")
    return(n_features[1])
  }
  stability_scores <- sapply(1:(n - window_size + 1), function(i) {
    window <- auc_vals[i:(i + window_size - 1)]
    max(window) - min(window)
  })
  stable_idx <- which(stability_scores <= tolerance)[1]
  if (is.na(stable_idx)) {
    min_score_idx <- which.min(stability_scores)
    stable_idx <- min_score_idx
    warning("No window met stability tolerance; using window with minimum AUC range")
  }
  optimal_idx <- stable_idx + window_size - 1
  return(n_features[optimal_idx])
}

# ------------------------------------------------------------------------------
# Main exported functions
# ------------------------------------------------------------------------------
#' Run Backward Feature Elimination for Multiple Models
#'
#' Performs backward feature elimination on a set of caret models using
#' cross-validation. At each step, the least important feature is removed
#' based on variable importance from the current model. The function returns
#' a structured object containing performance curves and selected feature sets
#' for each algorithm.
#'
#' @param object A \code{Train_Model} or \code{Stat} S4 object containing the
#'   training data with features and a binary outcome column.
#' @param models Character vector of caret model names to evaluate.
#'   Default is \code{c("rf", "gbm", "glmboost")}. Any classification model
#'   that supports variable importance can be used.
#' @param metric Evaluation metric for model performance. Default is \code{"ROC"}.
#'   Other options include \code{"Accuracy"}, \code{"Kappa"}, etc., as supported
#'   by \code{caret::train}.
#' @param number Number of cross-validation folds. Default is \code{5}.
#' @param seed Random seed for reproducibility. Default is \code{825}.
#' @param min_features Minimum number of features to retain during elimination.
#'   The process stops when this many features remain. Default is \code{2}.
#' @param smooth_span Smoothing span for loess curve applied to the performance
#'   trajectory (0 = no smoothing). Default is \code{0}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'   Default is \code{TRUE}.
#' @param parallel Logical. If \code{TRUE}, parallel processing is enabled
#'   using \code{doParallel}. Default is \code{FALSE}.
#' @param n_cores Integer. Number of cores for parallel processing. If
#'   \code{NULL} and \code{parallel = TRUE}, uses \code{detectCores() - 1}.
#'   Default is \code{NULL}.
#'
#' @return An object of S3 class \code{"FeatureElimination"} containing:
#'   \itemize{
#'     \item \code{results}: A data frame with columns \code{n_features},
#'       \code{AUC_mean}, \code{AUC_se}, \code{AUC_smooth}, and \code{Algorithm}.
#'     \item \code{step_features}: A list of feature sets retained at each step,
#'       organized by algorithm.
#'     \item \code{models}: Character vector of successful algorithms.
#'     \item \code{metric}, \code{number}, \code{min_features}, \code{smooth_span}:
#'       Parameters used.
#'     \item \code{call}: The matched call.
#'   }
#'
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item For each algorithm, the model is trained on all features.
#'   \item Variable importance is extracted and the least important feature
#'     is removed.
#'   \item The model is retrained on the reduced feature set, and performance
#'     is recorded.
#'   \item Steps 2–3 are repeated until \code{min_features} is reached.
#' }
#' Models that fail to train or provide importance at any step are skipped.
#' The performance trajectory can be smoothed using \code{smooth_span}.
#'
#' @examples
#' \dontrun{
#' # Assuming 'model_obj' is a Train_Model object
#' elim <- run_feature_elimination(
#'   object = model_obj,
#'   models = c("rf", "gbm"),
#'   metric = "ROC",
#'   number = 5,
#'   min_features = 3,
#'   verbose = TRUE
#' )
#'
#' # Access results
#' head(elim$results)
#' # Select optimal features
#' best <- select_elbow(elim, elbow_method = "difference")
#' print(best$best_features)
#' }
#'
#' @importFrom caret train trainControl varImp
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach registerDoSEQ
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   elim <- run_feature_elimination(model, models = c("rf"), number = 3, min_features = 2)
#'   best <- select_elbow(elim, elbow_method = "difference")
#'   print(best$optimal_counts)
#' }
#' }
run_feature_elimination <- function(object,
                                    models = c("rf", "gbm", "glmboost"),
                                    metric = "ROC",
                                    number = 5,
                                    seed = 825,
                                    min_features = 2,
                                    smooth_span = 0,
                                    verbose = TRUE,
                                    parallel = FALSE,
                                    n_cores = NULL) {
  .check_fs_packages()
  set.seed(seed)
  if (!inherits(object, c("Train_Model", "Stat"))) stop("'object' must be Train_Model or Stat")
  if (!is.character(models) || length(models) == 0) stop("'models' must be non-empty character vector")
  
  # Exclude models without built-in varImp in the first step
  support_df <- check_varImp_availability(models)
  unsupported <- support_df$Model[!support_df$Has_BuiltIn]
  if (length(unsupported) > 0) {
    if (verbose) {
      cat("Note: The following models do NOT have built-in varImp and will be skipped:\n")
      cat(paste("  -", unsupported), sep = "\n")
    }
    models <- setdiff(models, unsupported)
  }
  if (length(models) == 0) {
    stop("No valid models remaining after filtering out unsupported models via check_varImp_availability().")
  }
  
  xy <- .extract_xy(object)
  full_features <- colnames(xy$x)
  if (length(full_features) < min_features) stop(sprintf("Only %d features, but min_features=%d", length(full_features), min_features))
  levels(xy$y) <- make.names(levels(xy$y))
  df <- cbind(xy$x, group = xy$y)
  df$group <- factor(df$group, levels = unique(df$group))
  if (length(levels(df$group)) != 2) stop("Response must be binary for metric ROC")
  
  all_results <- list()
  all_step_features <- list()
  
  for (algo in models) {
    if (verbose) cat(sprintf("\n====== Algorithm: %s ======\n", algo))
    model_check <- .validate_model(algo)
    if (!model_check$valid) {
      warning(sprintf("Skipping '%s': %s", algo, model_check$message), call. = FALSE)
      next
    }
    algo_out <- .run_algorithm_sensitivity(
      df = df, full_features = full_features, algo = algo,
      metric = metric, number = number, min_features = min_features,
      verbose = verbose
    )
    if (is.null(algo_out) || nrow(algo_out$step_results) == 0) {
      if (verbose) cat(sprintf("  No valid results for %s\n", algo))
      next
    }
    algo_out$step_results <- .apply_smoothing(algo_out$step_results, smooth_span)
    algo_out$step_results$Algorithm <- algo
    all_results[[algo]] <- algo_out$step_results
    all_step_features[[algo]] <- algo_out$step_features
  }
  if (length(all_results) == 0) stop("No algorithms succeeded")
  combined_results <- do.call(rbind, all_results)
  rownames(combined_results) <- NULL
  structure(
    list(
      results = combined_results,
      step_features = all_step_features,
      models = names(all_results),
      metric = metric,
      number = number,
      min_features = min_features,
      smooth_span = smooth_span,
      call = match.call()
    ),
    class = "FeatureElimination"
  )
}

#' Select Optimal Feature Count from Elimination Results
#'
#' Determines the optimal number of features to retain for each algorithm
#' based on the performance curve from backward elimination. Four elbow
#' detection methods are available, each identifying the point where
#' performance begins to plateau or drop.
#'
#' @param elim_obj A \code{FeatureElimination} object returned by
#'   \code{\link{run_feature_elimination}}.
#' @param elbow_method Character string specifying the elbow detection method.
#'   Options:
#'   \itemize{
#'     \item \code{"difference"}: Maximum absolute drop in AUC between
#'       consecutive steps.
#'     \item \code{"ratio"}: Maximum relative drop in AUC.
#'     \item \code{"perf_tolerance"}: Minimum features with AUC within
#'       \code{tol} of the maximum AUC.
#'     \item \code{"stability_window"}: First stable window of length
#'       \code{window_size} with AUC range <= \code{stability_tol}.
#'   }
#'   Default is \code{"difference"}.
#' @param tol Numeric. Tolerance for \code{"perf_tolerance"} method. The
#'   threshold is defined as \code{max_AUC * (1 - tol)}. Default is \code{0.05}.
#' @param window_size Integer. Window size for \code{"stability_window"}
#'   method. Default is \code{5}.
#' @param stability_tol Numeric. Maximum allowed AUC range within a stable
#'   window for \code{"stability_window"} method. Default is \code{0.01}.
#' @param min_features Integer. Minimum number of features to retain.
#'   If \code{NULL}, uses the value from \code{elim_obj$min_features}.
#'   Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages and
#'   adjustments. Default is \code{TRUE}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{best_features}: A named list where each element is a
#'       character vector of feature names selected for the corresponding
#'       algorithm.
#'     \item \code{optimal_counts}: A named integer vector of the optimal
#'       number of features for each algorithm.
#'   }
#'
#' @details
#' The four elbow detection methods work as follows:
#' \describe{
#'   \item{\code{"difference"}}{Calculates the absolute change in AUC between
#'     consecutive steps and selects the feature count before the largest drop.}
#'   \item{\code{"ratio"}}{Similar to \code{"difference"}, but uses the
#'     relative change (drop divided by the current AUC).}
#'   \item{\code{"perf_tolerance"}}{Finds the smallest feature set whose AUC
#'     is within \code{tol * 100\%} of the maximum AUC achieved.}
#'   \item{\code{"stability_window"}}{Identifies the first window of
#'     \code{window_size} consecutive steps where the AUC range is
#'     \code{<= stability_tol}, indicating that performance has stabilised.}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'elim_obj' is a FeatureElimination object
#'
#' # Use the default difference method
#' result <- select_elbow(elim_obj)
#' print(result$optimal_counts)
#'
#' # Use performance tolerance method
#' result_tol <- select_elbow(elim_obj,
#'   elbow_method = "perf_tolerance",
#'   tol = 0.02
#' )
#'
#' # Use stability window method with custom parameters
#' result_stab <- select_elbow(elim_obj,
#'   elbow_method = "stability_window",
#'   window_size = 4,
#'   stability_tol = 0.005
#' )
#'
#' # Access selected features for each algorithm
#' str(result$best_features)
#' }
select_elbow <- function(elim_obj,
                         elbow_method = c("difference", "ratio", "perf_tolerance", "stability_window"),
                         tol = 0.05,
                         window_size = 5,
                         stability_tol = 0.01,
                         min_features = NULL,
                         verbose = TRUE) {
  if (!inherits(elim_obj, "FeatureElimination")) stop("'elim_obj' must be a FeatureElimination object")
  elbow_method <- match.arg(elbow_method)
  if (is.null(min_features)) min_features <- elim_obj$min_features
  if (elbow_method == "perf_tolerance" && (tol < 0 || tol > 1)) {
    stop("'tol' must be between 0 and 1 for perf_tolerance method")
  }
  if (elbow_method == "stability_window" && (window_size < 2)) {
    stop("'window_size' must be at least 2 for stability_window method")
  }
  
  best_sets <- list()
  algo_names <- unique(elim_obj$results$Algorithm)
  
  for (algo in algo_names) {
    sub <- elim_obj$results[elim_obj$results$Algorithm == algo, ]
    sub <- sub[order(sub$n_features, decreasing = TRUE), ]
    n <- nrow(sub)
    auc_vals <- sub$AUC_smooth
    n_feats <- sub$n_features
    
    optimal_n <- switch(
      elbow_method,
      "difference" = .elbow_difference(auc_vals, n_feats, n),
      "ratio" = .elbow_ratio(auc_vals, n_feats, n),
      "perf_tolerance" = .elbow_perf_tolerance(auc_vals, n_feats, tol),
      "stability_window" = .elbow_stability_window(auc_vals, n_feats, window_size, stability_tol)
    )
    
    if (optimal_n < min_features) {
      if (verbose) cat(sprintf("  Optimal n=%d adjusted to min_features=%d\n", optimal_n, min_features))
      optimal_n <- min_features
    }
    if (!optimal_n %in% n_feats) {
      warning(sprintf("Optimal n=%d not found; using closest", optimal_n), call. = FALSE)
      optimal_n <- n_feats[which.min(abs(n_feats - optimal_n))]
    }
    best_sets[[algo]] <- elim_obj$step_features[[algo]][[as.character(optimal_n)]]
  }
  list(best_features = best_sets,
       optimal_counts = sapply(best_sets, length))
}

#' Get selected features for multiple elbow methods at once
#' @param elim_obj A \code{FeatureElimination} object.
#' @param methods Character vector of elbow methods to apply.
#' @param ... Additional arguments passed to \code{select_elbow} (e.g., tol, window_size, stability_tol).
#' @return A list where each element is the result of \code{select_elbow} for a given method.
#' @export
get_selected_features <- function(elim_obj, methods = c("difference", "ratio", "perf_tolerance", "stability_window"), ...) {
  if (!inherits(elim_obj, "FeatureElimination")) stop("'elim_obj' must be a FeatureElimination object")
  results <- list()
  for (m in methods) {
    results[[m]] <- select_elbow(elim_obj, elbow_method = m, ...)
  }
  return(results)
}

#' Legacy Wrapper: Feature Elimination + Elbow Selection in One Call
#'
#' This is a convenience wrapper that combines \code{\link{run_feature_elimination}}
#' and \code{\link{select_elbow}} into a single function call. It performs
#' backward feature elimination on the specified models, then selects the optimal
#' feature count using the chosen elbow method, returning a structured result
#' object ready for downstream analysis and plotting.
#'
#' @param object A \code{Train_Model} or \code{Stat} object containing the data.
#' @param models Character vector of caret model names to evaluate.
#'   Default: \code{c("rf", "gbm", "glmboost")}.
#' @param metric Evaluation metric for model performance. Default: \code{"ROC"}.
#' @param number Number of cross-validation folds. Default: \code{5}.
#' @param seed Random seed for reproducibility. Default: \code{825}.
#' @param elbow_method Elbow detection method. One of \code{"difference"},
#'   \code{"ratio"}, \code{"perf_tolerance"}, or \code{"stability_window"}.
#'   Default: \code{"difference"}.
#' @param tol Tolerance for \code{"perf_tolerance"} method. Default: \code{0.05}.
#' @param window_size Window size for \code{"stability_window"} method. Default: \code{5}.
#' @param stability_tol Stability tolerance for \code{"stability_window"} method.
#'   Default: \code{0.01}.
#' @param min_features Minimum number of features to retain. Default: \code{2}.
#' @param smooth_span Smoothing span for loess curve (0 = no smoothing).
#'   Default: \code{0}.
#' @param verbose Logical. Print progress messages. Default: \code{TRUE}.
#' @param parallel Logical. Use parallel processing? Default: \code{FALSE}.
#' @param n_cores Number of cores for parallel processing. If \code{NULL},
#'   uses \code{parallel::detectCores() - 1}. Default: \code{NULL}.
#'
#' @return An object of class \code{FeatureSensitivityAnalysis} (inherits from list)
#'   with two components:
#'   \itemize{
#'     \item \code{results}: The full elimination results data frame.
#'     \item \code{best_features}: A named list of selected features per model.
#'   }
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   sens <- FeatureSensitivityAnalysis(model, models = c("rf"), number = 3)
#'   print(sens$best_features)
#' }
#' }
FeatureSensitivityAnalysis <- function(object,
                                       models = c("rf", "gbm", "glmboost"),
                                       metric = "ROC",
                                       number = 5,
                                       seed = 825,
                                       elbow_method = c("difference", "ratio", "perf_tolerance", "stability_window"),
                                       tol = 0.05,
                                       window_size = 5,
                                       stability_tol = 0.01,
                                       min_features = 2,
                                       smooth_span = 0,
                                       verbose = TRUE,
                                       parallel = FALSE,
                                       n_cores = NULL) {
  elim <- run_feature_elimination(object = object, models = models, metric = metric,
                                  number = number, seed = seed, min_features = min_features,
                                  smooth_span = smooth_span, verbose = verbose,
                                  parallel = parallel, n_cores = n_cores)
  elbow <- select_elbow(elim, elbow_method = elbow_method, tol = tol,
                        window_size = window_size, stability_tol = stability_tol,
                        min_features = min_features, verbose = verbose)
  result <- list(results = elim$results, best_features = elbow$best_features)
  class(result) <- c("FeatureSensitivityAnalysis", "list")
  result
}

# ------------------------------------------------------------------------------
# Plotting with low-saturation palette
# ------------------------------------------------------------------------------
.get_low_sat_palette <- function(n) {
  if (n <= 8) {
    pal <- RColorBrewer::brewer.pal(min(n, 8), "Set3")
    if (n < 8) pal <- pal[1:n]
    return(pal)
  } else {
    return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(n))
  }
}

#' Plot Elbow Curve from Feature Elimination Results
#'
#' Visualizes the performance curve (AUC) as a function of the number of
#' features retained during backward elimination. The plot can display
#' confidence intervals (error bars or ribbons), highlight optimal feature
#' counts, and optionally facet by algorithm.
#'
#' @param elim_obj A \code{FeatureElimination} object returned by
#'   \code{\link{run_feature_elimination}}.
#' @param best_features A named list of selected features per algorithm,
#'   typically from \code{\link{select_elbow}}. If provided, optimal points
#'   are highlighted on the plot. Default is \code{NULL}.
#' @param facet Logical. If \code{TRUE}, separate panels are created for each
#'   algorithm. Default is \code{FALSE}.
#' @param palette_name Character string specifying the Wes Anderson palette
#'   name for algorithm colors. Default is \code{"Darjeeling1"}.
#' @param use_low_sat Logical. If \code{TRUE}, uses a low-saturation color
#'   palette for better print readability. Default is \code{TRUE}.
#' @param ci_style Character string specifying the style for confidence
#'   intervals. Options are \code{"errorbar"}, \code{"ribbon"}, or
#'   \code{"none"}. Default is \code{"errorbar"}.
#' @param save_plot Logical. If \code{TRUE}, the plot is saved to a PDF file.
#'   Default is \code{FALSE}.
#' @param save_dir Character string specifying the directory to save the plot.
#'   If \code{NULL} and \code{save_plot = TRUE}, the plot is saved to
#'   \code{"./Sensitivity/"}. Default is \code{NULL}.
#' @param width Numeric. Width of the saved plot in inches. If \code{facet = TRUE},
#'   defaults to \code{12}, otherwise \code{8}.
#' @param height Numeric. Height of the saved plot in inches. Default is \code{5}.
#' @param ncol Integer. Number of columns for faceted plots. Ignored if
#'   \code{facet = FALSE}. Default is \code{NULL}.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}, such as
#'   \code{subtitle} or \code{caption}.
#'
#' @return Invisibly returns the \code{ggplot} object. The plot is drawn on
#'   the current graphics device and optionally saved to PDF.
#' @export
#' @examples
#' \dontrun{
#' # Basic plot
#' plot_elbow(elim_obj)
#'
#' # With optimal features highlighted and faceted
#' best <- select_elbow(elim_obj, elbow_method = "difference")
#' plot_elbow(elim_obj, best_features = best$best_features, facet = TRUE)
#'
#' # Custom confidence interval style and save
#' plot_elbow(elim_obj,
#'   ci_style = "ribbon",
#'   save_plot = TRUE,
#'   save_dir = "./figures",
#'   width = 10,
#'   height = 6
#' )
#' }
plot_elbow <- function(elim_obj,
                       best_features = NULL,
                       facet = FALSE,
                       palette_name = "Darjeeling1",
                       use_low_sat = TRUE,
                       ci_style = c("errorbar", "ribbon", "none"),
                       save_plot = FALSE,
                       save_dir = NULL,
                       width = if(facet) 12 else 8,
                       height = 5,
                       ncol = NULL,
                       ...) {
  if (!inherits(elim_obj, "FeatureElimination")) {
    stop("'elim_obj' must be a FeatureElimination object from run_feature_elimination()")
  }
  ci_style <- match.arg(ci_style)
  df <- elim_obj$results
  if (nrow(df) == 0) stop("No results to plot")
  
  opts <- NULL
  if (!is.null(best_features) && length(best_features) > 0) {
    best_n <- sapply(best_features, length)
    algo_names <- names(best_n)
    opts <- data.frame(Algorithm = algo_names, n_features = best_n, stringsAsFactors = FALSE)
    if (all(c("Algorithm", "n_features") %in% colnames(df))) {
      opts <- merge(opts, df, by = c("Algorithm", "n_features"), all.x = TRUE, sort = FALSE)
    }
  }
  
  n_algos <- length(unique(df$Algorithm))
  if (use_low_sat) {
    cols <- .get_low_sat_palette(n_algos)
  } else {
    cols <- wesanderson::wes_palette(palette_name, n_algos, type = "discrete")
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n_features, y = AUC_mean, color = Algorithm))
  if (ci_style == "ribbon") {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = AUC_mean - AUC_se, ymax = AUC_mean + AUC_se, fill = Algorithm),
                                  alpha = 0.15, colour = NA) +
      ggplot2::scale_fill_manual(values = cols, guide = "none")
  } else if (ci_style == "errorbar") {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = AUC_mean - AUC_se, ymax = AUC_mean + AUC_se),
                                    width = 0.3, alpha = 0.6)
  }
  p <- p + ggplot2::geom_line(linewidth = 1.2) + ggplot2::geom_point(size = 2)
  
  if (!is.null(opts) && nrow(opts) > 0) {
    p <- p + ggplot2::geom_point(data = opts, size = 4, shape = 8, color = "red") +
      ggrepel::geom_text_repel(data = opts, ggplot2::aes(label = paste0("n=", n_features)),
                               color = "red", size = 4, fontface = "bold",
                               nudge_x = -0.5, nudge_y = 0.02)
  }
  
  if (facet) {
    p <- p + ggplot2::facet_wrap(~ Algorithm, scales = "free_y", ncol = ncol) +
      ggplot2::scale_color_manual(values = cols, guide = "none")
  } else {
    p <- p + ggplot2::scale_color_manual(values = cols)
  }
  
  p <- p + ggplot2::labs(title = "Feature Sensitivity Analysis (Backward Elimination)",
                         x = "Number of Features", y = "AUC (mean +/- SE)", ...) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom")
  
  print(p)
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- "./Sensitivity/"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    filename <- if(facet) "sensitivity_elbow_faceted.pdf" else "sensitivity_elbow.pdf"
    ggplot2::ggsave(file.path(save_dir, filename), plot = p, width = width, height = height, dpi = 300)
  }
  invisible(p)
}

# Deprecated plotting functions (to avoid confusion)
PlotSensitivityElbow <- function(sens_result, ...) {
  stop("PlotSensitivityElbow is deprecated. Use plot_elbow(elim_obj, best_features = sens_result$best_features, ...)")
}
PlotSensitivityElbowFaceted <- function(sens_result, ...) {
  stop("PlotSensitivityElbowFaceted is deprecated. Use plot_elbow(elim_obj, best_features = sens_result$best_features, facet = TRUE, ...)")
}

# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------
#' @export
print.FeatureElimination <- function(x, ...) {
  cat("Feature Elimination Object\n")
  cat("==========================\n")
  cat(sprintf("Algorithms: %s\n", paste(unique(x$results$Algorithm), collapse=", ")))
  cat(sprintf("CV folds: %d\n", x$number))
  cat(sprintf("Min features: %d\n", x$min_features))
  cat(sprintf("Smoothing span: %.2f\n", x$smooth_span))
  cat("\nUse summary() for detailed results.\n")
  invisible(x)
}

#' @export
print.FeatureSensitivityAnalysis <- function(x, ...) {
  cat("Feature Sensitivity Analysis Results\n")
  cat("=====================================\n\n")
  n_algos <- length(unique(x$results$Algorithm))
  cat(sprintf("Algorithms analyzed: %d\n", n_algos))
  if (length(x$best_features) > 0) {
    cat("\nOptimal feature counts:\n")
    for (algo in names(x$best_features)) {
      n_feat <- length(x$best_features[[algo]])
      cat(sprintf("  %s: %d features\n", algo, n_feat))
    }
  }
  cat("\nUse summary() for detailed results\n")
  invisible(x)
}

#' @export
summary.FeatureElimination <- function(object, ...) {
  cat("Feature Elimination Summary\n")
  cat("===========================\n")
  print(object$results)
  invisible(object)
}

#' @export
summary.FeatureSensitivityAnalysis <- function(object, ...) {
  cat("Feature Sensitivity Analysis Summary\n")
  cat("=====================================\n\n")
  cat("Results:\n")
  print(object$results)
  cat("\nOptimal Feature Sets:\n")
  cat("---------------------\n")
  for (algo in names(object$best_features)) {
    cat(sprintf("\n%s (%d features):\n", algo, length(object$best_features[[algo]])))
    cat(paste0("  ", object$best_features[[algo]]), sep = "\n")
  }
  invisible(object)
}

# Helper
`%||%` <- function(a, b) if (!is.null(a)) a else b
