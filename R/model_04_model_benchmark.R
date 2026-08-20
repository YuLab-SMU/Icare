#' Preprocessing Benchmark for Binary Classification
#'
#' Evaluates the impact of different preprocessing strategies (imputation and
#' normalization) on the performance of multiple classification algorithms.
#' 
#' @param data A data.frame, Train_Model, or Stat object containing the data.
#' @param group_col Character string indicating the column name of the response variable (binary).
#' @param algorithms Character vector of caret model names to evaluate.
#' @param impute_methods Character vector of imputation methods: "none", "median", "knn", "bag".
#' @param norm_methods List of normalization methods passed to caret::preProcess.
#' @param outer_folds Number of outer cross-validation folds. Default: 5.
#' @param inner_cv_number Number of folds for inner repeated CV. Default: 5.
#' @param inner_repeats Number of repetitions for inner CV. Default: 2.
#' @param tuneLength Length of tuning parameter grid. Default: 3.
#' @param seed Random seed for reproducibility. Default: 123.
#' @param verbose Logical; if TRUE, print progress messages. Default: TRUE.
#' @param metrics Character vector of performance metrics to compute.
#' 
#' @return A data frame with class "PreprocessingBenchmark".
#' @export
PreprocessingBenchmark <- function(
    data,
    group_col = NULL,
    algorithms = c('glm',"lda", "knn","nb","rpart","rf","svmRadial"),
    impute_methods = c("none", "median", "knn", "bag"),
    norm_methods = list("none", c("center", "scale"), c("center", "scale", "YeoJohnson")),
    outer_folds = 5,
    inner_cv_number = 5,
    inner_repeats = 2,
    tuneLength = 3,
    seed = 123,
    verbose = TRUE,
    metrics = c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1", "Kappa")
) {
  
  # 1. Handle S4 objects and data extraction
  if (methods::is(data, "Train_Model")) {
    df <- data@clean.df
    if (is.null(group_col)) group_col <- data@group_col
  } else if (methods::is(data, "Stat")) {
    df <- data@clean.data
    if (is.null(group_col)) group_col <- data@group_col
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("data must be a data.frame, Train_Model, or Stat object.")
  }
  
  if (is.null(group_col)) stop("group_col must be provided.")
  if (!group_col %in% names(df)) stop("group_col not found in data.")
  
  # 2. Dependency checks
  required_pkgs <- c("caret", "pROC", "dplyr")
  missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))
  
  if ("knn" %in% impute_methods && !requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for 'knn' imputation.")
  }
  if ("bag" %in% impute_methods && !requireNamespace("ipred", quietly = TRUE)) {
    stop("Package 'ipred' is required for 'bag' imputation.")
  }
  
  # 3. Safe binary outcome handling
  df[[group_col]] <- as.factor(df[[group_col]])
  if (length(levels(df[[group_col]])) != 2) stop("Only binary classification supported.")
  levels(df[[group_col]]) <- make.names(levels(df[[group_col]]))
  positive_class <- levels(df[[group_col]])[2]
  
  # 4. Validate metrics
  valid_metrics <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1", "Kappa")
  metrics <- intersect(metrics, valid_metrics)
  if (length(metrics) == 0) stop("No valid metrics selected.")
  
  # 5. Inner trainControl
  inner_control <- caret::trainControl(
    method = "repeatedcv",
    number = inner_cv_number,
    repeats = inner_repeats,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  # 6. Outer folds
  set.seed(seed)
  folds <- caret::createFolds(df[[group_col]], k = outer_folds, list = TRUE)
  
  results_list <- list()
  error_log <- list()
  
  method_map <- c(median = "medianImpute", knn = "knnImpute", bag = "bagImpute")
  feature_cols <- setdiff(names(df), group_col)
  
  # 7. Main Loops
  for (algo in algorithms) {
    if (verbose) cat("Processing", algo, "...\n")
    
    for (imp in impute_methods) {
      for (norm in norm_methods) {
        norm_name <- if (length(norm) == 1 && norm == "none") "none" else paste(norm, collapse = "_")
        
        fold_metrics_list <- list()
        
        for (k in seq_along(folds)) {
          test_idx <- folds[[k]]
          train_df <- df[-test_idx, , drop = FALSE]
          test_df  <- df[ test_idx, , drop = FALSE]
          
          train_x <- train_df[, feature_cols, drop = FALSE]
          test_x  <- test_df[, feature_cols, drop = FALSE]
          train_y <- train_df[[group_col]]
          test_y  <- test_df[[group_col]]
          
          # Construct preProcess steps
          pp_steps <- c()
          if (imp != "none") pp_steps <- c(pp_steps, method_map[imp])
          if (!(length(norm) == 1 && norm == "none")) pp_steps <- c(pp_steps, norm)
          
          # Apply preprocessing
          if (length(pp_steps) > 0) {
            pp_obj <- tryCatch({
              caret::preProcess(train_x, method = pp_steps)
            }, error = function(e) {
              error_log[[length(error_log) + 1]] <<- data.frame(
                Algorithm = algo, Imputation = imp, Normalization = norm_name,
                Fold = k, Error = paste("preProcess:", e$message), stringsAsFactors = FALSE
              )
              NULL
            })
            if (is.null(pp_obj)) next
            
            train_x_proc <- predict(pp_obj, train_x)
            test_x_proc  <- predict(pp_obj, test_x)
          } else {
            train_x_proc <- train_x
            test_x_proc  <- test_x
          }
          
          train_final <- data.frame(train_x_proc, Class = train_y)
          test_final  <- data.frame(test_x_proc,  Class = test_y)
          
          # Train model
          train_args <- list(
            form = Class ~ ., data = train_final, method = algo,
            metric = "ROC", maximize = TRUE, trControl = inner_control, tuneLength = tuneLength
          )
          if (algo %in% c("glm", "glmStepAIC")) train_args$family <- "binomial"
          
          fit <- tryCatch({
            suppressWarnings(suppressMessages(do.call(caret::train, train_args)))
          }, error = function(e) {
            error_log[[length(error_log) + 1]] <<- data.frame(
              Algorithm = algo, Imputation = imp, Normalization = norm_name,
              Fold = k, Error = as.character(e$message), stringsAsFactors = FALSE
            )
            NULL
          })
          if (is.null(fit)) next
          
          # Predict probabilities & classes
          eval_res <- tryCatch({
            pred_prob <- predict(fit, newdata = test_final, type = "prob")
            prob_val  <- if (positive_class %in% colnames(pred_prob)) pred_prob[, positive_class] else pred_prob[, 2]
            
            pred_class <- predict(fit, newdata = test_final, type = "raw")
            
            list(prob = prob_val, class = pred_class)
          }, error = function(e) {
            error_log[[length(error_log) + 1]] <<- data.frame(
              Algorithm = algo, Imputation = imp, Normalization = norm_name,
              Fold = k, Error = paste("Prediction:", e$message), stringsAsFactors = FALSE
            )
            NULL
          })
          if (is.null(eval_res)) next
          
          # Compute AUC
          auc_val <- tryCatch({
            roc_obj <- pROC::roc(response = test_y, predictor = eval_res$prob, levels = levels(test_y), quiet = TRUE)
            as.numeric(roc_obj$auc)
          }, error = function(e) NA)
          
          # Compute Confusion Matrix and metrics via caret
          cm_obj <- tryCatch({
            caret::confusionMatrix(eval_res$class, test_y, positive = positive_class)
          }, error = function(e) NULL)
          
          if (is.null(cm_obj)) next
          all_metrics <- c(
            AUC         = auc_val,
            Accuracy    = as.numeric(cm_obj$overall["Accuracy"]),
            Sensitivity = as.numeric(cm_obj$byClass["Sensitivity"]),
            Specificity = as.numeric(cm_obj$byClass["Specificity"]),
            F1          = as.numeric(cm_obj$byClass["F1"]),
            Kappa       = as.numeric(cm_obj$overall["Kappa"])
          )
          
          fold_metrics_list[[k]] <- all_metrics[metrics]
        }
        
        # Summarize results
        if (length(fold_metrics_list) > 0) {
          fold_df <- do.call(rbind, fold_metrics_list)
          fold_df <- fold_df[apply(!is.na(fold_df), 1, any), , drop = FALSE]
          n_success <- nrow(fold_df)
          
          if (n_success > 0) {
            mean_vals <- colMeans(fold_df, na.rm = TRUE)
            sd_vals   <- apply(fold_df, 2, sd, na.rm = TRUE)
            
            result_row <- data.frame(
              Algorithm     = algo,
              Imputation    = imp,
              Normalization = norm_name,
              n_success     = n_success,
              failure_rate  = (outer_folds - n_success) / outer_folds,
              stringsAsFactors = FALSE
            )
            
            for (m in metrics) {
              m_mean <- mean_vals[m]
              m_sd   <- sd_vals[m]
              se     <- m_sd / sqrt(n_success)
              ci_lower <- m_mean - stats::qt(0.975, df = n_success - 1) * se
              ci_upper <- m_mean + stats::qt(0.975, df = n_success - 1) * se
              
              result_row[[paste0("Mean_", m)]]     <- m_mean
              result_row[[paste0("SD_", m)]]       <- m_sd
              result_row[[paste0("CI_lower_", m)]] <- ci_lower
              result_row[[paste0("CI_upper_", m)]] <- ci_upper
            }
            results_list[[length(results_list) + 1]] <- result_row
          }
        }
      }
    }
  }
  
  results_df <- dplyr::bind_rows(results_list)
  attr(results_df, "error_log") <- dplyr::bind_rows(error_log)
  class(results_df) <- c("PreprocessingBenchmark", class(results_df))
  
  if (verbose) cat("Benchmark completed.\n")
  return(results_df)
}

#' Plot Forest Plot of Benchmark Results
#'
#' Creates a faceted forest plot to compare the performance of different
#' preprocessing strategies across algorithms. The plot shows point estimates
#' with confidence intervals for a chosen metric.
#'
#' @param benchmark_result A data frame returned by \code{\link{PreprocessingBenchmark}}.
#' @param metric Character string specifying which metric to plot.
#'   Must match one of the columns in the benchmark result, e.g.,
#'   \code{"AUC"}, \code{"Accuracy"}, \code{"Sensitivity"}, \code{"Specificity"},
#'   \code{"F1"}, or \code{"Kappa"}. Default: \code{"AUC"}.
#' @param palette_name Character string naming the Wes Anderson palette to use
#'   for algorithm colors. Default: \code{"Zissou1"}.
#' @param remove_na Logical; if \code{TRUE}, removes rows with missing values
#'   for the selected metric. Default: \code{TRUE}.
#' @param global_median Logical; if \code{TRUE}, adds a dashed vertical line
#'   indicating the global median of the selected metric across all facets.
#'   Default: \code{FALSE}.
#' @param save_plot Logical; if \code{TRUE}, saves the plot to a PDF file.
#'   Default: \code{FALSE}.
#' @param save_dir Character string; directory where the plot will be saved.
#'   Required if \code{save_plot = TRUE}. Default: \code{NULL}.
#'
#' @return Invisibly returns the \code{ggplot} object. The plot is drawn on
#'   the current graphics device.
#'
#' @examples
#' \dontrun{
#'   res <- PreprocessingBenchmark(my_data, group_col = "Class")
#'   PlotBenchmarkForest(res, metric = "Accuracy", global_median = TRUE)
#' }
#' @export
PlotBenchmarkForest <- function(benchmark_result,
                                metric = "AUC",
                                palette_name = "Zissou1",
                                remove_na = TRUE,
                                global_median = FALSE,
                                save_plot = FALSE,
                                save_dir = NULL) {
  
  # Check required packages
  required_pkgs <- c("ggplot2", "dplyr", "wesanderson", "ggprism")
  missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages for plotting: ", paste(missing, collapse = ", "))
  }
  
  # Validate input
  if (!inherits(benchmark_result, "PreprocessingBenchmark") && !is.data.frame(benchmark_result)) {
    stop("Input must be the output from PreprocessingBenchmark()")
  }
  
  df <- as.data.frame(benchmark_result)
  
  # Construct column names for the selected metric
  mean_col <- paste0("Mean_", metric)
  lower_col <- paste0("CI_lower_", metric)
  upper_col <- paste0("CI_upper_", metric)
  
  # Check if the metric exists
  if (!mean_col %in% colnames(df)) {
    stop("Metric '", metric, "' not found in benchmark results. Available: ",
         paste(grep("^Mean_", colnames(df), value = TRUE), collapse = ", "))
  }
  
  # Remove rows with NA for the selected metric if requested
  if (remove_na) {
    df <- df[!is.na(df[[mean_col]]), ]
    if (nrow(df) == 0) stop("No valid results after removing NAs for metric '", metric, "'.")
  }
  
  # Ensure all required columns exist; if CI columns missing, set to mean
  if (!lower_col %in% colnames(df)) {
    df[[lower_col]] <- df[[mean_col]]
    df[[upper_col]] <- df[[mean_col]]
    ci_label <- "No CI"
  } else {
    ci_label <- "95% CI"
  }
  
  # Factor imputation methods in a sensible order
  df$Imputation <- factor(df$Imputation,
                          levels = c("none", "median", "knn", "bag"))
  
  # Order algorithms by median of the selected metric (descending)
  algo_order <- df %>%
    dplyr::group_by(Algorithm) %>%
    dplyr::summarise(med = median(.data[[mean_col]], na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(med)) %>%
    dplyr::pull(Algorithm)
  df$Algorithm <- factor(df$Algorithm, levels = algo_order)
  
  # Color palette for algorithms
  n_algos <- length(unique(df$Algorithm))
  my_cols <- wesanderson::wes_palette(palette_name, n = n_algos, type = "continuous")
  
  # Facet medians (median per facet)
  facet_stats <- df %>%
    dplyr::group_by(Imputation, Normalization) %>%
    dplyr::summarise(Facet_Median = median(.data[[mean_col]], na.rm = TRUE), .groups = "drop")
  
  # Global median (if requested)
  global_med <- median(df[[mean_col]], na.rm = TRUE)
  
  # Build base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[mean_col]], y = Algorithm,
                                        color = Algorithm)) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = .data[[lower_col]],
                                          xmax = .data[[upper_col]]),
                             linewidth = 0.85, alpha = 0.9) +
    ggplot2::facet_grid(Imputation ~ Normalization, scales = "free_y") +
    ggplot2::scale_color_manual(values = my_cols) +
    ggplot2::geom_text(data = facet_stats,
                       ggplot2::aes(x = -Inf, y = Inf,
                                    label = paste0("Median: ", round(Facet_Median, 3))),
                       vjust = 1.6, hjust = -0.05, size = 3.3,
                       color = "black", fontface = "bold", inherit.aes = FALSE) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::labs(x = paste("Mean", metric, "+/-", ci_label),
                  y = NULL,
                  title = paste("Preprocessing Benchmark -", metric, "Performance")) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )+
    ggplot2::geom_vline(xintercept = global_med, linetype = "dashed",
                        color = "gray40", alpha = 0.7)
  
  # Add global median line if requested
  if (global_median) {
    p <- p +ggplot2::annotate("text", x = global_med, y = -Inf,
                        label = paste("Global Median:", round(global_med, 3)),
                        vjust = -0.8, hjust = 0.5, size = 3.8, color = "gray30")
  }
  
  # Print plot
  print(p)
  
  # Save if requested
  if (save_plot) {
    if (is.null(save_dir)) stop("save_dir must be specified when save_plot = TRUE")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    filename <- file.path(save_dir, paste0("Preprocessing_Benchmark_Forest_", metric, ".pdf"))
    ggplot2::ggsave(filename, plot = p, width = 14, height = 10, dpi = 300)
    cat("Plot saved:", filename, "\n")
  }
  
  invisible(p)
}


#' Logistic Regression Diagnostic Benchmark (Univariate + Combined)
#'
#' @description
#' Establishes a baseline diagnostic performance benchmark for a binary
#' classification problem using \strong{logistic regression only}: one
#' univariate \code{glm(family = "binomial")} model per feature, plus one
#' combined multivariate model using all features together. Both are trained
#' via \code{caret::train()} with a user-selectable resampling scheme
#' (k-fold CV or LOOCV), and ROC curves are built from pooled
#' out-of-fold (cross-validated) predictions rather than from a single
#' train/test split, so univariate and combined curves are on the same
#' footing and directly comparable.
#'
#' This is meant as a reference/baseline step before comparing more complex
#' models (rf, gbm, xgbTree, ...): "how much signal is in each single
#' biomarker on its own, and how much do we gain by combining all of them
#' with plain logistic regression?"
#'
#' @param object A \code{Train_Model} or \code{Stat} S4 object (or any object
#'   supported by the internal \code{.extract_xy()} helper) containing a
#'   binary outcome and a set of numeric/factor predictors.
#' @param cv_method Resampling scheme: \code{"cv"} (k-fold, default) or
#'   \code{"LOOCV"} (leave-one-out). Repeated CV is intentionally not
#'   supported here, since pooling out-of-fold predictions across repeats
#'   would give each observation multiple, non-independent ROC contributions.
#' @param number Number of folds when \code{cv_method = "cv"} (e.g. 5 or 10).
#'   Ignored for LOOCV. Default \code{10}.
#' @param seed Random seed for reproducibility. Default \code{825}.
#' @param top_label_n How many top univariate features (by AUC) to label
#'   directly on the plot. Default \code{5}. Set to \code{0} to disable
#'   labels.
#' @param verbose Logical; print per-feature progress. Default \code{TRUE}.
#' @param save_plot Logical; save the ROC plot to PDF. Default \code{FALSE}.
#' @param save_dir Directory to save the plot to when \code{save_plot = TRUE}.
#'   Defaults to \code{"./figures/"}.
#' @param width,height PDF dimensions in inches. Defaults \code{8, 8}.
#' @param preProcess Character vector of preprocessing methods to apply
#'   to the predictors before model fitting (e.g., \code{c("center", "scale")},
#'   \code{"BoxCox"}). Passed directly to \code{caret::train()}. Default
#'   \code{NULL} (no preprocessing).
#' @param smote Logical. If \code{TRUE}, Synthetic Minority Over-sampling
#'   Technique (SMOTE) is applied within each cross-validation fold to
#'   handle class imbalance. Default \code{FALSE}.
#'
#' @return An object of S3 class \code{"LogisticBenchmark"} (a list) with:
#' \describe{
#'   \item{performance_table}{data.frame with columns \code{Feature}, \code{Type}
#'     (\code{"Univariate"} or \code{"Combined"}), \code{AUC}, \code{AUC_low},
#'     \code{AUC_high} (DeLong 95\% CI), \code{n_features}, sorted by AUC
#'     descending.}
#'   \item{roc_curves}{named list of \code{pROC::roc} objects, one per
#'     feature plus \code{"Combined (all features)"}.}
#'   \item{models}{named list of the underlying \code{caret::train} fit
#'     objects.}
#'   \item{plot}{the \code{ggplot} ROC object (also drawn on the current
#'     device).}
#'   \item{errors}{named list of features that failed to train, with reasons.}
#'   \item{call}{the matched call.}
#' }
#'
#' @details
#' \strong{Why out-of-fold predictions instead of \code{fit$results$ROC}:}
#' \code{fit$results$ROC} is the mean of per-fold AUCs, which is the right
#' number for comparing models but cannot be drawn as a single ROC curve.
#' Pooling \code{fit$pred} (each observation's prediction from the one fold
#' where it was held out) gives one prediction per subject, from which a
#' single, plottable ROC curve can be computed — this mirrors what
#' \code{caret}'s own \code{twoClassSummary} does internally to get the ROC
#' number, just applied to the pooled set instead of averaged per-fold.
#'
#' \strong{Positive class convention:} following \code{caret::twoClassSummary},
#' the first level of the outcome factor (\code{levels(y)[1]}) is treated as
#' the class whose probability column drives the ROC curve, for consistency
#' with any other caret-based comparison in the same pipeline.
#'
#' @examples
#' \dontrun{
#' bench <- LogisticDiagnosticBenchmark(
#'   object = train_obj,
#'   cv_method = "cv",
#'   number = 10,
#'   top_label_n = 5,
#'   preProcess = c("center", "scale"),
#'   smote = TRUE
#' )
#' bench$performance_table
#' bench$plot
#' }
#'
#' @importFrom caret train trainControl twoClassSummary
#' @importFrom pROC roc auc ci.auc
#' @export
LogisticDiagnosticBenchmark <- function(object,
                                        cv_method = c("cv", "LOOCV"),
                                        number = 10,
                                        seed = 825,
                                        top_label_n = 5,
                                        verbose = TRUE,
                                        save_plot = FALSE,
                                        save_dir = NULL,
                                        width = 8,
                                        height = 8,
                                        preProcess = NULL,
                                        smote = FALSE) {
  
  cv_method <- match.arg(cv_method)
  .check_fs_packages()
  if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  set.seed(seed)
  if (smote && !requireNamespace("themis", quietly = TRUE)) {
    stop("Package 'themis' is required when smote = TRUE. Please install it.")
  }
  # ---------------------------------------------------------------------
  # 1. Extract & clean data
  # ---------------------------------------------------------------------
  xy <- .extract_xy(object)
  if (!is.factor(xy$y)) xy$y <- as.factor(xy$y)
  xy$y <- droplevels(xy$y)
  if (nlevels(xy$y) != 2) {
    stop("LogisticDiagnosticBenchmark() requires exactly 2 classes. Found ",
         nlevels(xy$y), ".")
  }
  levels(xy$y) <- make.names(levels(xy$y))
  df <- cbind(xy$x, group = xy$y)
  
  if (anyNA(df)) {
    n_bad <- sum(!stats::complete.cases(df))
    warning("Removing ", n_bad, " rows with missing values.")
    df <- df[stats::complete.cases(df), ]
  }
  if (nrow(df) == 0) stop("No complete rows left after NA removal.")
  
  char_cols <- names(df)[vapply(df, is.character, logical(1))]
  for (col in char_cols) df[[col]] <- as.factor(df[[col]])
  
  full_features <- setdiff(colnames(df), "group")
  if (length(full_features) == 0) stop("No predictor columns available.")
  
  lev <- levels(df$group)  # lev[1] = positive class, matching caret::twoClassSummary convention
  
  # ---------------------------------------------------------------------
  # 2. Resampling control (with optional SMOTE)
  # ---------------------------------------------------------------------
  n_samples <- nrow(df)
  if (cv_method == "cv" && n_samples < number) {
    warning("n_samples (", n_samples, ") < number (", number, "). Switching to LOOCV.")
    cv_method <- "LOOCV"
  }
  
  # Build trainControl, adding sampling if SMOTE is requested
  ctrl <- caret::trainControl(
    method          = cv_method,
    number          = if (cv_method == "cv") number else NULL,
    classProbs      = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    verboseIter     = FALSE,
    allowParallel   = FALSE,
    sampling        = if (smote) "smote" else NULL   # SMOTE handling
  )
  
  # ---------------------------------------------------------------------
  # 3. Helper: fit one glm (uni- or multi-variate) and pool CV predictions
  # ---------------------------------------------------------------------
  .fit_and_pool <- function(formula_rhs, label) {
    form <- stats::as.formula(paste("group ~", formula_rhs))
    saw_missing_metric <- FALSE
    
    fit <- tryCatch(
      withCallingHandlers(
        caret::train(form, data = df, method = "glm", family = "binomial",
                     trControl = ctrl,
                     metric = "ROC",
                     preProcess = preProcess),   # 应用预处理
        warning = function(w) {
          if (grepl("all the .* metric values are missing", conditionMessage(w))) {
            saw_missing_metric <<- TRUE
          }
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit) || saw_missing_metric) {
      return(list(
        fit = NULL, roc = NULL,
        reason = if (saw_missing_metric) {
          "All CV folds produced NA performance (likely separation for this feature)."
        } else {
          "Training error."
        }
      ))
    }
    
    pred_df <- fit$pred
    if (is.null(pred_df) || nrow(pred_df) == 0 || !lev[1] %in% colnames(pred_df)) {
      return(list(fit = fit, roc = NULL, reason = "No usable out-of-fold predictions."))
    }
    
    roc_obj <- tryCatch(pROC::roc(pred_df$obs, pred_df[[lev[1]]], 
                                  direction = "auto", levels = lev, quiet = TRUE), error = function(e) NULL)
    if (is.null(roc_obj)) {
      return(list(fit = fit, roc = NULL, reason = "ROC computation failed on pooled predictions."))
    }
    
    list(fit = fit, roc = roc_obj, reason = NULL)
  }
  
  # ---------------------------------------------------------------------
  # 4. Univariate models
  # ---------------------------------------------------------------------
  models     <- list()
  roc_curves <- list()
  error_log  <- list()
  perf_rows  <- list()
  
  for (f in full_features) {
    if (verbose) cat(sprintf("Univariate glm: %-40s ", f))
    out <- .fit_and_pool(sprintf("`%s`", f), f)
    if (is.null(out$roc)) {
      if (verbose) cat("FAILED (", out$reason, ")\n")
      error_log[[f]] <- out$reason
      next
    }
    auc_val <- as.numeric(pROC::auc(out$roc))
    ci_val  <- tryCatch(as.numeric(pROC::ci.auc(out$roc, quiet = TRUE)), error = function(e) c(NA, NA, NA))
    if (verbose) cat(sprintf("AUC = %.4f\n", auc_val))
    
    models[[f]]     <- out$fit
    roc_curves[[f]] <- out$roc
    perf_rows[[f]]  <- data.frame(
      Feature = f, Type = "Univariate", AUC = auc_val,
      AUC_low = ci_val[1], AUC_high = ci_val[3], n_features = 1,
      stringsAsFactors = FALSE
    )
  }
  
  # ---------------------------------------------------------------------
  # 5. Combined (all-features) model
  # ---------------------------------------------------------------------
  if (verbose) cat(sprintf("Combined glm (%d features): ", length(full_features)))
  combined_label <- "Combined (all features)"
  out <- .fit_and_pool(".", combined_label)
  if (is.null(out$roc)) {
    if (verbose) cat("FAILED (", out$reason, ")\n")
    error_log[[combined_label]] <- out$reason
  } else {
    auc_val <- as.numeric(pROC::auc(out$roc))
    ci_val  <- tryCatch(as.numeric(pROC::ci.auc(out$roc, quiet = TRUE)), error = function(e) c(NA, NA, NA))
    if (verbose) cat(sprintf("AUC = %.4f\n", auc_val))
    
    models[[combined_label]]     <- out$fit
    roc_curves[[combined_label]] <- out$roc
    perf_rows[[combined_label]]  <- data.frame(
      Feature = combined_label, Type = "Combined", AUC = auc_val,
      AUC_low = ci_val[1], AUC_high = ci_val[3], n_features = length(full_features),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(perf_rows) == 0) {
    stop("No logistic regression models could be trained. See errors for details.")
  }
  
  performance_table <- do.call(rbind, perf_rows)
  rownames(performance_table) <- NULL
  performance_table <- performance_table[order(-performance_table$AUC), ]
  
  # ---------------------------------------------------------------------
  # 6. Plot: combined = red, univariate = grayscale ordered by AUC
  # ---------------------------------------------------------------------
  uni_order <- performance_table$Feature[performance_table$Type == "Univariate"]
  n_uni <- length(uni_order)
  
  roc_df_list <- lapply(names(roc_curves), function(nm) {
    r <- roc_curves[[nm]]
    data.frame(
      Sensitivity = r$sensitivities,
      Specificity = r$specificities,
      Model = nm,
      stringsAsFactors = FALSE
    )
  })
  roc_df <- do.call(rbind, roc_df_list)
  
  # Grayscale palette: darker = higher AUC, ordered by performance_table.
  # Kept within grey30 (dark) to grey80 (light) so lines stay visible
  # without competing with the red combined-model line.
  gray_vals <- if (n_uni > 0) grDevices::grey(seq(0.25, 0.80, length.out = n_uni)) else character(0)
  names(gray_vals) <- uni_order
  
  color_map <- c(gray_vals, setNames("red", combined_label))
  # Only keep entries that actually got a curve (in case some failed)
  color_map <- color_map[names(color_map) %in% names(roc_curves)]
  
  roc_df$Model <- factor(roc_df$Model, levels = names(color_map))
  size_map  <- setNames(rep(0.6, length(color_map)), names(color_map))
  if (combined_label %in% names(size_map)) size_map[combined_label] <- 1.6
  alpha_map <- setNames(rep(0.85, length(color_map)), names(color_map))
  
  p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - Specificity, y = Sensitivity,
                                            color = Model, linewidth = Model, alpha = Model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::scale_linewidth_manual(values = size_map, guide = "none") +
    ggplot2::scale_alpha_manual(values = alpha_map, guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Logistic Regression Diagnostic Benchmark",
      subtitle = sprintf("%s (%s%s) | Red = combined model, greyscale = univariate features",
                         if (cv_method == "LOOCV") "LOOCV" else sprintf("%d-fold CV", number),
                         cv_method, if (cv_method == "cv") paste0(", n=", number) else ""),
      x = "1 - Specificity", y = "Sensitivity"
    ) +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                   legend.text = ggplot2::element_text(size = 7),
                   plot.subtitle = ggplot2::element_text(size = 9))
  
  if (top_label_n > 0 && n_uni > 0) {
    top_feats <- utils::head(uni_order, top_label_n)
    label_pts <- do.call(rbind, lapply(top_feats, function(nm) {
      r <- roc_curves[[nm]]
      best_idx <- which.max(r$sensitivities + r$specificities)
      data.frame(x = 1 - r$specificities[best_idx], y = r$sensitivities[best_idx],
                 label = sprintf("%s (%.3f)", nm, as.numeric(pROC::auc(r))),
                 stringsAsFactors = FALSE)
    }))
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = label_pts,
                                        ggplot2::aes(x = x, y = y, label = .data$label),
                                        inherit.aes = FALSE, size = 2.8, color = "grey30",
                                        max.overlaps = 20)
    }
  }
  
  print(p)
  
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- "./figures/"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "logistic_diagnostic_benchmark.pdf"),
                    plot = p, width = width, height = height, dpi = 300)
  }
  
  structure(
    list(
      performance_table = performance_table,
      roc_curves        = roc_curves,
      models            = models,
      plot              = p,
      errors            = error_log,
      call              = match.call()
    ),
    class = "LogisticBenchmark"
  )
}
