# =============================================================================
# model_feature_selection.R
# Comprehensive Feature Selection Module for Train_Model objects
# Integrates RFE, GA, SA, and SBF methods from caret with automated 
# cross-validation and visualization.
# =============================================================================

# 0. Package check ------------------------------------------------------------
.check_fs_packages <- function() {
  required <- c("caret", "ggplot2", "wesanderson", "ggprism", "dplyr", "tidyr",
                "reshape2", "gridExtra", "doParallel", "foreach")
  missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Install them before using feature selection.")
  }
  invisible(TRUE)
}

# -- Internal helpers ---------------------------------------------------------
.get_output_dir <- function(...) {
  get_output_dir(...)
}

.safe_dir <- function(save_dir) {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  invisible(save_dir)
}

.extract_xy <- function(object) {
  # Extract predictors (x) and outcome (y) from a Train_Model object
  if (inherits(object, "Train_Model")) {
    cd <- object@clean.df
    gc <- object@group_col
  } else if (inherits(object, "Stat")) {
    cd <- object@clean.data
    gc <- object@group_col
  } else {
    stop("Object must be Train_Model or Stat.")
  }
  x <- cd[, setdiff(colnames(cd), gc), drop = FALSE]
  y <- cd[[gc]]
  list(x = x, y = y, group_col = gc)
}

# -- 1. RFE: Recursive Feature Elimination ------------------------------------

#' Recursive Feature Elimination (RFE) for Train_Model objects
#'
#' @param object   A Train_Model or Stat S4 object.
#' @param sizes    Numeric vector of subset sizes to evaluate.
#' @param rfe_func A caret RFE function list (e.g., \code{rfFuncs}).
#' @param method   External resampling: "cv", "repeatedcv", "boot".
#' @param number   Number of folds / resampling iterations.
#' @param repeats  For repeatedcv only.
#' @param metric   Evaluation metric: "Accuracy", "Kappa", "ROC".
#' @param allowParallel Logical.
#' @param save_plot Logical. Save RFE profile plot?
#' @param save_dir  Directory to save outputs.
#' @param seed      Random seed.
#' @return A list with elements: \code{result} (rfe object), 
#'   \code{opt_vars} (optimal features), \code{plot} (ggplot).
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   rfe_res <- FeatureSelectRFE(model, sizes = c(2,4,6), method = "cv", number = 2)
#'   print(rfe_res$opt_vars)
#' }
#' }
FeatureSelectRFE <- function(object,
                             sizes    = NULL,
                             rfe_func = caret::rfFuncs,
                             method   = "repeatedcv",
                             number   = 5,
                             repeats  = 2,
                             metric   = "Accuracy",
                             allowParallel = FALSE,
                             save_plot = FALSE,
                             save_dir  = NULL,
                             seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  if (is.null(sizes))
    sizes <- unique(round(seq(2, min(30, ncol(xy$x)), length.out = 10)))
  
  ctrl <- caret::rfeControl(
    functions = rfe_func,
    method    = method,
    number    = number,
    repeats   = repeats,
    verbose   = TRUE,
    allowParallel = allowParallel
  )
  
  cat(sprintf("Running RFE: %d sizes, %s (%d-fold, %d repeats)...\n",
              length(sizes), method, number, repeats))
  
  rfe_res <- caret::rfe(xy$x, xy$y, sizes = sizes, rfeControl = ctrl, metric = metric)
  
  opt_vars <- caret::predictors(rfe_res)
  cat(sprintf("RFE selected: %d features\n", length(opt_vars)))
  
  # Build plot
  p <- ggplot2::ggplot(rfe_res, metric = metric) + 
    ggplot2::theme_bw() +
    ggplot2::labs(title = "RFE - Recursive Feature Elimination",
                  subtitle = paste0("Optimal: ", length(opt_vars), " variables"))
  print(p)
  
  if (save_plot) {
    .safe_dir(save_dir)
    ggplot2::ggsave(file.path(save_dir, "rfe_profile.pdf"), p,
                    width = 7, height = 5, device = "pdf")
  }
  
  invisible(list(result = rfe_res, opt_vars = opt_vars, plot = p))
}


# -- 2. GA: Genetic Algorithm Feature Selection -------------------------------

#' Genetic Algorithm Feature Selection for Train_Model objects
#'
#' @param object   A Train_Model or Stat object.
#' @param iters    Number of GA generations.
#' @param popSize  Population size.
#' @param ga_func  A caret GA function list (e.g., \code{caretGA}, \code{rfGA}).
#' @param method   External resampling method.
#' @param number   Number of folds.
#' @param repeats  Repeats.
#' @param metric   Internal fitness metric.
#' @param allowParallel,genParallel Logical.
#' @param save_plot Logical.
#' @param save_dir  Directory.
#' @param seed      Random seed.
#' @return List with \code{result} (gafs object), \code{opt_vars}, \code{plot}.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   ga_res <- FeatureSelectGA(model, iters = 3, popSize = 5, method = "cv", number = 2)
#'   print(ga_res$opt_vars)
#' }
#' }
FeatureSelectGA <- function(object,
                            iters    = 10,
                            popSize  = 20,
                            ga_func  = caret::caretGA,
                            method   = "repeatedcv",
                            number   = 5,
                            repeats  = 2,
                            metric   = "Accuracy",
                            allowParallel = FALSE,
                            genParallel   = FALSE,
                            save_plot = FALSE,
                            save_dir  = NULL,
                            seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  ctrl <- caret::gafsControl(
    functions      = ga_func,
    method         = method,
    number         = number,
    repeats        = repeats,
    verbose        = TRUE,
    allowParallel  = allowParallel,
    genParallel    = genParallel
  )
  
  cat(sprintf("Running GA: %d iters, popSize %d, %s...\n",
              iters, popSize, method))
  
  ga_res <- caret::gafs(xy$x, xy$y, iters = iters, popSize = popSize,
                        gafsControl = ctrl, metric = metric)
  
  opt_vars <- ga_res$optVariables
  cat(sprintf("GA selected: %d features\n", length(opt_vars)))
  
  p <- plot(ga_res) + ggplot2::theme_bw() +
    ggplot2::labs(title = "GA - Genetic Algorithm Feature Selection")
  print(p)
  
  if (save_plot) {
    .safe_dir(save_dir)
    ggplot2::ggsave(file.path(save_dir, "ga_trace.pdf"), p,
                    width = 7, height = 5, device = "pdf")
  }
  
  invisible(list(result = ga_res, opt_vars = opt_vars, plot = p))
}


# -- 3. SA: Simulated Annealing Feature Selection -----------------------------

#' Simulated Annealing Feature Selection for Train_Model objects
#'
#' @param object   A Train_Model or Stat object.
#' @param iters    Number of SA iterations.
#' @param sa_func  A caret SA function list (e.g., \code{caretSA}, \code{rfSA}).
#' @param method   External resampling method.
#' @param number   Number of folds.
#' @param repeats  Repeats.
#' @param metric   Internal fitness metric.
#' @param improve  SA restart after `improve` iters without improvement.
#' @param allowParallel Logical.
#' @param save_plot Logical.
#' @param save_dir  Directory.
#' @param seed      Random seed.
#' @return List with \code{result} (safs object), \code{opt_vars}, \code{plot}.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   sa_res <- FeatureSelectSA(model, iters = 3, method = "cv", number = 2)
#'   print(sa_res$opt_vars)
#' }
#' }
FeatureSelectSA <- function(object,
                            iters    = 25,
                            sa_func  = caret::caretSA,
                            method   = "repeatedcv",
                            number   = 5,
                            repeats  = 2,
                            metric   = "Accuracy",
                            improve  = 3L,
                            allowParallel = FALSE,
                            save_plot = FALSE,
                            save_dir  = NULL,
                            seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  ctrl <- caret::safsControl(
    functions      = sa_func,
    method         = method,
    number         = number,
    repeats        = repeats,
    improve        = improve,
    verbose        = TRUE,
    allowParallel  = allowParallel
  )
  
  cat(sprintf("Running SA: %d iters, %s...\n", iters, method))
  
  sa_res <- caret::safs(xy$x, xy$y, iters = iters,
                        safsControl = ctrl, metric = metric)
  
  opt_vars <- sa_res$optVariables
  cat(sprintf("SA selected: %d features\n", length(opt_vars)))
  
  p <- plot(sa_res) + ggplot2::theme_bw() +
    ggplot2::labs(title = "SA - Simulated Annealing Feature Selection")
  print(p)
  
  if (save_plot) {
    .safe_dir(save_dir)
    ggplot2::ggsave(file.path(save_dir, "sa_trace.pdf"), p,
                    width = 7, height = 5, device = "pdf")
  }
  
  invisible(list(result = sa_res, opt_vars = opt_vars, plot = p))
}


# -- 4. SBF: Selection By Univariate Filter -----------------------------------

#' Univariate Filter Feature Selection for Train_Model objects
#'
#' @param object   A Train_Model or Stat object.
#' @param sbf_func A caret SBF function list (e.g., \code{rfSBF}, \code{caretSBF}).
#' @param method   External resampling method.
#' @param number   Number of folds.
#' @param repeats  Repeats.
#' @param metric   Evaluation metric.
#' @param allowParallel Logical.
#' @param save_plot Logical (saves variable count barplot).
#' @param save_dir  Directory.
#' @param seed      Random seed.
#' @return List with \code{result} (sbf object), \code{opt_vars}, \code{plot}.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   str(mtcars)
#'   mtcars$vs=as.numeric(mtcars$vs)
#'   mtcars$cyl=as.numeric(mtcars$cyl)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   sbf_res <- FeatureSelectSBF(model, method = "cv", number = 5)
#'   print(sbf_res$opt_vars)
#' }
#' }
FeatureSelectSBF <- function(object,
                             sbf_func = caret::rfSBF,
                             method   = "repeatedcv",
                             number   = 5,
                             repeats  = 2,
                             metric   = "Accuracy",
                             allowParallel = FALSE,
                             save_plot = FALSE,
                             save_dir  = NULL,
                             seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  ctrl <- caret::sbfControl(
    functions = sbf_func,
    method    = method,
    number    = number,
    repeats   = repeats,
    verbose   = TRUE
  )
  
  cat(sprintf("Running SBF: %s...\n", method))
  
  sbf_res <- caret::sbf(xy$x, xy$y, sbfControl = ctrl, metric = metric)
  
  opt_vars <- sbf_res$optVariables
  cat(sprintf("SBF selected: %d features\n", length(opt_vars)))
  
  p <- NULL
  tryCatch({
    p <- plot(sbf_res) + ggplot2::theme_bw() +
      ggplot2::labs(title = "SBF - Univariate Filter Feature Selection")
    print(p)
  }, error = function(e) {
    message("SBF plot could not be generated: ", e$message)
    p <- NULL
  })
  
  if (save_plot) {
    .safe_dir(save_dir)
    if (!is.null(p)) {
      ggplot2::ggsave(file.path(save_dir, "sbf_filter.pdf"), p,
                      width = 7, height = 5, device = "pdf")
    } else {
      message("SBF plot not saved due to plotting error.")
    }
  }
  
  invisible(list(result = sbf_res, opt_vars = opt_vars, plot = p))
}


# -- 5. Built-In Importance (from any caret train object) ---------------------
#' Explain Model Performance (ROC, Lift, Boxplot)
#'
#' @param explainer DALEX explainer object.
#' @param geom Plot type: "roc", "lift", or "boxplot".
#' @param save_plots Logical, save plot to file.
#' @param save_dir Directory to save the plot.
#' @param plot_width,plot_height Dimensions in inches.
#' @param ... Additional arguments passed to `ggplot2::ggsave`.
#' @return Invisibly, the DALEX model_performance object.
#' @export
ExplainModelPerformance <- function(explainer,
                                    geom        = c("roc", "lift", "boxplot"),
                                    save_plots  = FALSE,
                                    save_dir    = "ModelExplain",
                                    plot_width  = 6,
                                    plot_height = 5,
                                    ...) {
  if (!requireNamespace("DALEX", quietly = TRUE))
    stop("Package 'DALEX' is required.")
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Package 'pROC' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggprism", quietly = TRUE))
    stop("Package 'ggprism' is required.")
  
  geom <- match.arg(geom)
  
  cat("-- Model Performance ------------------------------------------------------\n")
  mp <- DALEX::model_performance(explainer)
  print(mp)
  auc_val <- mp$measures$auc[1]
  
  # ---- ROC curve (custom ggplot) ----
  if (geom == "roc") {
    probs <- explainer$predict_function(explainer$model, explainer$data)
    true  <- explainer$y
    if (is.factor(true)) true <- as.numeric(true) - 1
    if (!is.numeric(true)) true <- as.numeric(as.factor(true)) - 1
    roc_obj <- pROC::roc(true, probs, levels = c(0, 1), direction = "auto", quiet = TRUE)
    roc_df <- data.frame(
      Sensitivity = roc_obj$sensitivities,
      Specificity = roc_obj$specificities
    )
    
    p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - Specificity, y = Sensitivity)) +
      ggplot2::geom_line(color = "#006d2c", linewidth = 1.2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
      ggplot2::labs(
        title    = paste("Model Performance -", explainer$label),
        subtitle = paste("AUC =", round(auc_val, 3)),
        x        = "1 - Specificity",
        y        = "Sensitivity"
      ) +
      ggplot2::coord_equal() +
      ggprism::theme_prism(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    
    print(p)
    
    if (save_plots) {
      if (is.null(save_dir) || length(save_dir) == 0 || !nzchar(save_dir)) {
        save_dir <- "ModelExplain"
      }
      save_dir <- normalizePath(save_dir, mustWork = FALSE)
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      }
      file_path <- file.path(save_dir, "performance_roc.pdf")
      ggplot2::ggsave(filename = file_path, plot = p,
                      width = plot_width, height = plot_height, ...)
      cat("ROC plot saved to:", file_path, "\n")
    }
  } else {
    # Lift or boxplot using DALEX plot
    p <- plot(mp, geom = geom) +
      ggplot2::labs(title = paste("Model Performance -", explainer$label)) +
      ggprism::theme_prism(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    
    print(p)
    
    if (save_plots) {
      if (is.null(save_dir) || length(save_dir) == 0 || !nzchar(save_dir)) {
        save_dir <- "ModelExplain"
      }
      save_dir <- normalizePath(save_dir, mustWork = FALSE)
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      }
      file_path <- file.path(save_dir, paste0("performance_", geom, ".pdf"))
      ggplot2::ggsave(filename = file_path, plot = p,
                      width = plot_width, height = plot_height, ...)
      cat("Plot saved to:", file_path, "\n")
    }
  }
  
  invisible(mp)
}


# -- 6. Unified Wrapper: run all methods and benchmark ------------------------
#' Comprehensive Feature Selection with Method Benchmarking
#'
#' Runs up to four caret-based selection methods on the same `Train_Model`
#' object, returns the union or intersection of selected features, and
#' generates a comparative Upset plot.
#' @description
#' Integrates four robust feature selection paradigms from the \pkg{caret} package 
#' to identify optimal predictive biomarker subsets. This pipeline supports automated 
#' cross-validation, parallel execution, and visual benchmarking via UpSet plots
#' 
#' The module includes the following core strategies:
#' \itemize{
#'   \item \strong{RFE (Recursive Feature Elimination):} Backward selection wrapper that 
#'     iteratively removes weaker features. By default, it uses \code{caret::rfFuncs}, 
#'     leveraging Random Forest as the underlying evaluation engine
#'   \item \strong{GA (Genetic Algorithm):} Evolutionary global optimization wrapper. 
#'     By default, it utilizes \code{caret::caretGA} (a general wrapper typically 
#'     backed by decision trees like \pkg{rpart} unless explicitly overridden with 
#'     \code{rfGA}).
#'   \item \strong{SA (Simulated Annealing):} Probabilistic metaheuristic optimization wrapper. 
#'     By default, it utilizes \code{caret::caretSA} (general wrapper, typically 
#'     backed by \pkg{rpart} unless overridden with \code{rfSA})
#'   \item \strong{SBF (Selection By Filter):} Univariate filter-based feature selection. 
#'     By default, it uses \code{caret::rfSBF}, applying a Random Forest filter across 
#'     resampling folds
#' }
#' @param object   A Train_Model or Stat S4 object.
#' @param methods  Character vector: "rfe","ga","sa","sbf" (default all four).
#'                Can also pass a named list of pre-run selection results.
#' @param combine  How to combine across methods: "union", "intersect", or "none".
#' @param rfe_args List of extra args passed to \code{FeatureSelectRFE}.
#' @param ga_args  List of extra args passed to \code{FeatureSelectGA}.
#' @param sa_args  List of extra args passed to \code{FeatureSelectSA}.
#' @param sbf_args List of extra args passed to \code{FeatureSelectSBF}.
#' @param upset_plot Logical. Draw an Upset plot of selected features.
#' @param save_plot  Logical.
#' @param save_dir   Directory.
#' @param seed       Random seed.
#'
#' @return An invisible list with:
#'   \describe{
#'   \item{results}{Named list of selection results per method.}
#'   \item{selected_features}{Character vector of final features.}
#'   \item{feature_matrix}{Binary matrix used for Upset.}
#'   \item{upset_plot}{ggplot (if \code{upset_plot = TRUE}).}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   fs <- FeatureSelectionPipeline(model, methods = c("rfe", "sbf"), combine = "union",
#'                                  rfe_args = list(sizes = c(2,4), method = "cv", number = 2),
#'                                  upset_plot = FALSE)
#'   print(fs$selected_features)
#' }
#' }
FeatureSelectionPipeline <- function(object,
                                     methods     = c("rfe", "ga", "sa", "sbf"),
                                     combine     = c("union", "intersect", "none"),
                                     rfe_args    = list(),
                                     ga_args     = list(),
                                     sa_args     = list(),
                                     sbf_args    = list(),
                                     upset_plot  = TRUE,
                                     save_plot   = FALSE,
                                     save_dir    = NULL,
                                     seed        = 825) {
  .check_fs_packages()
  combine <- match.arg(combine)
  
  if (is.null(save_dir))
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  result_list <- list()
  
  # Run each requested method
  if ("rfe" %in% methods) {
    cat("\n===== RFE =====\n")
    result_list$RFE <- do.call(FeatureSelectRFE,
      c(list(object = object, save_dir = save_dir, seed = seed), rfe_args))
  }
  if ("ga" %in% methods) {
    cat("\n===== GA =====\n")
    result_list$GA <- do.call(FeatureSelectGA,
      c(list(object = object, save_dir = save_dir, seed = seed), ga_args))
  }
  if ("sa" %in% methods) {
    cat("\n===== SA =====\n")
    result_list$SA <- do.call(FeatureSelectSA,
      c(list(object = object, save_dir = save_dir, seed = seed), sa_args))
  }
  if ("sbf" %in% methods) {
    cat("\n===== SBF =====\n")
    result_list$SBF <- do.call(FeatureSelectSBF,
      c(list(object = object, save_dir = save_dir, seed = seed), sbf_args))
  }
  
  # Compile feature lists
  feature_lists <- lapply(result_list, `[[`, "opt_vars")
  
  # Combine
  if (combine == "union") {
    selected <- unique(unlist(feature_lists))
  } else if (combine == "intersect") {
    selected <- Reduce(intersect, feature_lists)
  } else {
    selected <- unlist(feature_lists)
  }
  
  cat(sprintf("\n%s of features across methods: %d\n", combine, length(selected)))
  
  # Upset plot
  p_upset <- NULL
  if (upset_plot && length(feature_lists) >= 2) {
    if (requireNamespace("UpSetR", quietly = TRUE)) {
      mat <- UpSetR::fromList(feature_lists)
      p_upset <- UpSetR::upset(mat, sets = names(feature_lists),
                               order.by = "freq",
                               text.scale = 1.2,
                               mainbar.y.label = "Number Intersected",
                               sets.x.label   = "Number Selected")
      print(p_upset)
      if (save_plot) {
        .safe_dir(save_dir)
        pdf(file.path(save_dir, "feature_upset.pdf"), width = 8, height = 5)
        print(p_upset)
        dev.off()
      }
    } else {
      cat("Install 'UpSetR' for the Upset plot.\n")
      # Fallback: simple venn-style barplot
      all_vars <- unique(unlist(feature_lists))
      bin_mat <- sapply(feature_lists, function(v) as.integer(all_vars %in% v))
      rownames(bin_mat) <- all_vars
      p_upset <- bin_mat
    }
  }
  
  invisible(list(
    results          = result_list,
    selected_features = selected,
    feature_matrix   = if (exists("bin_mat")) bin_mat else NULL,
    upset_plot       = p_upset
  ))
}


# -- 7. Quick train-and-select: built-in importance from multiple models ------
#' Feature selection using built-in variable importance from multiple models
#'
#' Trains one or more classification models and extracts their built-in variable
#' importance measures. The top features from each model are combined via union
#' or intersection to produce a final selected feature set.
#'
#' @param object A data object containing features and a binary response.
#'   Must be compatible with internal extractor `.extract_xy()`.
#' @param models Character vector of model names supported by \pkg{caret} that
#'   have a built-in `varImp` method. Default: `c("rf", "gbm")`.
#' @param method Resampling method for `trainControl`, e.g., `"repeatedcv"`.
#' @param number Number of folds or resampling iterations.
#' @param top_n Integer, number of top features to retain per model.
#' @param combine Character, either `"union"` (default) or `"intersection"`,
#'   defining how per‑model feature sets are combined.
#' @param seed Random seed for reproducibility.
#' @param metric Character string specifying the performance metric to optimize.
#'   Default `NULL` auto‑selects: `"ROC"` for binary classification, otherwise `"Accuracy"`.
#'   User‑supplied values override the automatic choice.
#' @param tuneLength Integer, passed to `caret::train()`; if `NULL`, model
#'   defaults are used.
#' @param class_importance Specifies which class's importance to use for multi-class models.
#'   - `NULL` (default): use the first column (original behavior).
#'   - Character string: name of the class (must match a column name in `varImp` output).
#'   - Integer: column index (1-based).
#'   - `"max"`: use the maximum importance across all classes.
#'   - `"min"`: use the minimum importance across all classes.
#'   - `"all"`: use row means across all classes (discouraged if you care about specific classes).
#'   For binary classification, this parameter is ignored (only one column exists).
#' @param ... Additional arguments passed to `caret::train()` (e.g., `ntree`, `n.trees`).
#'
#' @return An invisible list with components:
#'   \item{importance_table}{A data frame showing which features were selected
#'     by each model and the final `Selected` status.}
#'   \item{selected_features}{Character vector of the final chosen features.}
#'   \item{per_model}{A named list of feature vectors per model.}
#'
#' @details
#' Models lacking a built-in `varImp` are automatically skipped with a note.
#' For `xgbTree`, verbose C++ messages are suppressed to keep console clean.
#' The function requires the \pkg{caret} and respective modelling packages.
#' When `metric = "ROC"`, the summary function `twoClassSummary` is used,
#' which requires class probabilities (`classProbs = TRUE`).
#'
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   imp <- FeatureSelectBuiltin(model, models = c("rf", "glm"), method = "cv", number = 2, top_n = 3)
#'   print(imp$selected_features)
#' }
#' }
FeatureSelectBuiltin <- function(object, models = c("rf", "gbm"), method = "repeatedcv",
                                 number = 5, top_n = 15, combine = "union", seed = 825,
                                 metric = NULL, tuneLength = NULL, class_importance = NULL,
                                 ...) {
  .check_fs_packages()
  
  # Ensure gbm is fully attached to avoid missing S3 methods inside caret
  if ("gbm" %in% models) {
    if (!requireNamespace("gbm", quietly = TRUE)) {
      stop("Package 'gbm' is required when 'gbm' is included in `models`. ",
           "Please install it with install.packages('gbm').", call. = FALSE)
    }
  }
  
  set.seed(seed)
  xy <- .extract_xy(object)
  
  # ---- Ensure response is factor and fix levels ----
  y_factor <- as.factor(xy$y)
  levels(y_factor) <- make.names(levels(y_factor))
  is_binary <- length(levels(y_factor)) == 2
  
  # ---- Fix column name conflict for response ----
  resp_name <- ".group"  # safer name
  if (resp_name %in% colnames(xy$x)) {
    # If .group already exists, use a unique name
    resp_name <- make.names(c(colnames(xy$x), "response"))[length(colnames(xy$x)) + 1]
    warning("Column name '.group' already exists; using '", resp_name, "' as response column.")
  }
  df <- cbind(xy$x, resp = y_factor)
  colnames(df)[ncol(df)] <- resp_name
  form <- as.formula(paste(resp_name, "~ ."))
  
  # ---- Determine metric ----
  if (is.null(metric)) {
    final_metric <- if (is_binary) "ROC" else "Accuracy"
  } else {
    final_metric <- metric
    if (final_metric == "ROC" && !is_binary) {
      warning("'ROC' metric requested but response is not binary; falling back to 'Accuracy'.")
      final_metric <- "Accuracy"
    }
  }
  
  # ---- Build trainControl once (independent of model) ----
  trControl_args <- list(method = method, number = number, classProbs = TRUE)
  if (final_metric == "ROC" && is_binary) {
    trControl_args$summaryFunction <- caret::twoClassSummary
  }
  trControl <- do.call(caret::trainControl, trControl_args)
  
  # ---- Check built-in varImp availability ----
  support_df <- check_varImp_availability(models)
  unsupported <- support_df$Model[!support_df$Has_BuiltIn]
  if (length(unsupported) > 0) {
    cat("Note: The following models do NOT have built-in varImp and will be skipped:\n")
    cat(paste("  -", unsupported), sep = "\n")
  }
  
  imp_list <- list()
  for (m in models) {
    if (m %in% unsupported) 
      next
    model_info <- caret::getModelInfo(m, regex = FALSE)[[1]]
    needed_pkgs <- unique(c(model_info$library))
    missing_pkg <- needed_pkgs[!sapply(needed_pkgs, requireNamespace, 
                                       quietly = TRUE)]
    if (length(missing_pkg) > 0) {
      warning("Skipping model '", m, "' because package(s) missing: ", 
              paste(missing_pkg, collapse = ", "))
      next
    }
    cat(sprintf("Training %s for importance...\n", m))
    
    # ---- Assemble train() arguments ----
    train_args <- list(
      form = form,
      data = df,
      method = m,
      trControl = trControl,
      metric = final_metric
    )
    if (!is.null(tuneLength)) train_args$tuneLength <- tuneLength
    if (m == "glm" || m == "glmStepAIC") train_args$family <- "binomial"
    extra_args <- list(...)
    if (length(extra_args) > 0) train_args <- utils::modifyList(train_args, extra_args)
    
    fit <- tryCatch({
      if (m == "xgbTree") {
        result <- NULL
        invisible(capture.output(result <- do.call(caret::train, train_args)))
        result
      } else {
        do.call(caret::train, train_args)
      }
    }, error = function(e) {
      warning("Training failed for model '", m, "': ", e$message)
      return(NULL)
    })
    
    if (is.null(fit)) 
      next
    
    imp <- tryCatch({
      caret::varImp(fit, scale = TRUE)$importance
    }, error = function(e) {
      warning("Could not extract importance for model '", 
              m, "': ", e$message)
      return(NULL)
    })
    
    if (is.null(imp)) 
      next
    
    # ---- Ensure imp is a data.frame for consistent handling ----
    if (!is.data.frame(imp)) {
      imp <- as.data.frame(imp)
    }
    
    # ---- Extract importance scores based on class_importance ----
    if (ncol(imp) == 1) {
      imp_scores <- imp[, 1]
      names(imp_scores) <- rownames(imp)
    } else {
      # Multi-class: determine which column(s) to use
      if (is.null(class_importance)) {
        # Default: first column
        imp_scores <- imp[, 1]
        names(imp_scores) <- rownames(imp)
      } else if (is.character(class_importance) && length(class_importance) == 1) {
        if (class_importance %in% c("max", "min", "all")) {
          if (class_importance == "max") {
            imp_scores <- apply(imp, 1, max, na.rm = TRUE)
          } else if (class_importance == "min") {
            imp_scores <- apply(imp, 1, min, na.rm = TRUE)
          } else { # "all"
            imp_scores <- rowMeans(imp, na.rm = TRUE)
          }
          names(imp_scores) <- rownames(imp)
        } else {
          # User specified a class name
          if (class_importance %in% colnames(imp)) {
            imp_scores <- imp[, class_importance]
            names(imp_scores) <- rownames(imp)
          } else {
            warning("Class name '", class_importance, "' not found in importance columns for model '", m, 
                    "'. Available: ", paste(colnames(imp), collapse = ", "), 
                    ". Using first column as fallback.")
            imp_scores <- imp[, 1]
            names(imp_scores) <- rownames(imp)
          }
        }
      } else if (is.numeric(class_importance) && length(class_importance) == 1) {
        idx <- as.integer(class_importance)
        if (idx >= 1 && idx <= ncol(imp)) {
          imp_scores <- imp[, idx]
          names(imp_scores) <- rownames(imp)
        } else {
          warning("Column index ", idx, " out of range for model '", m, 
                  "' (1-", ncol(imp), "). Using first column as fallback.")
          imp_scores <- imp[, 1]
          names(imp_scores) <- rownames(imp)
        }
      } else {
        warning("Invalid 'class_importance' specification for model '", m, 
                "'; using first column as fallback.")
        imp_scores <- imp[, 1]
        names(imp_scores) <- rownames(imp)
      }
    }
    
    # Sort and take top_n
    feats <- head(names(sort(imp_scores, decreasing = TRUE)), top_n)
    imp_list[[m]] <- feats
  }
  
  if (length(imp_list) == 0) 
    stop("No models could be trained or provided importance.")
  
  selected <- if (combine == "union") 
    unique(unlist(imp_list))
  else Reduce(intersect, imp_list)
  
  all_feats <- unique(unlist(imp_list))
  comb_imp <- data.frame(Feature = all_feats, row.names = all_feats)
  for (m in names(imp_list)) {
    comb_imp[[m]] <- ifelse(all_feats %in% imp_list[[m]], "Yes", "-")
  }
  comb_imp$Selected <- ifelse(all_feats %in% selected, "[OK]", "")
  
  # drop = FALSE prevents single-model errors in rowSums()
  comb_imp <- comb_imp[order(-rowSums(comb_imp[, names(imp_list), drop = FALSE] == "Yes")), ]
  
  cat(sprintf("Built-in (%s): %d features\n", combine, length(selected)))
  invisible(list(importance_table = comb_imp, selected_features = selected, 
                 per_model = imp_list))
}


# -- 8. Helper: Update Train_Model with selected features ---------------------

#' Apply Selected Features to a Train_Model Object
#'
#' Subsets the Train_Model to keep only the chosen features. If a train/test
#' split already exists on \code{object} (i.e. \code{object@split.data}
#' already has non-empty \code{training}/\code{testing} elements), that same
#' row-level split is preserved -- only its columns are subset to the
#' selected features -- rather than being wiped out. This matters because
#' feature selection is expected to be run \emph{on the training portion
#' only} (e.g. via a separate \code{Train_Model} object built from
#' \code{object@split.data$training}); silently discarding an
#' already-established split here would force the caller to re-draw a new
#' random split afterwards, which is easy to forget and previously caused
#' downstream steps (e.g. \code{preProcess()} on
#' \code{object@split.data$training}) to fail with "0 rows"/NULL errors.
#'
#' \code{split.scale.data} is subset the same way if present. Trained
#' models and results are always cleared, since they were fit on the old
#' feature set and are no longer valid.
#'
#' @param object   A Train_Model object.
#' @param features Character vector of feature names to keep.
#' @return An updated Train_Model object.
#' @export
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' # Select first 3 numeric features as example
#' features <- colnames(mtcars)[1:3]
#' model <- ApplyFeatureSelection(model, features)
#' print(dim(model@clean.df))
#' }
ApplyFeatureSelection <- function(object, features) {
  if (!inherits(object, "Train_Model"))
    stop("Object must be Train_Model.")
  
  gc <- object@group_col
  keep <- unique(c(features, gc))
  keep <- intersect(keep, colnames(object@clean.df))
  
  object@clean.df <- object@clean.df[, keep, drop = FALSE]
  object@data.df  <- object@data.df[, intersect(keep, colnames(object@data.df)), drop = FALSE]
  
  # -- Preserve an existing split (row-level) by subsetting its columns,
  # instead of discarding it. This is the safe default: it never changes
  # *which* rows are training vs. testing (no new randomness / no risk of
  # accidentally re-mixing rows across the split), it only drops columns
  # that were not selected.
  .subset_split <- function(split_list, keep_cols) {
    if (length(split_list) == 0) return(list())
    lapply(split_list, function(df) {
      if (is.null(df) || !is.data.frame(df)) return(df)
      cols <- intersect(keep_cols, colnames(df))
      df[, cols, drop = FALSE]
    })
  }
  
  had_split <- length(object@split.data) > 0 &&
    all(c("training", "testing") %in% names(object@split.data)) &&
    !is.null(object@split.data$training) && !is.null(object@split.data$testing)
  
  if (had_split) {
    object@split.data <- .subset_split(object@split.data, keep)
    if (length(object@split.scale.data) > 0) {
      object@split.scale.data <- .subset_split(object@split.scale.data, keep)
    }
    cat("Existing train/test split preserved (rows unchanged); columns subset to selected features.\n")
  } else {
    # Nothing to preserve -- keep prior behaviour of resetting to empty lists
    # so a fresh split is drawn downstream.
    object@split.data       <- list()
    object@split.scale.data <- list()
  }
  
  # Trained models/results were fit on the previous feature set and are no
  # longer valid -- these are always cleared regardless of the split.
  object@train.models      <- list()
  object@all.results       <- list()
  object@best.model.result <- list()
  
  cat(sprintf("Applied feature selection: %d features kept.\n", length(keep) - 1))
  cat("Downstream slots (models, results) have been reset", 
      if (!had_split) "; split.data/split.scale.data reset as well (no prior split found).\n" else ".\n")
  return(object)
}


#' Check Built-in Variable Importance Availability for caret Models
#'
#' @param model_names Character vector of caret model names.
#' @return A data frame with model name and whether varImp is supported.
#' @export
check_varImp_availability <- function(model_names) {
  if (!requireNamespace("caret", quietly = TRUE)) stop("caret required.")
  result <- data.frame(Model = model_names, Has_BuiltIn = NA)
  for (i in seq_along(model_names)) {
    if (model_names[i] == "rpart") {
      result$Has_BuiltIn[i] <- FALSE
      next
    }
    info <- caret::getModelInfo(model_names[i], regex = FALSE)[[1]]
    has_imp <- !is.null(info$varImp)
    known <- c("rf", "C5.0", "glmnet")
    result$Has_BuiltIn[i] <- has_imp || (model_names[i] %in% known)
  }
  return(result)
}

