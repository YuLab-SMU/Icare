# ============================================================================
# Internal Utility Functions
# These are helper functions used across modules
# All start with dot (.) to indicate they are internal
# ============================================================================

#' @importFrom stats aggregate approx ave binomial complete.cases confint cor
#'   density glm IQR kmeans kruskal.test lm loess model.matrix pchisq prcomp
#'   reorder setNames
#' @importFrom utils capture.output install.packages packageVersion read.csv
#' @keywords internal
NULL

utils::globalVariables(c("BIC", "time"))
#' All Available caret Models Metadata
#'
#' A data frame containing metadata for all classification and regression models
#' available in the caret package, including model labels, required libraries,
#' task type, tuning parameters, and tags.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{model}{Character, caret model method name (e.g., "rf", "glmnet").}
#'   \item{label}{Character, human-readable model label.}
#'   \item{library}{Character, comma-separated required packages.}
#'   \item{type}{Character, either "Classification" or "Regression".}
#'   \item{parameters}{Character, comma-separated tuning parameter names.}
#'   \item{tags}{Character, comma-separated tags for model characteristics.}
#' }
#'
#' @source Generated from the caret package model registry.
#'
#' @examples
#' \dontrun{
#' data("allmodel")
#' head(allmodel)
#' }
"allmodel"

#' @keywords internal
.check_class <- function(object, allowed) {
  if (!inherits(object, allowed))
    stop("Object must be one of: ", paste(allowed, collapse = ", "),
         ". Got: ", class(object)[1])
}

#' @keywords internal
.require_pkgs <- function(pkgs, hint = NULL) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    inst <- sprintf("install.packages('%s')", miss)
    stop(sprintf(
      "This functionality requires the package%s: %s. Please install %s%s",
      if (length(miss) > 1) "s" else "",
      paste(sprintf("'%s'", miss), collapse = ", "),
      paste(inst, collapse = " and "),
      if (!is.null(hint)) sprintf(" (%s)", hint) else ""
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.subset_df <- function(df, rows = NULL, cols = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!is.null(rows)) df <- df[rownames(df) %in% rows, , drop = FALSE]
  if (!is.null(cols)) df <- df[, colnames(df) %in% cols,  drop = FALSE]
  df
}

#' @keywords internal
.safe_rbind <- function(df1, df2) {
  if (is.null(df1) || nrow(df1) == 0) return(df2)
  if (is.null(df2) || nrow(df2) == 0) return(df1)
  shared <- intersect(colnames(df1), colnames(df2))
  if (length(shared) == 0) {
    warning("No shared columns; returning first data frame only.")
    return(df1)
  }
  if (length(shared) < max(ncol(df1), ncol(df2)))
    message(length(shared), " shared columns kept.")
  out <- rbind(df1[, shared, drop = FALSE], df2[, shared, drop = FALSE])
  rownames(out) <- make.unique(rownames(out))
  out
}

#' @keywords internal
.clean_symbol_values <- function(data) {
  for(col in colnames(data)) {
    if(is.character(data[[col]])) {
      if(any(grepl("^[<>]", data[[col]]), na.rm = TRUE)) {
        data[[col]] <- gsub("[<>]", "", data[[col]])
        data[[col]] <- as.numeric(data[[col]])
        warning(paste("Removed >/< symbols from column", col))
      }
    }
  }
  return(data)
}

#' @keywords internal
.prepare_data <- function(data, data_name) {
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame())
  }
  
  if (anyDuplicated(rownames(data))) {
    warning(paste("Duplicate row names in", data_name, "made unique."))
    rownames(data) <- make.unique(rownames(data))
  }
  
  if (anyDuplicated(colnames(data))) {
    warning(paste("Duplicate column names in", data_name, "made unique."))
    colnames(data) <- make.unique(colnames(data))
  }
  
  if (is.null(colnames(data))){
    stop(paste(data_name, "is missing column names."))
  }
  
  return(data)
}

#' @keywords internal
.ensure_numeric_data <- function(data, data_name, convert_factors = TRUE) {
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame())
  }
  
  row_names <- rownames(data)
  
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x)) {
      if (convert_factors) {
        levels <- levels(x)
        if (length(levels) == 2) {
          return(as.numeric(x) - 1)
        } else {
          warning("Factor with >2 levels in ", data_name)
          return(as.numeric(x))
        }
      } else {
        return(x)
      }
    } else if (is.character(x)) {
      if (convert_factors) {
        x <- as.factor(x)
        levels <- levels(x)
        if (length(levels) == 2) {
          return(as.numeric(x) - 1)
        } else {
          return(as.numeric(x))
        }
      } else {
        return(x)
      }
    } else if (is.numeric(x)) {
      return(x)
    } else {
      warning("Column in ", data_name, " cannot be converted to numeric.")
      return(x)
    }
  }))
  
  rownames(data) <- row_names
  
  # Remove all-NA columns
  na_columns <- colnames(data)[apply(data, 2, function(x) all(is.na(x)))]
  if (length(na_columns) > 0) {
    warning("Removing all-NA columns from ", data_name, ": ", 
            paste(na_columns, collapse = ", "))
    data <- data[, !colnames(data) %in% na_columns]
  }
  
  return(data)
}

#' Remove constant (zero-variance) columns from a data frame or matrix
#'
#' Identifies and removes columns that have only one distinct value (including NA).
#'
#' @param data A data frame or matrix.
#' @param na.rm Logical. If TRUE, NA values are ignored when checking constancy
#'   (i.e., a column with only one non-NA value and the rest NA is considered constant).
#'   If FALSE (default), NA is treated as a distinct value.
#' @param verbose Logical. Print messages about removed columns (default TRUE).
#' @param allow_empty Logical. If TRUE, allow returning an empty data frame when
#'   all columns are constant (default FALSE, which raises an error).
#'
#' @return A data frame or matrix (same class as input) with constant columns removed.
#' @export
#'
#' @examples
#' df <- data.frame(a = 1:5, b = rep(2, 5), c = c(1, NA, 1, 1, 1))
#' remove_constant_columns(df)                 # removes column 'b'
#' remove_constant_columns(df, na.rm = TRUE)   # also removes column 'c'
#' remove_constant_columns(df, verbose = FALSE)
remove_constant_columns <- function(data,
                                    na.rm   = FALSE,
                                    verbose = TRUE,
                                    allow_empty = FALSE) {
  
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix.")
  }
  
  # Preserve original class for return
  original_class <- class(data)
  
  # Convert to data frame for uniform handling
  if (is.matrix(data)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE)
    was_matrix <- TRUE
  } else {
    df <- data
    was_matrix <- FALSE
  }
  
  if (ncol(df) == 0) {
    if (verbose) message("No columns to check.")
    return(data)
  }
  
  # Check each column for constancy
  is_constant <- sapply(df, function(col) {
    if (na.rm) {
      non_na <- col[!is.na(col)]
      if (length(non_na) == 0) {
        # All NA -> considered constant
        return(TRUE)
      }
      length(unique(non_na)) == 1
    } else {
      length(unique(col)) == 1
    }
  })
  
  # If all columns are constant and empty not allowed, error
  if (all(is_constant) && !allow_empty) {
    stop("All columns are constant. Set allow_empty = TRUE to return an empty object.")
  }
  
  # Remove constant columns
  if (any(is_constant)) {
    removed <- names(is_constant)[is_constant]
    if (verbose) {
      cat("Removing constant columns:", paste(removed, collapse = ", "), "\n")
    }
    df <- df[, !is_constant, drop = FALSE]
  } else {
    if (verbose) cat("No constant columns found.\n")
  }
  
  # Restore original format (matrix or data frame)
  if (was_matrix && ncol(df) > 0) {
    result <- as.matrix(df)
    rownames(result) <- rownames(data)
    colnames(result) <- colnames(df)
  } else {
    result <- df
    if (is.data.frame(result) && !is.data.frame(data)) {
      result <- as.data.frame(result, stringsAsFactors = FALSE)
      rownames(result) <- rownames(data)
    }
  }
  
  return(result)
}


#' Example Stat Object for Testing
#'
#' A \code{Stat} S4 object created from the Bacteremia public dataset.
#' Used for demonstration and testing of statistical analysis functions.
#'
#' @details
#' This object was created from the first 1000 rows of the 
#' \href{https://zenodo.org/records/7554815}{Bacteremia public dataset}.
#' The \code{BloodCulture} column is used as the grouping variable.
#'
#' The object contains:
#' \itemize{
#'   \item \code{@raw.data}: Raw data with BloodCulture as grouping column
#'   \item \code{@clean.data}: Processed numeric data
#'   \item \code{@info.data}: Metadata including BloodCulture group
#'   \item \code{@group_col}: "BloodCulture"
#' }
#'
#' @format A \code{Stat} S4 object with slots:
#' \describe{
#'   \item{raw.data}{Original data frame}
#'   \item{clean.data}{Cleaned numeric data matrix}
#'   \item{info.data}{Metadata data frame}
#'   \item{group_col}{Character, grouping column name ("BloodCulture")}
#'   \item{...}{Additional slots for analysis results}
#' }
#'
#' @source \url{https://zenodo.org/records/7554815}
#' @keywords datasets
"stat_obj_test"


#' Global variables used in non-standard evaluation
#' @keywords internal
#' @noRd
utils::globalVariables(c(
  ".","label", "Prediction", "y_label", "auc_lower", "auc_upper", "allmodel",
  "::<-",
  "AUC",
  "AUC_mean",
  "AUC_se",
  "Actual",
  "Algorithm",
  "Best_Score",
  "C_Index",
  "CI_lower",
  "CI_upper",
  "Category",
  "Class",
  "Cluster",
  "Component",
  "Confidence",
  "Count",
  "Dataset",
  "Dim1",
  "Dim2",
  "Dimension 1",
  "Dimension 2",
  "Estimate",
  "FPR",
  "Facet_Median",
  "Feature",
  "Feature_Removed",
  "Fill",
  "Freq",
  "Frequency",
  "G",
  "Generation",
  "Group",
  "HR",
  "HR_95CI",
  "HR_label",
  "ID",
  "Importance",
  "Imputation",
  "Index",
  "Iteration",
  "Label",
  "Lower",
  "Mean",
  "Mean_AUC",
  "Method",
  "Method1",
  "Method2",
  "Metric",
  "Missing_Percentage",
  "Model",
  "NetBenefit",
  "New",
  "Normalization",
  "NumFeatures",
  "OUT_DIR",
  "Overall",
  "Overlap",
  "PC1",
  "PC2",
  "P_Value",
  "P_label",
  "P_value",
  "Parameter",
  "Pct",
  "Performance",
  "Performance_Drop",
  "Predicted",
  "Predicted_Subtype",
  "Probability",
  "Proportion",
  "Ref",
  "Risk",
  "SD",
  "SE",
  "Sample",
  "Score",
  "Selected",
  "Selected_Factor",
  "Selected_Status",
  "Sensitivity",
  "Sig",
  "Specificity",
  "Status",
  "Strategy",
  "Subgroup",
  "TPR",
  "Threshold",
  "Time",
  "TimePoint",
  "Upper",
  "Value",
  "Variable",
  "Variables",
  "bin",
  "change",
  "cluster_name",
  "contribution",
  "correct",
  "dif",
  "dropout_loss",
  "feat_scaled",
  "feature",
  "feature_value",
  "feature_value_num",
  "fill_col",
  "fill_group",
  "fpr",
  "group",
  "groups",
  "id",
  "logFC",
  "log_val",
  "lower",
  "mean_abs",
  "mean_pred",
  "mean_shap",
  "med",
  "mlr_learners",
  "n_features",
  "neg_log10p",
  "new_sens",
  "new_spec",
  "nri_type",
  "nri_value",
  "obs_rate",
  "observation",
  "observed",
  "outcome",
  "perc",
  "predicted",
  "prob",
  "rainbow",
  "ref_sens",
  "ref_spec",
  "sample_id",
  "se_loss",
  "shap",
  "shap_value",
  "significance",
  "sil_width",
  "stratum",
  "tAUC",
  "target_group",
  "threshold",
  "times",
  "tpr",
  "truth",
  "upper",
  "value",
  "variable",
  "x",
  "y",
  "y_position",
  "yhat","cluster", "n_high_risk_lower", "n_high_risk_upper",
  "n_high_risk", "n_events_lower", "n_events_upper",
  "n_events_high", "r"
))

#' Match factor levels of two data frames
#' @keywords internal
match_factor_levels <- function(data, ref) {
  for (col in intersect(colnames(data), colnames(ref))) {
    if (is.factor(ref[[col]])) {
      data[[col]] <- factor(data[[col]], levels = levels(ref[[col]]))
    }
  }
  data
}



.onAttach <- function(libname, pkgname) {
  required_version <- "1.7.7.1"
  installed_version <- tryCatch(
    as.character(utils::packageVersion("xgboost")),
    error = function(e) NA
  )
  
  if (is.na(installed_version) || installed_version != required_version) {
    packageStartupMessage(
      "This package requires xgboost == ", required_version,
      " (", if (is.na(installed_version)) "not found" else paste0("currently installed: ", installed_version), ").\n\n",
      "Please install the required version by running:\n\n",
      '  url <- "https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_',
      required_version, '.tar.gz"\n',
      '  dest <- file.path(tempdir(), "xgboost_', required_version, '.tar.gz")\n',
      "  download.file(url, destfile = dest, mode = \"wb\")\n",
      "  install.packages(dest, repos = NULL, type = \"source\")\n\n",
      "After installation, please restart R before loading this package."
    )
  }
}


#' Get model lists for benchmarking (with multi-index support)
#'
#' @param preset Character or character vector. Predefined sets:
#'   - Basic/Core/Extended/Advanced/All (mixed types)
#'   - all_linear, all_tree, all_ensemble, all_nn, all_svm,
#'     all_da, all_bayes, all_rule, all_proto, all_regularized,
#'     all_feature, all_robust, all_ordinal, all_special
#'   You can also pass a vector of model names directly.
#' @param include_tags,exclude_tags As before, for custom filtering.
#' @param return_df Logical; return data.frame instead of character vector.
#' @return Character vector or data.frame.
#' @export
get_models <- function(preset = NULL,
                       include_tags = NULL,
                       exclude_tags = NULL,
                       return_df = FALSE) {
  
  if (!exists("allmodel")) stop("Data 'allmodel' not found.")
  df <- allmodel
  
  # If preset is provided, handle both old and new preset names
  if (!is.null(preset)) {
    # Define all presets (old + new)
    presets_list <- list(
      # Old ones (mixed)
      basic = c("glm", "lda", "knn", "nb", "rpart"),
      core = c("glm", "lda", "knn", "nb", "rpart", "glmnet", "rf", "svmRadial", "C5.0"),
      extended = c("glm", "lda", "knn", "nb", "rpart", "glmnet", "rf", "xgbTree", "svmRadial",
                   "C5.0", "C5.0Rules", "JRip", "PART", "nnet", "treebag", "earth", "qda", "ridge", "pls"),
      advanced = c("glm", "lda", "knn", "nb", "rpart", "glmnet", "rf", "xgbTree", "svmRadial",
                   "C5.0", "C5.0Rules", "JRip", "PART", "nnet", "treebag", "earth", "qda", "ridge", "pls",
                   "avNNet", "bagEarth", "gbm", "fda", "sparseLDA", "pcaNNet", "rda", "sda", "rpart2", "cforest"),
      all = unique(df$model),
      
      # New tag-based presets (mapped to include_tags)
      all_linear = list(include = c("Linear Classifier", "Linear Regression", "Generalized Linear Model", "Logistic Regression")),
      all_tree = list(include = c("Tree-Based Model", "Random Forest", "Gradient Boosting Machines", "CART")),
      all_ensemble = list(include = c("Ensemble Model", "Bagging", "Boosting", "Random Forest")),
      all_nn = list(include = "Neural Network"),
      all_svm = list(include = c("Support Vector Machines", "Kernel Method")),
      all_da = list(include = c("Discriminant Analysis", "Linear Discriminant Analysis", "Quadratic Discriminant Analysis")),
      all_bayes = list(include = c("Bayesian Model", "Gaussian Process")),
      all_rule = list(include = "Rule-Based Model"),
      all_proto = list(include = c("Prototype Models", "Nearest Neighbors")),
      all_regularized = list(include = c("L1 Regularization", "L2 Regularization", "Ridge", "Lasso")),
      all_feature = list(include = c("Implicit Feature Selection", "Feature Extraction", "Feature Selection Wrapper")),
      all_robust = list(include = c("Robust Model", "Quantile Regression", "Robust Methods")),
      all_ordinal = list(include = c("Ordinal Outcomes", "Two Class Only", "Cost Sensitive Learning")),
      all_special = list(include = c("String Kernel", "Text Mining", "Self-Organising Maps"))
    )
    
    if (length(preset) == 1 && preset %in% names(presets_list)) {
      p <- presets_list[[preset]]
      if (is.character(p)) {
        # Old preset: direct model names
        df <- df[df$model %in% p, ]
      } else if (is.list(p) && "include" %in% names(p)) {
        # New preset: based on tags
        include_tags <- c(include_tags, p$include)
        # We'll apply include_tags filtering later
      }
    } else {
      # Treat preset as a character vector of model names
      df <- df[df$model %in% preset, ]
    }
  }
  
  # Apply tag filters (if any)
  if (!is.null(include_tags)) {
    keep <- sapply(df$tags, function(tag_str) {
      tags_vec <- strsplit(tag_str, ",\\s*")[[1]]
      all(include_tags %in% tags_vec)
    })
    df <- df[keep, ]
  }
  if (!is.null(exclude_tags)) {
    drop <- sapply(df$tags, function(tag_str) {
      tags_vec <- strsplit(tag_str, ",\\s*")[[1]]
      any(exclude_tags %in% tags_vec)
    })
    df <- df[!drop, ]
  }
  
  if (nrow(df) == 0) {
    warning("No models match the specified criteria.")
    return(if (return_df) df else character(0))
  }
  if (return_df) return(df) else return(unique(df$model))
}

#' List all available presets (both old and new)
#'
#' @return Prints a table with preset names, type, and description.
#' @export
list_presets <- function() {
  old_presets <- c("basic", "core", "extended", "advanced", "all")
  new_presets <- c("all_linear", "all_tree", "all_ensemble", "all_nn", "all_svm",
                   "all_da", "all_bayes", "all_rule", "all_proto", "all_regularized",
                   "all_feature", "all_robust", "all_ordinal", "all_special")
  descriptions <- c(
    "basic" = "5 baseline models (glm, lda, knn, nb, rpart)",
    "core" = "10 key models covering major paradigms",
    "extended" = "19 models, adding rules, MARS, QDA, etc.",
    "advanced" = "29 models with more ensembles and regularized variants",
    "all" = "All classification models in the dataset",
    "all_linear" = "All models with Linear Classifier/Regression tags",
    "all_tree" = "All tree-based models (CART, RF, GBM, etc.)",
    "all_ensemble" = "All ensemble models (bagging, boosting, RF)",
    "all_nn" = "All neural network models",
    "all_svm" = "All support vector machines and kernel methods",
    "all_da" = "All discriminant analysis models",
    "all_bayes" = "All Bayesian models",
    "all_rule" = "All rule-based models",
    "all_proto" = "All prototype/distance-based models",
    "all_regularized" = "All models with L1/L2 regularization",
    "all_feature" = "All models with feature selection/extraction",
    "all_robust" = "All robust/quantile regression models",
    "all_ordinal" = "All models handling ordinal outcomes or cost-sensitive learning",
    "all_special" = "Special models (text kernels, SOM, etc.)"
  )
  cat("Available presets:\n\n")
  cat("Old (mixed) presets:\n")
  for (p in old_presets) {
    cat(sprintf("  %-12s: %s\n", p, descriptions[p]))
  }
  cat("\nNew (tag-based) presets:\n")
  for (p in new_presets) {
    cat(sprintf("  %-12s: %s\n", p, descriptions[p]))
  }
  invisible(list(old = old_presets, new = new_presets))
}


#' Create a subdirectory under the global output root
#'
#' This function creates a subdirectory (or nested subdirectories) under the
#' global output root set by \code{\link{set_output_root}}. If the root has not
#' been set, it defaults to a temporary directory for the current R session.
#'
#' @param ... Character strings specifying subdirectory path components.
#' @return The full path to the created directory (invisibly).
#'
#' @examples
#' \dontrun{
#' set_output_root("./my_results")
#' sub_dir("figures", "final")
#' sub_dir("data", "processed")
#' }
#'
#' @seealso \code{\link{set_output_root}}, \code{\link{get_output_root}}
#' @export
sub_dir <- function(...) {
  base <- get_output_root()
  if (is.null(base)) {
    # Fallback to a session‑specific temporary directory
    base <- file.path(tempdir(), "icare_output")
    dir.create(base, showWarnings = FALSE, recursive = TRUE)
  }
  d <- file.path(base, ...)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  invisible(d)
}
