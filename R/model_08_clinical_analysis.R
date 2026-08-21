# =============================================================================
# clinical_analysis.R
# Comprehensive Clinical Analysis Module
# Subgroups - Confounders - Thresholds - Decision Curves
# =============================================================================
# Supports: Train_Model, caret train, ensemble, fine-tuned models

#' @keywords internal
.check_clinical_pkgs <- function() {
  required <- c("ggplot2", "dplyr", "tidyr", "wesanderson", "ggprism",
                "pROC", "caret", "reshape2", "scales", "ggrepel",
                "nricens", "corrplot")
  missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Please install them.")
  }
  invisible(TRUE)
}

#' @keywords internal
.predict_probs <- function(model_obj, newdata, model_name = NULL,
                           positive_class = NULL) {
  if (inherits(model_obj, "Train_Model")) {
    if (is.null(model_name)) {
      best_model <- model_obj@best.model.result$model
      if (is.null(best_model)) best_model <- model_obj@train.models[[1]]
    } else if (model_name == "ensemble") {
      ens <- model_obj@best.model.result$ensemble
      if (is.null(ens)) stop("No ensemble found.")
      return(ens$predict_fn(newdata))
    } else {
      best_model <- model_obj@train.models[[model_name]]
    }
    prob_mat <- predict(best_model, newdata, type = "prob")

    if (is.matrix(prob_mat) || is.data.frame(prob_mat)) {
      if (ncol(prob_mat) == 1) {
        return(prob_mat[, 1])
      }
      if (is.null(positive_class)) {
        positive_class <- colnames(prob_mat)[2]
        warning("positive_class not specified; using second column: ", positive_class)
        return(prob_mat[, positive_class])
      }
      if (!positive_class %in% colnames(prob_mat)) {
        warning("Column '", positive_class, "' not found. Using second column.")
        return(prob_mat[, 2])
      }
      return(prob_mat[, positive_class])
    } else {
      return(prob_mat)
    }
    
  } else if (inherits(model_obj, "train")) {
    prob_mat <- predict(model_obj, newdata, type = "prob")

    if (is.matrix(prob_mat) || is.data.frame(prob_mat)) {
      if (ncol(prob_mat) == 1) return(prob_mat[, 1])
      if (is.null(positive_class)) {
        positive_class <- colnames(prob_mat)[2]
        warning("positive_class not specified; using second column: ", positive_class)
        return(prob_mat[, positive_class])
      }
      if (!positive_class %in% colnames(prob_mat)) {
        warning("Column '", positive_class, "' not found. Using second column.")
        return(prob_mat[, 2])
      }
      return(prob_mat[, positive_class])
    } else {
      return(prob_mat)
    }
  } else {
    stop("model_obj must be a Train_Model or caret train object.")
  }
}

#' Attach Clinical Data and External Validation Set to Train_Model Object
#'
#' @param object            A Train_Model object.
#' @param clinical_data     A data frame with clinical variables (rownames must match sample IDs).
#' @param external_data     Optional matrix/data.frame for external validation dataset (expression/feature matrix).
#' @param external_clinical Optional data frame for external validation dataset clinical metadata.
#' @return The updated Train_Model object.
#' @export
AttachClinicalData <- function(object, 
                               clinical_data     = NULL, 
                               external_data     = NULL, 
                               external_clinical = NULL) {
  if (!inherits(object, "Train_Model"))
    stop("object must be a Train_Model.")
  
  if (!is.null(clinical_data)) {
    if (!is.data.frame(clinical_data)) stop("clinical_data must be a data frame.")
    object@process.info$clinical_data <- clinical_data
    cat("Clinical data attached to model_obj@process.info$clinical_data.\n")
  }
  
  if (!is.null(external_data)) {
    object@process.info$external_data <- as.data.frame(external_data)
    cat("External validation data attached to model_obj@process.info$external_data.\n")
  }
  
  if (!is.null(external_clinical)) {
    if (!is.data.frame(external_clinical)) stop("external_clinical must be a data frame.")
    object@process.info$external_clinical <- external_clinical
    cat("External clinical data attached to model_obj@process.info$external_clinical.\n")
  }
  
  return(object)
}

#' Plot Clinical Correlation with Model Predictions
#'
#' Computes and visualizes associations between model predictions
#' (probabilities) and clinical variables from a given dataset.
#'
#' Unlike a naive approach that coerces ALL clinical variables to numeric
#' via \code{as.numeric(as.factor(...))} before computing Spearman
#' correlation, this function distinguishes variable types:
#'
#' \itemize{
#'   \item \strong{Continuous / ordinal variables} (numeric, or factors
#'         explicitly marked \code{ordered = TRUE}) are included in the
#'         Spearman correlation matrix, since a monotonic relationship is
#'         a meaningful concept for them.
#'   \item \strong{Unordered categorical (nominal) variables} (character or
#'         unordered factor, e.g. sex, tumor site, treatment arm) are
#'         analyzed separately using a group-difference test
#'         (Wilcoxon rank-sum for 2 groups, Kruskal-Wallis for 3+ groups)
#'         comparing the predicted probability across categories, and
#'         visualized with boxplots rather than forced into a correlation
#'         coefficient.
#' }
#'
#' This avoids assigning nominal categories an artificial numeric order,
#' which would make the resulting correlation coefficient and p-value
#' uninterpretable.
#'
#' @param model_obj A trained model object compatible with \code{predict} method
#'   (e.g., \code{train}, \code{glm}, \code{randomForest}) that supports
#'   \code{type = "prob"} for classification.
#' @param clinical_data A data frame containing clinical variables and a sample
#'   identifier column that matches the row names of \code{newdata}. All columns
#'   except the identifier will be used as clinical covariates.
#' @param newdata A data frame containing the feature data (predictors) used to
#'   generate predictions. The row names (or a column) must match those in
#'   \code{clinical_data} for alignment.
#' @param model_name Optional character string naming the model (used in plot
#'   subtitle). If \code{NULL}, the model's \code{method} or class is used.
#' @param dataset_type Character string indicating the dataset origin for plot
#'   labelling. Must be one of \code{"testing"}, \code{"training"}, or
#'   \code{"external"}. Default is \code{"testing"}.
#' @param ordinal_vars Character vector of column names in \code{clinical_data}
#'   that are ordinal (i.e. have a meaningful order) even though they may be
#'   stored as character/unordered factor. These will be treated as ordered
#'   and included in the Spearman correlation panel. Default \code{NULL}.
#' @param save_plot Logical; if \code{TRUE}, saves the plot(s) as PDF to
#'   \code{save_dir}. Default is \code{FALSE}.
#' @param save_dir Character string specifying the directory where the PDF(s)
#'   should be saved. Required if \code{save_plot = TRUE}.
#' @param palette_name Character string naming the Wes Anderson palette to use
#'   for the color gradient / boxplot fills. Default is \code{"Royal1"}. See
#'   \code{wesanderson::wes_palette()} for available options.
#' @param ... Additional arguments passed to \code{corrplot::corrplot()} for
#'   fine-tuning the correlation plot appearance (e.g., \code{tl.cex}).
#'
#' @return Invisibly returns a list with two elements:
#'   \describe{
#'     \item{\code{spearman}}{A list with \code{matrix} (the Spearman
#'       correlation matrix among Prediction + continuous/ordinal variables)
#'       and \code{p.values} (the corresponding p-value matrix). \code{NULL}
#'       if no continuous/ordinal variables are present.}
#'     \item{\code{categorical}}{A data frame summarizing, for each nominal
#'       variable, the test used (Wilcoxon or Kruskal-Wallis), the test
#'       statistic, and the p-value. \code{NULL} if no nominal variables are
#'       present.}
#'   }
#'
#' @section Required Packages:
#' \pkg{corrplot}, \pkg{wesanderson}, \pkg{ggplot2}, and internal helper
#' functions (\code{.check_clinical_pkgs()}, \code{.align_clinical_and_newdata()},
#' \code{.predict_probs()}). These dependencies are checked automatically.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   set.seed(123)
#'   idx <- sample(1:nrow(mtcars), 20)
#'   model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#'   trained <- ModelTrainAnalysis(model, methods = c("glm"), 
#'   control = list(method = "cv", number = 3), save_plots = FALSE)
#'   clin <- data.frame(row.names = rownames(mtcars), age = runif(nrow(mtcars), 30, 80))
#'   PlotClinicalCorrelation(trained, clinical_data = clin,
#'    newdata = mtcars, dataset_type = "testing", save_plot = FALSE)
#' }
#' }
PlotClinicalCorrelation <- function(model_obj,
                                    clinical_data = NULL,
                                    newdata       = NULL,
                                    model_name    = NULL,
                                    dataset_type  = c("testing", "training", "external"),
                                    ordinal_vars  = NULL,
                                    save_plot     = FALSE,
                                    save_dir      = NULL,
                                    palette_name  = "Royal1",
                                    ...) {
  .check_clinical_pkgs()
  dataset_type <- match.arg(dataset_type)
  
  if (save_plot && (is.null(save_dir) || !nzchar(save_dir))) {
    stop("`save_dir` must be provided when `save_plot = TRUE`.")
  }
  
  aligned <- .align_clinical_and_newdata(model_obj, clinical_data, newdata, dataset_type)
  clinical_data <- aligned$clinical_data
  newdata       <- aligned$newdata
  
  probs <- .predict_probs(model_obj, newdata, model_name)
  
  df <- clinical_data
  
  # --- Classify variables ---------------------------------------------
  is_numeric_col   <- vapply(df, is.numeric, logical(1))
  is_ordered_col   <- vapply(df, is.ordered, logical(1))
  is_marked_ordinal <- names(df) %in% ordinal_vars
  
  keep_for_spearman <- is_numeric_col | is_ordered_col | is_marked_ordinal
  nominal_vars <- names(df)[!keep_for_spearman]
  ordinal_numeric_vars <- names(df)[keep_for_spearman]
  
  if (length(nominal_vars) > 0) {
    message(
      "The following variables are treated as unordered categorical and ",
      "excluded from the Spearman matrix (analyzed via group-difference ",
      "test instead): ", paste(nominal_vars, collapse = ", ")
    )
  }
  
  result <- list(spearman = NULL, categorical = NULL)
  
  # --- Part 1: Spearman correlation for continuous/ordinal variables ---
  if (length(ordinal_numeric_vars) > 0) {
    df_cont <- df[, ordinal_numeric_vars, drop = FALSE]
    # Ordered factors -> integer codes reflecting their defined level order
    for (v in names(df_cont)) {
      if (is.ordered(df_cont[[v]])) {
        df_cont[[v]] <- as.numeric(df_cont[[v]])
      } else if (!is.numeric(df_cont[[v]])) {
        # marked ordinal but stored as plain character/factor: warn and coerce
        # using factor level order as given (user is responsible for order)
        warning(
          "Variable '", v, "' was listed in `ordinal_vars` but is not an ",
          "ordered factor or numeric. Coercing using its current level ",
          "order via as.factor(); verify this order is correct."
        )
        df_cont[[v]] <- as.numeric(as.factor(df_cont[[v]]))
      }
    }
    df_cont$Prediction <- probs
    
    cor_res <- corrplot::cor.mtest(df_cont, conf.level = 0.95)
    M <- cor(df_cont, method = "spearman", use = "pairwise.complete.obs")
    
    cols <- colorRampPalette(
      c(wesanderson::wes_palette(palette_name, 2, type = "discrete")[1],
        "white",
        wesanderson::wes_palette(palette_name, 2, type = "discrete")[2])
    )(200)
    
    draw_corr <- function() {
      corrplot::corrplot(
        M, method = "square", type = "lower", tl.col = "black", tl.cex = 0.8,
        diag = FALSE, p.mat = cor_res$p, sig.level = c(0.001, 0.01, 0.05),
        insig = "label_sig", pch.cex = 1, pch.col = "grey20", col = cols,
        title = paste0("Prediction vs Continuous/Ordinal Clinical Variables (", dataset_type, ")"),
        mar = c(0, 0, 2, 0), ...
      )
    }
    
    draw_corr()
    
    if (save_plot) {
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      pdf(file.path(save_dir, paste0("clinical_correlation_spearman_", dataset_type, ".pdf")),
          width = 8, height = 8)
      draw_corr()
      dev.off()
    }
    
    result$spearman <- list(matrix = M, p.values = cor_res$p)
  } else {
    message("No continuous/ordinal clinical variables found; skipping Spearman panel.")
  }
  
  # --- Part 2: Group-difference test + boxplots for nominal variables --
  if (length(nominal_vars) > 0) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required to plot nominal-variable boxplots.")
    }
    
    pal_vec <- wesanderson::wes_palette(palette_name)
    
    cat_results <- lapply(nominal_vars, function(v) {
      grp <- as.factor(df[[v]])
      plot_df <- data.frame(Prediction = probs, Group = grp)
      plot_df <- plot_df[!is.na(plot_df$Group), ]
      
      n_groups <- nlevels(droplevels(plot_df$Group))
      if (n_groups < 2) {
        return(data.frame(variable = v, test = NA, statistic = NA, p.value = NA))
      }
      
      if (n_groups == 2) {
        test_name <- "Wilcoxon rank-sum"
        test_res <- wilcox.test(Prediction ~ Group, data = plot_df)
      } else {
        test_name <- "Kruskal-Wallis"
        test_res <- kruskal.test(Prediction ~ Group, data = plot_df)
      }
      
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = Prediction, fill = Group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.85) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
        ggplot2::scale_fill_manual(values = colorRampPalette(pal_vec)(n_groups)) +
        ggplot2::labs(
          title = paste0("Prediction by ", v, " (", dataset_type, ")"),
          subtitle = paste0(
            test_name, " p = ",
            formatC(unname(test_res$p.value), format = "g", digits = 3)
          ),
          x = v, y = "Predicted probability"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
      
      print(p)
      
      if (save_plot) {
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
        ggplot2::ggsave(
          filename = file.path(save_dir, paste0("clinical_boxplot_", v, "_", dataset_type, ".pdf")),
          plot = p, width = 6, height = 5
        )
      }
      
      data.frame(
        variable  = v,
        test      = test_name,
        statistic = unname(test_res$statistic),
        p.value   = unname(test_res$p.value)
      )
    })
    
    result$categorical <- do.call(rbind, cat_results)
  } else {
    message("No unordered categorical clinical variables found; skipping group-difference panel.")
  }
  
  invisible(result)
}

#' Subgroup Performance Forest Plot (table style, via forestploter)
#'
#' Computes AUC with 95% CI within each subgroup defined by one or more
#' categorical clinical variables, and renders a single publication-style
#' forest table (variable | N | forest graphic | estimate (95% CI) | p)
#' similar to clinical hazard-ratio tables, using the \pkg{forestploter}
#' package. When \code{subgroup_var} has length > 1, all variables are
#' stacked into ONE figure, each with its own bold section header row.
#'
#' @section Statistical validity of the p-value column:
#' Comparing a subgroup's AUC to the AUC of the *entire* sample (which
#' necessarily contains that subgroup) violates the independence assumption
#' underlying both DeLong's test and the naive paired/unpaired bootstrap
#' comparison implemented by \code{pROC::roc.test()}: the two ROC curves
#' are computed on overlapping, non-independent data. The resulting p-value
#' is not statistically valid and tends to be anti-conservative (spuriously
#' significant).
#'
#' To keep the p-value meaningful, this function instead compares each
#' subgroup's AUC against the AUC of the complementary set of subjects
#' (\emph{everyone not in that subgroup}) when \code{compare_method} is
#' \code{"delong_vs_rest"} or \code{"bootstrap_vs_rest"}. These two sets are
#' disjoint, so the independence assumption holds and the test answers a
#' well-posed question: "does the model discriminate differently in this
#' subgroup than in the rest of the sample?"
#'
#' By default (\code{compare_method = "none"}), no p-value is computed and
#' the column is left blank.
#'
#' @param model_obj Fitted model object (\code{Train_Model} or caret \code{train}).
#' @param clinical_data Data frame with clinical/subgroup variables.
#' @param subgroup_var Character vector defining subgroup column names.
#' @param newdata Data frame used for model prediction.
#' @param model_name Optional model name within \code{Train_Model}.
#' @param var_labels Optional named character vector to relabel section headers.
#' @param level_order Optional named list specifying the display order of factor levels.
#' @param indent Character prefix used to indent subgroup levels. Default four spaces.
#' @param dataset_type Dataset origin label (\code{"testing"}, \code{"training"}, or \code{"external"}).
#' @param ref_line Reference vertical line position on graphic (default 0.5).
#' @param xlim Numeric length-2 vector for x-axis limits of forest column.
#' @param ticks_at Axis tick positions for forest plot.
#' @param min_n Minimum subgroup size required for evaluation. Default 5.
#' @param digits Number of decimal places for AUC/CI display.
#' @param p_digits Number of decimal places for p-value display.
#' @param box_col Color for point estimate box in forest graphic.
#' @param ci_col Color for CI whiskers.
#' @param overall_row Logical; whether to prepend an "Overall" dataset row. Default \code{TRUE}.
#' @param overall_label Character label for the overall row.
#' @param compare_method Statistical test for subgroup comparison (\code{"none"}, \code{"delong_vs_rest"}, or \code{"bootstrap_vs_rest"}).
#' @param n_boot Number of bootstrap replicates for bootstrap comparison method.
#' @param theme_args List of visual theme overrides for \code{forestploter::forest_theme()}.
#' @param arrow_labels Length-2 character vector for axis arrows, or \code{NULL}.
#' @param save_plot Logical; whether to save plot as a PDF file.
#' @param save_dir Directory path where plot PDF will be saved.
#' @param file_name Output file name for saved plot.
#' @param width Saved plot width in inches.
#' @param height Saved plot height in inches.
#'
#' @return An invisible list containing the \code{forestploter} graphic and summary \code{data.frame}.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("forestploter", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   set.seed(123)
#'   idx <- sample(1:nrow(mtcars), 20)
#'   model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#'   trained <- ModelTrainAnalysis(model, methods = c("glm"), 
#'   control = list(method = "cv", number = 3), save_plots = FALSE)
#'   clin <- data.frame(row.names = rownames(mtcars), cyl = mtcars$cyl)
#'   PlotSubgroupForest(trained, clinical_data = clin, 
#'   subgroup_var = "cyl", newdata = mtcars, save_plot = FALSE)
#' }
#' }
PlotSubgroupForest <- function(model_obj,
                               clinical_data  = NULL,
                               subgroup_var,
                               newdata        = NULL,
                               model_name     = NULL,
                               var_labels     = NULL,
                               level_order    = NULL,
                               indent         = "    ",
                               ref_line       = 0.5,
                               dataset_type   = c("testing", "training", "external"),
                               xlim           = c(0.4, 1.0),
                               ticks_at       = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                               min_n          = 5,
                               digits         = 3,
                               p_digits       = 3,
                               box_col        = "#377eb8",
                               ci_col         = "black",
                               overall_row    = TRUE,
                               overall_label  = "Overall",
                               compare_method = c("none", "delong_vs_rest", "bootstrap_vs_rest"),
                               n_boot         = 2000,
                               theme_args     = list(),
                               arrow_labels   = NULL,
                               save_plot      = FALSE,
                               save_dir       = NULL,
                               file_name      = NULL,
                               width          = 8,
                               height         = 6) {
  
  .check_clinical_pkgs()
  if (length(subgroup_var) == 0) stop("subgroup_var cannot be empty.")
  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required for this style of plot.\n",
         "Install it with: install.packages('forestploter')")
  }
  
  dataset_type <- match.arg(dataset_type)
  aligned <- .align_clinical_and_newdata(model_obj, clinical_data, newdata, dataset_type)
  clinical_data  <- aligned$clinical_data
  newdata        <- aligned$newdata
  compare_method <- match.arg(compare_method)
  
  probs    <- .predict_probs(model_obj, newdata, model_name)
  true     <- factor(newdata[[model_obj@group_col]])
  positive <- levels(true)[2]
  
  # Safeguarded overall ROC computation
  full_roc <- tryCatch(
    pROC::roc(true, probs, levels = c(levels(true)[1], positive), direction = "auto", quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(full_roc)) {
    stop("Failed to calculate overall ROC curve. Check if outcome labels contain both classes.")
  }
  
  # Helper: Safely calculate ROC object for a subset of samples
  # Prevents crash when sample size < min_n OR outcome contains only one class level
  .roc_for <- function(idx) {
    if (length(idx) < min_n) return(NULL)
    sub_true <- true[idx]
    
    # Check if subgroup contains at least two distinct outcome levels
    if (length(unique(sub_true[!is.na(sub_true)])) < 2) return(NULL)
    
    tryCatch(
      pROC::roc(sub_true, probs[idx],
                levels = c(levels(true)[1], positive),
                direction = "auto", quiet = TRUE),
      error = function(e) NULL
    )
  }
  
  # Helper: Calculate comparison p-value against the complementary "rest" set
  .subgroup_pval <- function(sub_idx) {
    if (compare_method == "none") return(NA_real_)
    
    rest_idx <- setdiff(seq_along(true), sub_idx)
    if (length(rest_idx) < min_n) return(NA_real_)
    
    roc_sub  <- .roc_for(sub_idx)
    roc_rest <- .roc_for(rest_idx)
    
    # If either subgroup or rest set lacks two class levels, return NA
    if (is.null(roc_sub) || is.null(roc_rest)) return(NA_real_)
    
    if (compare_method == "delong_vs_rest") {
      test <- tryCatch(pROC::roc.test(roc_sub, roc_rest, method = "delong"),
                       error = function(e) NULL)
    } else {
      test <- tryCatch(
        pROC::roc.test(roc_sub, roc_rest, method = "bootstrap", boot.n = n_boot),
        error = function(e) NULL
      )
    }
    if (is.null(test)) NA_real_ else test$p.value
  }
  
  # Helper: Generate data frame row for an individual subgroup level
  .auc_row <- function(label, idx) {
    if (length(idx) < min_n) return(NULL)
    roc_obj <- .roc_for(idx)
    if (is.null(roc_obj)) return(NULL) # Skip single-class subgroups safely
    
    ci <- tryCatch(pROC::ci.auc(roc_obj), error = function(e) NULL)
    if (is.null(ci)) return(NULL)
    
    data.frame(Label   = label,
               N       = length(idx),
               AUC     = as.numeric(roc_obj$auc),
               Lower   = ci[1],
               Upper   = ci[3],
               p       = .subgroup_pval(idx),
               Header  = FALSE,
               stringsAsFactors = FALSE)
  }
  
  # Helper: Generate section header row
  .header_row <- function(label, n_total) {
    data.frame(Label = label, N = n_total, AUC = NA_real_, Lower = NA_real_,
               Upper = NA_real_, p = NA_real_, Header = TRUE,
               stringsAsFactors = FALSE)
  }
  
  rows <- list()
  
  # Prepend overall row if requested
  if (overall_row) {
    ci_all <- tryCatch(pROC::ci.auc(full_roc), error = function(e) c(NA, NA, NA))
    rows[[length(rows) + 1]] <- data.frame(
      Label = overall_label, N = length(true), AUC = as.numeric(full_roc$auc),
      Lower = ci_all[1], Upper = ci_all[3], p = NA_real_, Header = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  # Loop over defined subgroup variables
  for (var in subgroup_var) {
    col <- clinical_data[[var]]
    if (is.factor(col)) {
      lvls <- levels(col)
      lvls <- lvls[lvls %in% as.character(col)]
    } else {
      lvls <- sort(unique(as.character(col[!is.na(col)])))
    }
    if (!is.null(level_order) && !is.null(level_order[[var]])) {
      lvls <- level_order[[var]][level_order[[var]] %in% lvls]
    }
    
    header_label <- if (!is.null(var_labels) && !is.na(var_labels[var])) var_labels[var] else var
    rows[[length(rows) + 1]] <- .header_row(header_label, sum(!is.na(col)))
    
    for (lv in lvls) {
      idx <- which(as.character(col) == lv)
      r <- .auc_row(paste0(indent, lv), idx)
      if (!is.null(r)) rows[[length(rows) + 1]] <- r
    }
  }
  
  forest_df <- do.call(rbind, rows)
  if (is.null(forest_df) || nrow(forest_df) == 0) {
    stop("No valid rows to plot - check `subgroup_var`, `min_n`, and class balance across subgroups.")
  }
  
  # Format numeric output strings for forestploter table display
  fmt <- paste0("%.", digits, "f")
  forest_df[["AUC (95% CI)"]] <- ifelse(
    forest_df$Header | is.na(forest_df$AUC), "",
    sprintf(paste0(fmt, " (", fmt, ", ", fmt, ")"),
            forest_df$AUC, forest_df$Lower, forest_df$Upper)
  )
  forest_df[["p value"]] <- ifelse(
    forest_df$Header | is.na(forest_df$p), "",
    formatC(forest_df$p, digits = p_digits, format = "f")
  )
  forest_df[[" "]] <- paste(rep(" ", 20), collapse = " ")
  
  display_df <- forest_df[, c("Label", "N", " ", "AUC (95% CI)", "p value")]
  names(display_df)[1] <- " "
  
  is_summary <- forest_df$Header
  
  # Configure forestploter theme settings
  base_theme <- list(
    base_size    = 10,
    ci_pch       = 15,
    ci_col       = ci_col,
    ci_fill      = box_col,
    ci_lty       = 1,
    ci_lwd       = 1.5,
    ci_Theight   = 0.2,
    refline_lwd  = 1,
    refline_lty  = "dashed",
    refline_col  = "grey40",
    footnote_col = "grey50"
  )
  if (!is.null(arrow_labels)) {
    base_theme$arrow_lab   <- arrow_labels
    base_theme$arrow_type  <- "open"
    base_theme$arrow_col   <- "black"
  }
  theme_args_final <- utils::modifyList(base_theme, theme_args)
  tm <- do.call(forestploter::forest_theme, theme_args_final)
  
  # Render publication-style forest plot
  p <- forestploter::forest(
    display_df,
    est        = forest_df$AUC,
    lower      = forest_df$Lower,
    upper      = forest_df$Upper,
    ci_column  = 3,
    is_summary = is_summary,
    ref_line   = ref_line,
    xlim       = xlim,
    ticks_at   = ticks_at,
    theme      = tm
  )
  
  print(p)
  
  # Optionally save plot to PDF
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- getwd()
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    if (is.null(file_name))
      file_name <- paste0("subgroup_", paste(subgroup_var, collapse = "_"), ".pdf")
    grDevices::pdf(file.path(save_dir, file_name), width = width, height = height)
    plot(p)
    grDevices::dev.off()
  }
  
  return(invisible(list(plot = p, data = forest_df)))
}
#' Plot a Confounder Adjustment Forest Plot
#'
#' Fits a logistic regression of the binary outcome on the model's predicted
#' logit plus a set of clinical covariates, and plots -log10(p-value) for each
#' term as a horizontal bar chart. This is used to check whether the model's
#' prediction remains significantly associated with the outcome after
#' adjusting for potential confounders (e.g. age, sex, stage).
#'
#' For \code{dataset_type = "training"} or \code{"testing"}, the feature data
#' is taken from \code{model_obj@filtered.set} (falling back to
#' \code{model_obj@split.data}), and clinical data is taken automatically from
#' \code{model_obj@process.info$clinical_data} (attached once via
#' \code{AttachClinicalData()} for the full dataset) by intersecting rownames
#' with the feature data — no manual \code{clinical_data} argument is needed
#' in this case.
#'
#' For \code{dataset_type = "external"}, feature data is taken from
#' \code{model_obj@process.info$external_data} and clinical data from
#' \code{model_obj@process.info$external_clinical}, since external samples
#' have no rowname lineage to the training data and must be attached
#' separately via \code{AttachClinicalData(object, external_data = ...,
#' external_clinical = ...)}.
#'
#' @param model_obj A \code{Train_Model} S4 object.
#' @param dataset_type One of \code{"testing"}, \code{"training"}, or
#'   \code{"external"}. Determines which feature/clinical data are used.
#' @param clinical_data Optional \code{data.frame} to override the clinical
#'   data that would otherwise be pulled automatically from \code{model_obj}.
#'   Rownames must match the feature data's rownames. Leave \code{NULL} for
#'   the default (recommended) behavior described above.
#' @param outcome_var Character. Name of the binary outcome column, expected
#'   to be present in either the feature data or the clinical data.
#' @param adjust_vars Optional character vector of clinical column names to
#'   adjust for. If \code{NULL}, all clinical columns (except
#'   \code{outcome_var}) are used as covariates.
#' @param model_name Optional character. Name of a specific trained model in
#'   \code{model_obj@train.models} to use for prediction, or \code{"ensemble"}
#'   to use the ensemble predictor. Defaults to the best model.
#' @param positive_class The outcome level treated as the "positive" class.
#'   Used both to select the corresponding column of predicted probabilities
#'   (passed to \code{.predict_probs()}) and to encode \code{outcome_bin}
#'   (\code{outcome == positive_class}), so the two are guaranteed to refer to
#'   the same class. Defaults to \code{"1"}. Coerced to character for
#'   comparison against outcome levels. If the model's predicted-probability
#'   matrix does not have a column matching \code{positive_class} exactly
#'   (common when outcome labels are numeric-like, since caret/many
#'   classifiers apply \code{make.names()} when training, e.g. turning
#'   \code{"1"} into \code{"X1"}), a \code{make.names()}-adjusted variant is
#'   tried automatically before falling back to \code{.predict_probs()}'s
#'   default column choice.
#' @param save_plot Logical. If \code{TRUE}, saves the plot as a PDF.
#' @param save_dir Directory to save the plot to when \code{save_plot = TRUE}.
#'   Defaults to the current working directory.
#'
#' @return A \code{ggplot} object (invisibly printed as a side effect).
#'
#' @export
PlotConfounderForest <- function (model_obj, dataset_type = c("testing", "training",
                                                              "external"), clinical_data = NULL, outcome_var, adjust_vars = NULL,
                                  model_name = NULL, positive_class = "1", save_plot = FALSE, save_dir = NULL)
{
  #.check_clinical_pkgs()
  dataset_type <- match.arg(dataset_type)
  if (missing(outcome_var) || is.null(outcome_var)) {
    stop("`outcome_var` is required: please specify the name of the outcome column ",
         "(e.g. outcome_var = \"group\").")
  }
  if (dataset_type == "testing") {
    newdata <- model_obj@filtered.set$testing
    if (is.null(newdata))
      newdata <- model_obj@split.data$testing
  }
  else if (dataset_type == "training") {
    newdata <- model_obj@filtered.set$training
    if (is.null(newdata))
      newdata <- model_obj@split.data$training
  }
  else if (dataset_type == "external") {
    newdata <- model_obj@process.info$external_data
  }
  if (is.null(newdata)) {
    if (dataset_type == "external") {
      stop("No external feature data found (process.info$external_data is NULL). ",
           "Please attach external data first (e.g. via AttachExternalData()).")
    }
    else {
      stop("No feature/outcome data found for dataset_type = '",
           dataset_type, "'. Please make sure the model has been split ",
           "(split.data / filtered.set is populated).")
    }
  }
  if (!is.null(clinical_data)) {
    if (!is.data.frame(clinical_data)) {
      stop("`clinical_data` must be a data.frame with rownames matching the feature data. ",
           "Got an object of class '", paste(class(clinical_data), collapse = "/"), "' instead. ",
           "(Did you mean to pass this value to `dataset_type` instead?)")
    }
    message("Using user-provided clinical_data (override). Ensure rownames match the feature data.")
    clinical <- clinical_data
  }
  else {
    if (dataset_type == "external") {
      clinical <- model_obj@process.info$external_clinical
    }
    else {
      clinical <- model_obj@process.info$clinical_data
    }
    if (is.null(clinical)) {
      if (dataset_type == "external") {
        stop("No external clinical data found (process.info$external_clinical is NULL). ",
             "External clinical data is not derived automatically from the training ",
             "data, so please provide it via the `clinical_data` argument, or attach it ",
             "beforehand (e.g. via AttachExternalData(..., clinical_data = ...)).")
      }
      else {
        stop("No internal clinical_data found for dataset_type = '", dataset_type, "'. ",
             "clinical_data is attached once for the full dataset (training + testing ",
             "share the same rownames lineage as data.df/clean.df), so please call ",
             "AttachClinicalData() on this model object before plotting, or pass ",
             "clinical_data manually for a one-off override.")
      }
    }
  }
  common <- intersect(rownames(newdata), rownames(clinical))
  if (length(common) == 0) {
    if (dataset_type == "external") {
      stop("Zero overlapping samples between external feature data and external clinical ",
           "data. Check that rownames of `clinical_data`/external_clinical match the ",
           "rownames of the external feature data exactly.")
    }
    else {
      stop("Zero overlapping samples between feature data and clinical data. Check rownames.")
    }
  }
  if (length(common) < nrow(newdata)) {
    warning("Subsetting newdata to ", length(common), " samples (out of ",
            nrow(newdata), ") that have matching clinical data; the remaining ",
            nrow(newdata) - length(common), " sample(s) were dropped.")
  }
  # Always re-index both by `common` (even when counts already match) so that
  # rows are guaranteed to be aligned in the same order before cbind() below,
  # which binds by position, not by rowname.
  newdata <- newdata[common, , drop = FALSE]
  clinical <- clinical[common, , drop = FALSE]
  if (!outcome_var %in% colnames(newdata) && outcome_var %in%
      colnames(clinical)) {
    newdata[[outcome_var]] <- clinical[[outcome_var]]
    message("Outcome column '", outcome_var, "' was copied from clinical_data to newdata.")
  }
  if (!outcome_var %in% colnames(newdata)) {
    stop("Outcome column '", outcome_var, "' not found in newdata or clinical_data.")
  }
  outcome_levels <- levels(factor(newdata[[outcome_var]]))
  positive_class <- as.character(positive_class)
  if (!positive_class %in% outcome_levels) {
    stop("`positive_class` = '", positive_class, "' is not among the observed levels of '",
         outcome_var, "' (", paste(outcome_levels, collapse = ", "), "). ",
         "Please set `positive_class` to one of these values.")
  }
  all_clin_vars <- colnames(clinical)
  if (is.null(adjust_vars)) {
    covars <- setdiff(all_clin_vars, outcome_var)
  }
  else {
    if (!is.character(adjust_vars)) {
      stop("adjust_vars must be a character vector or NULL.")
    }
    covars <- setdiff(adjust_vars, outcome_var)
    missing_vars <- setdiff(covars, all_clin_vars)
    if (length(missing_vars) > 0) {
      warning("The following variables in adjust_vars are not found in clinical_data and will be ignored: ",
              paste(missing_vars, collapse = ", "))
      covars <- intersect(covars, all_clin_vars)
    }
    if (length(covars) == 0) {
      message("No valid covariates remain after filtering. Model will contain only the prediction logit.")
    }
  }
  # `.predict_probs()` matches `positive_class` against the column names of
  # the model's predicted-probability matrix. Those column names are often
  # NOT the raw outcome labels: if the outcome is numeric-like (e.g. "0"/"1"),
  # caret/most classifiers apply make.names() when training, turning "1" into
  # "X1". Try the raw label first, then a make.names()-adjusted fallback,
  # instead of silently trusting `.predict_probs()`'s own "use 2nd column"
  # default (which is only correct by coincidence of column ordering).
  .try_predict_probs <- function(pc) {
    mismatched <- FALSE
    result <- withCallingHandlers(
      .predict_probs(model_obj, newdata, model_name, positive_class = pc),
      warning = function(w) {
        if (grepl("not found", conditionMessage(w), fixed = TRUE)) {
          mismatched <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    )
    list(probs = result, mismatched = mismatched)
  }
  attempt <- .try_predict_probs(positive_class)
  if (attempt$mismatched) {
    alt_class <- make.names(positive_class)
    if (!identical(alt_class, positive_class)) {
      attempt <- .try_predict_probs(alt_class)
    }
    if (attempt$mismatched) {
      warning("Could not find a probability column matching positive_class = '", positive_class,
              "' (also tried '", alt_class, "'). Falling back to the second column as chosen by ",
              ".predict_probs(); please verify this is really the class you intend, or pass a ",
              "`positive_class` value that matches the model's actual probability column names.")
    }
  }
  probs <- attempt$probs
  eps <- 1e-07
  probs_safe <- pmin(pmax(probs, eps), 1 - eps)
  logit_probs <- log(probs_safe/(1 - probs_safe))
  df <- cbind(clinical, prediction_logit = logit_probs)
  df$outcome <- factor(newdata[[outcome_var]])
  df$outcome_bin <- as.numeric(as.character(df$outcome) == positive_class)
  if (length(covars) == 0) {
    frm <- as.formula("outcome_bin ~ prediction_logit")
  }
  else {
    frm <- as.formula(paste("outcome_bin ~ prediction_logit +",
                            paste(covars, collapse = " + ")))
  }
  fit <- glm(frm, data = df, family = binomial)
  tem <- summary(fit)$coefficients
  tem <- as.data.frame(tem)
  colnames(tem) <- c("Estimate", "StdError", "zValue", "pVal")
  tem <- tem[rownames(tem) != "(Intercept)", , drop = FALSE]
  tem$pVal_safe <- pmax(tem$pVal, 1e-15)
  tem$`-log10P` <- -log10(tem$pVal_safe)
  tem$Variable <- rownames(tem)
  tem$Variable[tem$Variable == "prediction_logit"] <- "Model Prediction (Logit)"
  tem$Variable <- factor(tem$Variable, levels = tem$Variable[order(tem$`-log10P`)])
  p <- ggplot2::ggplot(tem, ggplot2::aes(x = .data[["-log10P"]],
                                         y = .data[["Variable"]], fill = .data[["-log10P"]])) +
    ggplot2::geom_col() + ggplot2::scale_fill_gradientn(colours = c("#f9ddda",
                                                                    "#eda8bd", "#ce78b3", "#9955a8", "#573b88"), name = expression(-log[10]("P"))) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                        color = "black", linewidth = 1) + ggplot2::annotate("text",
                                                                            x = -log10(0.05) + 0.02, y = 1, label = "p = 0.05", hjust = 0,
                                                                            size = 3.5, color = "black") + ggplot2::labs(title = paste0("Confounder Adjustment (",
                                                                                                                                        dataset_type, ")"), x = expression(-log[10](p - value)),
                                                                                                                         y = NULL) + ggprism::theme_prism(base_size = 13) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                                                                                                                                                                                              face = "bold"))
  print(p)
  if (save_plot) {
    if (is.null(save_dir))
      save_dir <- getwd()
    if (!dir.exists(save_dir))
      dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("confounder_bar_",
                                               dataset_type, ".pdf")), plot = p, width = 7, height = 5,
                    dpi = 300)
  }
  return(p)
}

#' Calculate Multi-Threshold Metrics with Custom Targets
#'
#' Computes a full table of Accuracy, PPV, NPV across all unique thresholds,
#' then finds:
#'   - Youden index (Se + Sp - 1)
#'   - Thresholds that reach a target Sensitivity, Specificity, PPV, or NPV
#'   - Threshold that maximises Accuracy (if requested)
#'
#' @section A note on where to call this:
#' This function \strong{searches} for a threshold (Youden index, "reach
#' target Se/Sp/PPV/NPV", "maximise Accuracy") over whatever data you pass
#' in. If you call it on the testing (or external) set and then report that
#' threshold's Accuracy/Sensitivity/etc. as "test-set performance", you are
#' reporting an optimistic, overfit number: the threshold was chosen
#' specifically because it did well on that same data. Prefer determining
#' the threshold on the \strong{training} set (\code{dataset_type =
#' "training"}), then applying that fixed numeric value to the testing/
#' external set yourself (e.g. via \code{ApplyThreshold(..., custom_threshold
#' = your_value)}) for a genuine out-of-sample evaluation.
#'
#' To keep this visible, calling this function with \code{dataset_type =
#' "testing"} or \code{"external"} emits a warning every time. The warning
#' does not block execution -- e.g. you may legitimately want to inspect the
#' full metrics table for exploratory purposes -- but any threshold selected
#' under that warning should not be reported as if it were fixed in advance.
#'
#' @param model_obj,newdata,model_name  As elsewhere.
#' @param dataset_type Character string indicating the dataset origin for plot
#'   labelling. Must be one of \code{"testing"}, \code{"training"}, or
#'   \code{"external"}. Default is \code{"testing"}.
#' @param target_se      Target Sensitivity (e.g., 0.9). Default NULL.
#' @param target_sp      Target Specificity (e.g., 0.9). Default NULL.
#' @param target_ppv     Target PPV (e.g., 0.95). Default NULL.
#' @param target_npv     Target NPV (e.g., 0.95). Default NULL.
#' @param target_acc     Logical; if TRUE, find threshold maximising Accuracy.
#' @return A list with \code{thresholds}, \code{metrics_df}, \code{probabilities},
#'   \code{true}, \code{positive}, \code{negative}.
#' @export
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' set.seed(123)
#' idx <- sample(1:nrow(mtcars), 20)
#' model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#' trained <- ModelTrainAnalysis(model, methods = c("glm"), control = list(method = "cv", number = 3),
#' save_plots = FALSE)
#' # Calculate thresholds on training set (preferred)
#' thresh <- CalculateThresholds(trained, newdata = trained@filtered.set$training, 
#' dataset_type = "training")
#' print(thresh$thresholds)
#' # Apply to test set
#' apply_res <- ApplyThreshold(thresh, which_threshold = "Youden", 
#' newdata = trained@filtered.set$testing)
#' print(apply_res$metrics)
#' # ClinicalThreshold wrapper (plots if interactive)
#' if (interactive()) {
#'   ClinicalThreshold(model_obj = trained, newdata = trained@filtered.set$testing, 
#'   save_plot = FALSE)
#' }
#' }
CalculateThresholds <- function(model_obj,
                                newdata      = NULL,
                                model_name   = NULL,
                                dataset_type = c("testing", "training", "external"),
                                target_se    = NULL,
                                target_sp    = NULL,
                                target_ppv   = NULL,
                                target_npv   = NULL,
                                target_acc   = TRUE) {
  .check_clinical_pkgs()
  dataset_type <- match.arg(dataset_type)
  
  if (dataset_type %in% c("testing", "external")) {
    warning(
      "CalculateThresholds() is searching for a threshold (Youden / target ",
      "Se-Sp-PPV-NPV / max Accuracy) directly on the '", dataset_type, "' ",
      "set. Metrics reported at that threshold on this same data are ",
      "optimistically biased (threshold overfitting) and should not be ",
      "reported as generalizable performance. Prefer determining the ",
      "threshold on the training set, then applying that fixed value to ",
      "this set (e.g. via ApplyThreshold(..., custom_threshold = ...)).",
      call. = FALSE
    )
  }
  
  if (is.null(newdata) && inherits(model_obj, "Train_Model")) {
    aligned <- .align_clinical_and_newdata(model_obj, newdata = newdata, dataset_type = dataset_type)
    newdata <- aligned$newdata
  }
  stopifnot(!is.null(newdata))
  
  probs <- .predict_probs(model_obj, newdata, model_name)
  true  <- factor(newdata[[model_obj@group_col]])
  positive <- levels(true)[2]
  negative <- levels(true)[1]
  
  roc_obj <- pROC::roc(true, probs, levels = c(negative, positive),
                       direction = "auto", quiet = TRUE)
  coords_all <- pROC::coords(roc_obj, "all", ret = c("threshold", "se", "sp"))
  
  youden <- coords_all[which.max(coords_all$se + coords_all$sp - 1), ]
  
  uniq_thr <- sort(unique(round(probs, 4)), decreasing = TRUE)
  met_list <- lapply(uniq_thr, function(t) {
    pred_class <- factor(ifelse(probs > t, positive, negative),
                         levels = c(negative, positive))
    cm <- caret::confusionMatrix(pred_class, true, positive = positive)
    data.frame(Threshold = t,
               Sensitivity = cm$byClass["Sensitivity"],
               Specificity = cm$byClass["Specificity"],
               Accuracy    = cm$overall["Accuracy"],
               PPV         = cm$byClass["Pos Pred Value"],
               NPV         = cm$byClass["Neg Pred Value"],
               F1          = cm$byClass["F1"],
               Precision   = cm$byClass["Precision"],
               Recall      = cm$byClass["Recall"])
  })
  metrics_df <- do.call(rbind, met_list)
  
  .find_thresh <- function(metric_col, target, larger = FALSE) {
    if (larger) {
      idx <- which(metrics_df[[metric_col]] >= target)
      if (length(idx) == 0) return(NULL)
      metrics_df[idx[which.min(metrics_df$Threshold[idx])], ]
    } else {
      dist <- abs(metrics_df[[metric_col]] - target)
      metrics_df[which.min(dist), ]
    }
  }
  
  thresholds <- c(Youden = youden$threshold)
  if (!is.null(target_se)) {
    se_row <- .find_thresh("Sensitivity", target_se, larger = TRUE)
    if (!is.null(se_row)) thresholds <- c(thresholds, Se_Target = se_row$Threshold)
  }
  if (!is.null(target_sp)) {
    sp_row <- .find_thresh("Specificity", target_sp, larger = TRUE)
    if (!is.null(sp_row)) thresholds <- c(thresholds, Sp_Target = sp_row$Threshold)
  }
  if (!is.null(target_ppv)) {
    ppv_row <- .find_thresh("PPV", target_ppv, larger = TRUE)
    if (!is.null(ppv_row)) thresholds <- c(thresholds, PPV_Target = ppv_row$Threshold)
  }
  if (!is.null(target_npv)) {
    npv_row <- .find_thresh("NPV", target_npv, larger = TRUE)
    if (!is.null(npv_row)) thresholds <- c(thresholds, NPV_Target = npv_row$Threshold)
  }
  if (target_acc) {
    acc_max <- metrics_df[which.max(metrics_df$Accuracy), ]
    thresholds <- c(thresholds, MaxAcc = acc_max$Threshold)
  }
  
  list(thresholds    = thresholds,
       metrics_df    = metrics_df,
       probabilities = probs,
       true          = true,
       positive      = positive,
       negative      = negative)
}
#' Calculate Multi-Threshold Metrics from External Probabilities
#'
#' Accepts predicted probabilities and true labels directly, without requiring a
#' model object. This is useful when predictions come from an external source,
#' such as literature or other software. The function computes classification
#' metrics at each unique threshold, including Sensitivity, Specificity,
#' Accuracy, PPV, NPV, and F1, and returns selected thresholds based on
#' user-defined targets.
#'
#' @param probs Numeric vector of predicted probabilities (range 0-1).
#' @param true Factor vector of true binary labels.
#' @param positive Character string specifying the positive class
#'   (e.g., \code{"yes"}). Must be one of the levels of \code{true}.
#' @param target_se Numeric. Target Sensitivity threshold. If specified, the
#'   function returns the threshold that achieves at least this Sensitivity
#'   (closest lower threshold). Default \code{NULL}.
#' @param target_sp Numeric. Target Specificity threshold. If specified, returns
#'   the threshold that achieves at least this Specificity. Default \code{NULL}.
#' @param target_ppv Numeric. Target Positive Predictive Value threshold.
#'   If specified, returns the threshold that achieves at least this PPV.
#'   Default \code{NULL}.
#' @param target_npv Numeric. Target Negative Predictive Value threshold.
#'   If specified, returns the threshold that achieves at least this NPV.
#'   Default \code{NULL}.
#' @param target_acc Logical. If \code{TRUE}, returns the threshold that
#'   maximises Accuracy. Default \code{TRUE}.
#'
#' @export
#'
#' @return A list with the same structure as \code{\link{CalculateThresholds}},
#'   containing:
#'   \itemize{
#'     \item \code{thresholds}: Named numeric vector of selected thresholds
#'       (Youden, targets, MaxAcc).
#'     \item \code{metrics_df}: Data frame of all metrics at each unique threshold.
#'     \item \code{probabilities}: Input probabilities.
#'     \item \code{true}: Input true labels.
#'     \item \code{positive}: Positive class name.
#'     \item \code{negative}: Negative class name.
#'   }
#' @examples
#' \dontrun{
#' library(caret)
#' 
#' # Simulate example data
#' set.seed(123)
#' true_labels <- factor(sample(c("yes", "no"), 100, replace = TRUE),
#'                       levels = c("no", "yes"))
#' predicted_probs <- runif(100)
#' 
#' # Basic usage with Youden and Max Accuracy thresholds
#' result <- CalculateThresholdsFromProbs(
#'   probs = predicted_probs,
#'   true = true_labels,
#'   positive = "yes"
#' )
#' 
#' # View selected thresholds
#' print(result$thresholds)
#' 
#' # View first few rows of metrics table
#' head(result$metrics_df)
#' 
#' # With target specificity and sensitivity
#' result2 <- CalculateThresholdsFromProbs(
#'   probs = predicted_probs,
#'   true = true_labels,
#'   positive = "yes",
#'   target_se = 0.9,
#'   target_sp = 0.8,
#'   target_ppv = 0.85
#' )
#' 
#' print(result2$thresholds)
#' }
CalculateThresholdsFromProbs <- function(probs,
                                         true,
                                         positive,
                                         target_se  = NULL,
                                         target_sp  = NULL,
                                         target_ppv = NULL,
                                         target_npv = NULL,
                                         target_acc = TRUE) {
  .check_clinical_pkgs()
  if (!is.factor(true)) true <- factor(true)
  if (!positive %in% levels(true))
    stop("'positive' must be one of the levels of 'true'.")
  negative <- setdiff(levels(true), positive)[1]
  
  roc_obj <- pROC::roc(true, probs, levels = c(negative, positive),
                       direction = "auto", quiet = TRUE)
  coords_all <- pROC::coords(roc_obj, "all", ret = c("threshold", "se", "sp"))
  
  youden <- coords_all[which.max(coords_all$se + coords_all$sp - 1), ]
  
  uniq_thr <- sort(unique(round(probs, 4)), decreasing = TRUE)
  met_list <- lapply(uniq_thr, function(t) {
    pred_class <- factor(ifelse(probs > t, positive, negative),
                         levels = c(negative, positive))
    cm <- caret::confusionMatrix(pred_class, true, positive = positive)
    data.frame(Threshold = t,
               Sensitivity = cm$byClass["Sensitivity"],
               Specificity = cm$byClass["Specificity"],
               Accuracy    = cm$overall["Accuracy"],
               PPV         = cm$byClass["Pos Pred Value"],
               NPV         = cm$byClass["Neg Pred Value"],
               F1          = cm$byClass["F1"],
               Precision   = cm$byClass["Precision"],
               Recall      = cm$byClass["Recall"])
  })
  metrics_df <- do.call(rbind, met_list)
  
  .find_thresh <- function(metric_col, target, larger = FALSE) {
    if (larger) {
      idx <- which(metrics_df[[metric_col]] >= target)
      if (length(idx) == 0) return(NULL)
      metrics_df[idx[which.min(metrics_df$Threshold[idx])], ]
    } else {
      dist <- abs(metrics_df[[metric_col]] - target)
      metrics_df[which.min(dist), ]
    }
  }
  
  thresholds <- c(Youden = youden$threshold)
  if (!is.null(target_se)) {
    se_row <- .find_thresh("Sensitivity", target_se, larger = TRUE)
    if (!is.null(se_row)) thresholds <- c(thresholds, Se_Target = se_row$Threshold)
  }
  if (!is.null(target_sp)) {
    sp_row <- .find_thresh("Specificity", target_sp, larger = TRUE)
    if (!is.null(sp_row)) thresholds <- c(thresholds, Sp_Target = sp_row$Threshold)
  }
  if (!is.null(target_ppv)) {
    ppv_row <- .find_thresh("PPV", target_ppv, larger = TRUE)
    if (!is.null(ppv_row)) thresholds <- c(thresholds, PPV_Target = ppv_row$Threshold)
  }
  if (!is.null(target_npv)) {
    npv_row <- .find_thresh("NPV", target_npv, larger = TRUE)
    if (!is.null(npv_row)) thresholds <- c(thresholds, NPV_Target = npv_row$Threshold)
  }
  if (target_acc) {
    acc_max <- metrics_df[which.max(metrics_df$Accuracy), ]
    thresholds <- c(thresholds, MaxAcc = acc_max$Threshold)
  }
  
  list(thresholds    = thresholds,
       metrics_df    = metrics_df,
       probabilities = probs,
       true          = true,
       positive      = positive,
       negative      = negative)
}

#' Accuracy/PPV/NPV vs Threshold Plot with Custom Threshold Markers
#'
#' @param thresh_result Output from \code{CalculateThresholds}.
#' @param save_plot,save_dir  Output options.
#' @return A ggplot object.
#' @export
PlotThresholdAccuracy <- function(thresh_result,
                                  save_plot = FALSE,
                                  save_dir  = NULL) {
  df <- thresh_result$metrics_df
  df_long <- tidyr::pivot_longer(df, cols = c("Accuracy", "PPV", "NPV"),
                                 names_to = "Metric", values_to = "Value")
  thr <- thresh_result$thresholds
  
  # One vertical reference line + one label PER THRESHOLD, instead of one
  # point/label per (threshold x metric-curve) combination. The old approach
  # placed 3 overlapping "asterisk + text" labels for every named threshold
  # (its Accuracy value, its PPV value, its NPV value), which becomes an
  # unreadable pile-up whenever two thresholds land close together (e.g.
  # Youden and MaxAcc). Thresholds that round to the same value are merged
  # into a single "Youden / MaxAcc" label so they don't stack either.
  thr_df <- data.frame(Threshold = round(as.numeric(thr), 4), Label = names(thr))
  thr_df <- do.call(rbind, lapply(split(thr_df, thr_df$Threshold), function(g) {
    data.frame(Threshold = g$Threshold[1], Label = paste(g$Label, collapse = " / "))
  }))
  thr_df$y_label <- 1.06
  
  cols <- c("Accuracy" = "#800000", "PPV" = "#767676", "NPV" = "#cc8214")
  lty  <- c("Accuracy" = "solid",  "PPV" = "dotted",   "NPV" = "dashed")
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Threshold, y = Value,
                                             color = Metric, linetype = Metric)) +
    ggplot2::geom_line(size = 1.1) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_linetype_manual(values = lty) +
    ggplot2::geom_vline(data = thr_df, ggplot2::aes(xintercept = Threshold),
                        color = "red", linetype = "dotted", alpha = 0.5, inherit.aes = FALSE) +
    ggrepel::geom_text_repel(data = thr_df,
                             ggplot2::aes(x = Threshold, y = y_label,
                                          label = paste0(Label, "\n(", signif(Threshold, 2), ")")),
                             inherit.aes = FALSE, color = "red", size = 3.2,
                             direction = "x", segment.size = 0.3, min.segment.length = 0,
                             max.overlaps = Inf) +
    ggplot2::scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::labs(title = "Accuracy / PPV / NPV vs Threshold",
                  x = "Threshold", y = "Value") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom")
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "threshold_accuracy.pdf"),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  return(p)
}

#' Decision Density with Threshold Zones
#'
#' Visualizes predicted probability density distributions stratified by true outcome,
#' overlaid with decision threshold lines and risk zone summary statistics (Low/High counts,
#' PPV, and NPV).
#'
#' @param thresh_result Output list from \code{CalculateThresholds}.
#' @param lower_threshold Numeric threshold defining the lower decision boundary.
#'   Defaults to \code{Se_Target} if available, otherwise falls back to \code{Youden}.
#' @param upper_threshold Numeric threshold defining the upper decision boundary.
#'   Defaults to \code{Sp_Target} if available, otherwise falls back to \code{Youden}.
#' @param save_plot Logical; whether to save the rendered density plot as a PDF file. Default is \code{FALSE}.
#' @param save_dir Directory path where the output PDF will be saved.
#'
#' @return A \code{ggplot} object showing probability density distributions and decision zones.
#' @export
PlotThresholdDensity <- function(thresh_result,
                                 lower_threshold = NULL,
                                 upper_threshold = NULL,
                                 save_plot = FALSE,
                                 save_dir  = NULL) {
  # Validate required packages
  .check_clinical_pkgs()
  
  # Extract predictions and true labels
  probs    <- thresh_result$probabilities
  true     <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  
  # Step 1: Safely resolve lower and upper thresholds with fallback options
  if (is.null(lower_threshold) || is.na(lower_threshold)) {
    lower_threshold <- thresh_result$thresholds["Se_Target"]
    if (is.null(lower_threshold) || is.na(lower_threshold)) {
      lower_threshold <- thresh_result$thresholds["Youden"]
    }
  }
  
  if (is.null(upper_threshold) || is.na(upper_threshold)) {
    upper_threshold <- thresh_result$thresholds["Sp_Target"]
    if (is.null(upper_threshold) || is.na(upper_threshold)) {
      upper_threshold <- thresh_result$thresholds["Youden"]
    }
  }
  
  df <- data.frame(prob = probs, group = true)
  total_n <- nrow(df)
  
  # Step 2: Compute counts across decision zones
  low_count  <- sum(df$prob <= lower_threshold)
  mid_count  <- sum(df$prob > lower_threshold & df$prob <= upper_threshold)
  high_count <- total_n - low_count - mid_count
  
  # Step 3: Compute NPV and PPV safely to prevent NaN when subset is empty
  low_npv_str <- if (low_count > 0) {
    sprintf("%.1f%%", mean(df$group[df$prob <= lower_threshold] == negative) * 100)
  } else {
    "N/A"
  }
  
  high_ppv_str <- if (high_count > 0) {
    sprintf("%.1f%%", mean(df$group[df$prob > upper_threshold] == positive) * 100)
  } else {
    "N/A"
  }
  
  cols <- c("#969696", "#fed9a6")
  max_density_y <- max(density(df$prob)$y, na.rm = TRUE)
  
  # Step 4: Render density plot with decision zone annotations
  p <- ggplot2::ggplot(df, ggplot2::aes(x = prob, fill = group)) +
    ggplot2::geom_density(alpha = 0.6, colour = NA) +
    ggplot2::scale_fill_manual(values = cols) +
    # Highlight Low-risk and High-risk decision zones
    ggplot2::annotate("rect", xmin = -Inf, xmax = lower_threshold,
                      ymin = -Inf, ymax = Inf,
                      fill = "#ffffb3", alpha = 0.1) +
    ggplot2::annotate("rect", xmin = upper_threshold, xmax = Inf,
                      ymin = -Inf, ymax = Inf,
                      fill = "#8dd3c7", alpha = 0.1) +
    ggplot2::geom_vline(xintercept = c(lower_threshold, upper_threshold),
                        linetype = "dashed", color = "grey40") +
    # Add Low-risk zone label
    ggplot2::annotate("text", x = lower_threshold / 2,
                      y = max_density_y * 0.9,
                      label = paste0("Low: ", low_count, " (", round(low_count / total_n * 100, 1),
                                     "%)\nNPV: ", low_npv_str), size = 3.5) +
    # Add High-risk zone label
    ggplot2::annotate("text", x = (upper_threshold + 1) / 2,
                      y = max_density_y * 0.9,
                      label = paste0("High: ", high_count, " (", round(high_count / total_n * 100, 1),
                                     "%)\nPPV: ", high_ppv_str), size = 3.5) +
    ggplot2::labs(title = "Prediction Density with Decision Zones",
                  x = "Predicted Probability", y = "Density") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = "top")
  
  # Display plot
  print(p)
  
  # Step 5: Save plot to disk if requested
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- getwd()
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "threshold_density.pdf"),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  
  return(p)
}

#' Waterfall Plot for Threshold Classification
#'
#' @param thresh_result Output from \code{CalculateThresholds}.
#' @param which_threshold Name of the threshold in \code{thresh_result$thresholds}.
#' @param save_plot,save_dir,colors  Output options.
#' @return A ggplot object.
#' @export
PlotThresholdWaterfall <- function(thresh_result,
                                   which_threshold = "Youden",
                                   save_plot = FALSE,
                                   save_dir  = NULL,
                                   colors = c("#f1a340", "#998ec3")) {
  probs <- thresh_result$probabilities
  true  <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  thr <- thresh_result$thresholds[which_threshold]
  
  df <- data.frame(id = seq_along(probs), prob = probs, truth = true)
  df$dif <- df$prob - thr
  df <- df[order(df$dif), ]
  df$predict <- ifelse(df$prob > thr, positive, negative)
  df$correct <- ifelse(df$predict == df$truth, "Correct", "Wrong")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(id, dif), y = dif, fill = correct)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(title = paste("Waterfall Plot --", which_threshold, "Threshold"),
                  x = "Samples (sorted)", y = "Difference from threshold") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x  = ggplot2::element_blank(),
                   plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("waterfall_", which_threshold, ".pdf")),
                    plot = p, width = 8, height = 5, dpi = 300)
  }
  return(p)
}

#' Confusion Matrix Heatmap (Customizable Colors)
#'
#' @param thresh_result Output from \code{CalculateThresholds}.
#' @param which_threshold Name of the threshold.
#' @param save_plot Logical. Save plot?
#' @param save_dir Output directory.
#' @param fill_colors Character vector of colors for gradient.
#'   Default \code{c("#d8b365", "#f5f5f5", "#5ab4ac")}.
#' @return A ggplot object.
#' @export
PlotThresholdConfusion <- function(thresh_result,
                                   which_threshold = "Youden",
                                   save_plot = FALSE,
                                   save_dir = NULL,
                                   fill_colors = c("#d8b365", "#f5f5f5", "#5ab4ac")) {
  probs <- thresh_result$probabilities
  true  <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  thr <- thresh_result$thresholds[which_threshold]
  
  pred_class <- factor(ifelse(probs > thr, positive, negative),
                       levels = c(negative, positive))
  cm <- caret::confusionMatrix(pred_class, true, positive = positive)
  tab <- as.data.frame(cm$table)
  colnames(tab) <- c("Predicted", "Actual", "Freq")
  tab <- tab %>%
    dplyr::group_by(Actual) %>%
    dplyr::mutate(Pct = round(Freq / sum(Freq) * 100, 1)) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = Actual, y = Predicted, fill = Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Freq, "\n(", Pct, "%)")),
                       size = 5, fontface = "bold") +
    ggplot2::scale_fill_gradientn(colours = fill_colors) +
    ggplot2::labs(title = paste("Confusion Matrix --", which_threshold, "Threshold"),
                  x = "Actual", y = "Predicted", fill = "Count") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("confusion_", which_threshold, ".pdf")),
                    plot = p, width = 5, height = 4.5, dpi = 300)
  }
  return(p)
}

#' Compare Threshold-Based Classifiers with Original Score (Final, Fixed)
#'
#' Evaluates and visualizes ROC curves overlaid with specific threshold operating
#' points (e.g., Youden index, custom targets). Optionally compares performance
#' against a second classifier model.
#'
#' @param thresh_result Output list from \code{CalculateThresholds}.
#' @param compare_model Optional second threshold result list for model comparison.
#' @param compare_label Character label for the comparison model. Default is \code{"Clinician"}.
#' @param save_plot Logical; whether to save the rendered ROC plot as a PDF. Default is \code{FALSE}.
#' @param save_dir Directory path where the output PDF will be saved.
#'
#' @return A \code{ggplot} object representing the ROC curve with marked threshold points.
#' @export
PlotThresholdROC <- function(thresh_result,
                             compare_model = NULL,
                             compare_label = "Clinician",
                             save_plot = FALSE,
                             save_dir  = NULL) {
  .require_pkgs(c("nricens", "timeROC"))
  # Validate required packages
  .check_clinical_pkgs()
  
  # Extract prediction probabilities and true binary outcome labels
  probs    <- thresh_result$probabilities
  true     <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  
  # Build primary ROC curve object
  roc_main <- pROC::roc(true, probs, levels = c(negative, positive),
                        direction = "auto", quiet = TRUE)
  auc_main <- round(as.numeric(pROC::auc(roc_main)), 3)
  
  thr_names <- names(thresh_result$thresholds)
  
  # Step 1: Safely extract coordinates (Sensitivity and Specificity) for each threshold
  sens_sp <- lapply(thr_names, function(nm) {
    tval <- thresh_result$thresholds[nm]
    co <- tryCatch({
      # Explicitly specify x = tval and input = "threshold" to prevent pROC from
      # mistaking numeric thresholds for point indices, and remove redundant best.method
      res <- pROC::coords(roc_main, x = tval, input = "threshold",
                          ret = c("sensitivity", "specificity"))
      se  <- if (is.data.frame(res) || is.list(res)) as.numeric(res$sensitivity[1]) else as.numeric(res[1])
      sp  <- if (is.data.frame(res) || is.list(res)) as.numeric(res$specificity[1]) else as.numeric(res[2])
      data.frame(Threshold = nm, Sensitivity = se, Specificity = sp)
    }, error = function(e) NULL)
    
    if (is.null(co) || any(is.na(co[, c("Sensitivity", "Specificity")]))) return(NULL)
    co
  })
  
  thr_df <- do.call(rbind, sens_sp)
  if (is.null(thr_df) || nrow(thr_df) == 0) {
    stop("No valid threshold coordinates could be extracted.")
  }
  
  # Calculate False Positive Rate (FPR = 1 - Specificity)
  thr_df$FPR <- 1 - thr_df$Specificity
  
  # Construct dataframe for main ROC curve
  roc_df <- data.frame(FPR = 1 - roc_main$specificities,
                       TPR = roc_main$sensitivities)
  
  n_thr <- nrow(thr_df)
  cols  <- wesanderson::wes_palette("Darjeeling1", max(4, n_thr + 1), type = "discrete")
  
  # Step 2: Render primary ROC curve and mark threshold operating points
  p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = cols[1], linewidth = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(data = thr_df, size = 4,
                        mapping = ggplot2::aes(x = FPR, y = Sensitivity),
                        color = cols[2:(n_thr + 1)]) +
    ggrepel::geom_text_repel(data = thr_df,
                             mapping = ggplot2::aes(x = FPR, y = Sensitivity, label = Threshold),
                             size = 3.5) +
    ggplot2::annotate("text", x = 0.75, y = 0.25,
                      label = paste0("AUC = ", auc_main), size = 4, color = cols[1]) +
    ggplot2::labs(title = "ROC with Threshold Operating Points",
                  x = "1 - Specificity", y = "Sensitivity") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  # Step 3: Overlay comparison model ROC curve if provided
  if (!is.null(compare_model)) {
    roc_comp <- pROC::roc(compare_model$true, compare_model$probabilities,
                          levels = c(negative, positive), direction = "auto", quiet = TRUE)
    auc_comp <- round(as.numeric(pROC::auc(roc_comp)), 3)
    comp_df  <- data.frame(FPR = 1 - roc_comp$specificities,
                           TPR = roc_comp$sensitivities)
    
    p <- p +
      ggplot2::geom_line(data = comp_df, color = cols[n_thr + 2], linewidth = 1.2, linetype = "dashed") +
      ggplot2::annotate("text", x = 0.75, y = 0.15,
                        label = paste0(compare_label, " AUC = ", auc_comp),
                        size = 4, color = cols[n_thr + 2])
  }
  
  # Display plot
  print(p)
  
  # Step 4: Save plot to disk if requested
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- getwd()
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "threshold_roc.pdf"),
                    plot = p, width = 7, height = 6, dpi = 300)
  }
  
  return(p)
}

#' Calculate NRI and IDI between two survival/mortality prediction models
#'
#' @description
#' Computes the Net Reclassification Improvement (NRI) and Integrated
#' Discrimination Improvement (IDI) for two sets of predicted probabilities,
#' using the \code{nricens} package. Results are displayed as a forest-style
#' plot with 95% confidence intervals, optionally saved to PDF.
#'
#' @param thresh_result1 A list from a threshold evaluation (must contain
#'   \code{probabilities}, \code{true}, and \code{positive} components).
#'   This is considered the "new" model.
#' @param thresh_result2 A similar list for the "standard" or "old" model.
#' @param label1 Character label for the new model (default: "Model 1").
#' @param label2 Character label for the standard model (default: "Model 2").
#' @param cutoffs Numeric vector of probability cutoffs for categorical NRI
#'   (default: \code{c(0.5)}). Passed to \code{nricens::nribin}.
#' @param save_plot Logical; if \code{TRUE}, saves the plot as a PDF.
#' @param save_dir Directory path where the PDF will be saved. Created
#'   recursively if it does not exist. Ignored if \code{save_plot = FALSE}.
#'
#' @return A list (invisibly) with three components:
#' \item{nri}{The full output from \code{nricens::nribin} (contains both NRI and IDI).}
#' \item{idi}{The IDI table extracted from the \code{nribin} result (for convenience).}
#' \item{plot}{The ggplot2 object representing the forest plot.}
#'
#' @details
#' The function uses \code{nricens::nribin} with \code{updown = "category"} to
#' compute categorical NRI (NRI+) and IDI simultaneously. Estimates and standard
#' errors are converted to percentages for plotting. If the estimation fails,
#' a placeholder plot is generated. The plot includes a grey reference region
#' between -5% and +5% to highlight clinically negligible changes.
#'
#' @seealso \code{\link[nricens]{nribin}}
#' @export
#' 
#' @examples
#' \dontrun{
#' # Assuming 'res1' and 'res2' are results from threshold evaluation
#' result <- CalculateNRI(res1, res2, label1 = "New Model", label2 = "Old Model")
#' print(result$plot)
#' }
CalculateNRI <- function(thresh_result1,
                         thresh_result2,
                         label1 = "Model 1",
                         label2 = "Model 2",
                         cutoffs = c(0.5),
                         save_plot = FALSE,
                         save_dir  = NULL) {
  .require_pkgs("nricens")
  # Ensure required packages are available (assumes .check_clinical_pkgs exists)
  .check_clinical_pkgs()
  
  # Extract probabilities and true labels
  probs1 <- thresh_result1$probabilities
  probs2 <- thresh_result2$probabilities
  true   <- thresh_result1$true
  binary <- as.numeric(true == thresh_result1$positive)
  
  # Compute NRI and IDI via nribin (single call)
  nri_obj <- tryCatch(
    nricens::nribin(event = binary, p.std = probs2, p.new = probs1,
                    cut = cutoffs, updown = "category"),
    error = function(e) NULL
  )
  
  # Extract estimates and standard errors for NRI+ and IDI
  if (!is.null(nri_obj)) {
    nri_plus <- tryCatch(nri_obj$nri["NRI+", "Estimate"],   error = function(e) NA_real_)
    nri_se   <- tryCatch(nri_obj$nri["NRI+", "Std.Error"],  error = function(e) NA_real_)
    idi_val  <- tryCatch(nri_obj$idi["IDI", "Estimate"],    error = function(e) NA_real_)
    idi_se   <- tryCatch(nri_obj$idi["IDI", "Std.Error"],   error = function(e) NA_real_)
  } else {
    nri_plus <- nri_se <- idi_val <- idi_se <- NA_real_
  }
  
  # Build data frame with percentages and 95% CIs
  df <- data.frame(
    Metric   = c("NRI+", "IDI"),
    Estimate = c(nri_plus * 100, idi_val * 100),
    Lower    = c((nri_plus - 1.96 * nri_se) * 100,
                 (idi_val  - 1.96 * idi_se) * 100),
    Upper    = c((nri_plus + 1.96 * nri_se) * 100,
                 (idi_val  + 1.96 * idi_se) * 100)
  )
  df <- df[!is.na(df$Estimate), ]
  
  # Create the plot
  if (nrow(df) == 0) {
    message("NRI/IDI could not be estimated.")
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "Not estimable", size = 6) +
      ggplot2::theme_void()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Estimate, y = Metric)) +
      # Vertical reference line at zero
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey70", linewidth = 0.8) +
      # Error bars (95% CI)
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = Lower, xmax = Upper),
        height = 0.15,
        linewidth = 1.0,
        color = "#2166ac"
      ) +
      # Point estimates
      ggplot2::geom_point(size = 3, color = "#2166ac", shape = 16) +
      # Labels and titles
      ggplot2::labs(
        title = paste("NRI/IDI:", label1, "vs", label2),
        subtitle = "Error bars represent 95% CI",
        x = "Value (%)",
        y = NULL
      ) +
      # Shaded region for clinical equivalence (+/-5%)
      ggplot2::annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf,
                        fill = "grey90", alpha = 0.3) +
      # Academic theme
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40", size = 10),
        axis.title.x  = ggplot2::element_text(face = "bold"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.background   = ggplot2::element_rect(fill = "white", color = NA),
        plot.background    = ggplot2::element_rect(fill = "white", color = NA)
      )
  }
  
  # Display the plot
  print(p)
  
  # Save to PDF if requested
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "nri_idi.pdf"),
                    plot = p, width = 6, height = 3, dpi = 300)
  }
  
  # Return components invisibly
  invisible(list(nri = nri_obj, idi = nri_obj$idi, plot = p))
}
#' Apply a Threshold to Predicted Probabilities
#'
#' Uses a specified threshold from \code{CalculateThresholds} (or a custom value)
#' to convert probabilities into binary predictions, then computes a
#' confusion matrix and common performance metrics.
#'
#' @param thresh_result Output of \code{CalculateThresholds}.
#' @param which_threshold Name of the threshold to use (e.g., "Youden", "PPV_Target").
#'   Ignored if \code{custom_threshold} is provided.
#' @param custom_threshold Numeric threshold value to override \code{which_threshold}.
#' @param newdata Optional data frame (if NULL, uses the true labels stored in thresh_result).
#' @param positive_class Positive class name (optional, auto-detected).
#' @return Invisibly returns a list with \code{predictions}, \code{conf_matrix}, \code{metrics}.
#' @export
ApplyThreshold <- function(thresh_result,
                           which_threshold = "Youden",
                           custom_threshold = NULL,
                           newdata = NULL,
                           positive_class = NULL) {
  prob <- thresh_result$probabilities
  true <- thresh_result$true
  positive <- if (!is.null(positive_class)) positive_class else thresh_result$positive
  negative <- setdiff(levels(true), positive)
  
  thr <- if (!is.null(custom_threshold)) custom_threshold else thresh_result$thresholds[which_threshold]
  
  pred_class <- factor(ifelse(prob > thr, positive, negative),
                       levels = c(negative, positive))
  # If newdata is provided, extract true labels from it (using same outcome column as model_obj)
  if (!is.null(newdata)) {
    # Assumes the group column is stored in thresh_result? Actually we need the outcome column.
    # Better: use thresh_result$true directly, but if newdata has different labels, we handle it.
    # For simplicity, we require that the outcome variable name be passed, but here we assume
    # the user wants to use the same true labels as the threshold calculation.
    # We'll just use thresh_result$true.
    warning("newdata argument not fully implemented; using original true labels.")
  }
  
  cm <- caret::confusionMatrix(pred_class, true, positive = positive)
  metrics <- data.frame(
    Threshold = thr,
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"],
    Accuracy = cm$overall["Accuracy"],
    F1 = cm$byClass["F1"]
  )
  cat(sprintf("--- Applied Threshold = %.4f ---\n", thr))
  print(cm$table)
  print(metrics)
  invisible(list(predictions = pred_class, conf_matrix = cm, metrics = metrics))
}

#' Compare Classification Performance at Specified Thresholds (Bar Plot)
#'
#' Applies user-chosen thresholds (or Youden) to two sets of predicted
#' probabilities, computes confusion matrices and common metrics, and displays
#' a side-by-side bar plot.
#'
#' @param thresh_result1 First threshold result (from \code{CalculateThresholds}).
#' @param thresh_result2 Second threshold result.
#' @param thr1 Threshold for model 1. If \code{NULL}, uses \code{"Youden"} from result1.
#' @param thr2 Threshold for model 2. If \code{NULL}, uses \code{"Youden"} from result2.
#' @param label1,label2 Model labels.
#' @param palette_name Wes Anderson palette for the two models (default \code{"Darjeeling1"}).
#' @param save_plot Save plot?
#' @param save_dir Output directory.
#' @return Invisible list with metrics and plots.
#' @export
CompareClassification <- function(thresh_result1,
                                  thresh_result2,
                                  thr1 = NULL,
                                  thr2 = NULL,
                                  label1 = "Model 1",
                                  label2 = "Model 2",
                                  palette_name = "Darjeeling1",
                                  save_plot = FALSE,
                                  save_dir = NULL) {
  .check_clinical_pkgs()
  # Resolve thresholds
  if (is.null(thr1)) thr1 <- thresh_result1$thresholds["Youden"]
  if (is.null(thr2)) thr2 <- thresh_result2$thresholds["Youden"]
  
  cat(sprintf("Threshold %s: %.4f\n", label1, thr1))
  cat(sprintf("Threshold %s: %.4f\n", label2, thr2))
  
  # Apply thresholds
  pred1 <- factor(ifelse(thresh_result1$probabilities > thr1, 
                         thresh_result1$positive, thresh_result1$negative),
                  levels = c(thresh_result1$negative, thresh_result1$positive))
  pred2 <- factor(ifelse(thresh_result2$probabilities > thr2,
                         thresh_result2$positive, thresh_result2$negative),
                  levels = c(thresh_result2$negative, thresh_result2$positive))
  true <- thresh_result1$true   # same labels for both
  
  cm1 <- caret::confusionMatrix(pred1, true, positive = thresh_result1$positive)
  cm2 <- caret::confusionMatrix(pred2, true, positive = thresh_result2$positive)
  
  # Extract metrics
  extract_metrics <- function(cm, model_name) {
    data.frame(Model = model_name,
               Sensitivity = cm$byClass["Sensitivity"],
               Specificity = cm$byClass["Specificity"],
               PPV = cm$byClass["Pos Pred Value"],
               NPV = cm$byClass["Neg Pred Value"],
               Accuracy = cm$overall["Accuracy"],
               F1 = cm$byClass["F1"],
               row.names = NULL)
  }
  metrics_df <- rbind(extract_metrics(cm1, label1),
                      extract_metrics(cm2, label2))
  
  # Long format for bar plot
  long_df <- tidyr::pivot_longer(metrics_df, -Model, names_to = "Metric", values_to = "Value")
  
  cols <- wesanderson::wes_palette(palette_name, 2, type = "discrete")
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = Metric, y = Value, fill = Model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", Value)),
      position = ggplot2::position_dodge(width = 0.8),
      vjust = -0.5, size = 3.5, colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::labs(title = "Classification Performance at Chosen Thresholds",
                  x = NULL, y = "Value") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "compare_classification.pdf"),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  
  invisible(list(metrics = metrics_df, cm1 = cm1, cm2 = cm2, plot = p))
}

#' Clinical Correlation Analysis (Standalone)
#'
#' Convenience wrapper around \code{\link{PlotClinicalCorrelation}}.
#' Automatically extracts clinical data from the model object if not provided,
#' and passes all additional arguments to \code{PlotClinicalCorrelation} to
#' compute and visualize associations between model predictions and clinical
#' variables.
#'
#' @param model_obj A \code{Train_Model} S4 object containing a trained model,
#'   and optionally clinical data attached via \code{\link{AttachClinicalData}}.
#' @param clinical_data Optional data frame with clinical variables. If
#'   \code{NULL} (default), the function attempts to use
#'   \code{model_obj@process.info$clinical_data}.
#' @param ... Additional arguments passed to \code{\link{PlotClinicalCorrelation}},
#'   such as:
#'   \describe{
#'     \item{newdata}{Data frame containing the feature data used for prediction.
#'       Required if not using the default from the model object.}
#'     \item{dataset_type}{One of \code{"testing"}, \code{"training"}, or
#'       \code{"external"}. Default is \code{"testing"}.}
#'     \item{model_name}{Character. Specific model name from
#'       \code{model_obj@train.models}, or \code{"ensemble"} to use an ensemble.
#'       Defaults to the best model.}
#'     \item{ordinal_vars}{Character vector of column names in \code{clinical_data}
#'       that should be treated as ordinal (included in Spearman correlation).}
#'     \item{palette_name}{Character. Wes Anderson palette name for boxplots.
#'       Default \code{"Royal1"}.}
#'     \item{save_plot}{Logical. Save plots as PDF? Default \code{FALSE}.}
#'     \item{save_dir}{Character. Directory to save plots if \code{save_plot = TRUE}.}
#'   }
#'
#' @return Invisibly returns a list with two components:
#'   \itemize{
#'     \item \code{spearman}: Spearman correlation matrix and p-values for
#'       continuous/ordinal variables (or \code{NULL} if none).
#'     \item \code{categorical}: Data frame of group-difference test results
#'       (Wilcoxon or Kruskal-Wallis) for nominal variables (or \code{NULL}).
#'   }
#'   The correlation heatmap and boxplots are drawn as side effects.
#'
#' @seealso \code{\link{PlotClinicalCorrelation}}, \code{\link{AttachClinicalData}}
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare binary classification model
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#'
#' # Attach clinical data (using other mtcars columns as "clinical")
#' clin <- data.frame(row.names = rownames(mtcars),
#'                    wt = mtcars$wt,
#'                    cyl = factor(mtcars$cyl))
#' model <- AttachClinicalData(model, clinical_data = clin)
#'
#' # Run correlation analysis
#' ClinicalCorrelation(model, dataset_type = "training", save_plot = FALSE)
#' }
ClinicalCorrelation <- function(model_obj, clinical_data = NULL, ...) {
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  PlotClinicalCorrelation(model_obj, clinical_data, ...)
}

#' Subgroup Analysis (Standalone)
#'
#' Convenience wrapper around \code{\link{PlotSubgroupForest}}.
#' Iterates over one or more categorical clinical variables and generates
#' a publication‑style forest plot showing subgroup‑specific AUCs with
#' 95% confidence intervals. Clinical data is automatically extracted
#' from the model object if not provided.
#'
#' @param model_obj A \code{Train_Model} S4 object containing a trained model,
#'   and optionally clinical data attached via \code{\link{AttachClinicalData}}.
#' @param clinical_data Optional data frame with clinical variables. If
#'   \code{NULL} (default), the function attempts to use
#'   \code{model_obj@process.info$clinical_data}.
#' @param subgroup_vars Character vector specifying the categorical column names
#'   in \code{clinical_data} to be used for subgroup stratification.
#'   Each unique level of each variable becomes a separate subgroup.
#'   This argument is \strong{required}.
#' @param ... Additional arguments passed to \code{\link{PlotSubgroupForest}},
#'   such as:
#'   \describe{
#'     \item{newdata}{Data frame used for model prediction. If \code{NULL},
#'       the default data from the model object is used.}
#'     \item{dataset_type}{One of \code{"testing"}, \code{"training"}, or
#'       \code{"external"}. Default is \code{"testing"}.}
#'     \item{model_name}{Character. Specific model or \code{"ensemble"}.
#'       Defaults to the best model.}
#'     \item{var_labels}{Named character vector to relabel section headers.}
#'     \item{level_order}{Named list specifying display order of factor levels.}
#'     \item{min_n}{Minimum subgroup size required for evaluation. Default \code{5}.}
#'     \item{xlim}{Numeric length-2 vector for x-axis limits.}
#'     \item{ref_line}{Reference vertical line (default \code{0.5} for AUC).}
#'     \item{save_plot}{Logical. Save plot as PDF? Default \code{FALSE}.}
#'     \item{save_dir}{Character. Directory to save the plot.}
#'   }
#'
#' @return Invisibly returns a list for the last generated forest plot,
#'   containing the \code{plot} (forestploter graphic), \code{data}
#'   (summary data frame), and \code{fits} (Cox regression fits, if applicable).
#'
#' @seealso \code{\link{PlotSubgroupForest}}, \code{\link{AttachClinicalData}}
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare model
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#'
#' # Attach clinical data (subgroup variables)
#' clin <- data.frame(row.names = rownames(mtcars),
#'                    cyl = factor(mtcars$cyl),
#'                    vs = factor(mtcars$vs))
#' model <- AttachClinicalData(model, clinical_data = clin)
#'
#' # Subgroup analysis by cyl and vs
#' ClinicalSubgroup(model, subgroup_vars = c("cyl", "vs"),
#'                  save_plot = FALSE)
#' }
ClinicalSubgroup <- function(model_obj, clinical_data = NULL,
                             subgroup_vars = NULL, ...) {
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  for (v in subgroup_vars)
    PlotSubgroupForest(model_obj, clinical_data, v, ...)
}

#' Confounder Adjustment Analysis (Standalone)
#'
#' This is a convenience wrapper around \code{\link{PlotConfounderForest}}.
#' It extracts clinical data from the model object if not provided, and passes
#' all additional arguments to \code{PlotConfounderForest} to fit a logistic
#' regression of the outcome on the model's predicted logit plus clinical
#' covariates, then plots the -log10(p-values) as a bar chart.
#'
#' @param model_obj A \code{Train_Model} S4 object containing a trained model,
#'   and optionally clinical data attached via \code{\link{AttachClinicalData}}.
#' @param clinical_data Optional data frame with clinical variables. If
#'   \code{NULL} (default), the function attempts to use
#'   \code{model_obj@process.info$clinical_data}.
#' @param ... Additional arguments passed to \code{\link{PlotConfounderForest}},
#'   such as:
#'   \describe{
#'     \item{dataset_type}{One of \code{"testing"}, \code{"training"}, or
#'       \code{"external"}. Determines which feature/clinical data are used.
#'       Default is \code{"testing"}.}
#'     \item{outcome_var}{Character. Name of the binary outcome column in the
#'       data. This argument is \strong{required}; must be supplied either here
#'       or via \code{...}.}
#'     \item{adjust_vars}{Character vector of clinical column names to adjust
#'       for. If \code{NULL} (default), all clinical columns except
#'       \code{outcome_var} are used.}
#'     \item{model_name}{Optional character. Specific model from
#'       \code{model_obj@train.models} to use for prediction, or
#'       \code{"ensemble"} to use an ensemble. Defaults to the best model.}
#'     \item{positive_class}{Character. The positive outcome level. Defaults
#'       to \code{"1"}.}
#'     \item{save_plot}{Logical. Whether to save the plot as PDF. Default
#'       \code{FALSE}.}
#'     \item{save_dir}{Character. Directory to save the plot if
#'       \code{save_plot = TRUE}.}
#'   }
#'
#' @return A \code{ggplot} object (invisibly printed). The plot shows
#'   -log10(p-values) for each covariate from the logistic regression model,
#'   with a dashed line at p = 0.05 for reference.
#'
#' @seealso \code{\link{PlotConfounderForest}}, \code{\link{AttachClinicalData}}
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#'
#' # Attach clinical data (here using other columns as "clinical")
#' clin <- data.frame(row.names = rownames(mtcars),
#'                    wt = mtcars$wt,
#'                    qsec = mtcars$qsec)
#' model <- AttachClinicalData(model, clinical_data = clin)
#'
#' # Run confounder analysis (outcome_var is required)
#' ClinicalConfounder(model, outcome_var = "am", save_plot = FALSE)
#' }
ClinicalConfounder <- function(model_obj, clinical_data = NULL, ...) {
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  PlotConfounderForest(model_obj, clinical_data, ...)
}

#' Threshold Analysis (Standalone) - Dual Mode
#'
#' Can be used in two ways:
#' 1. Provide \code{model_obj} (and optional \code{model_name}, \code{newdata}),
#'    thresholds are calculated internally and then visualized.
#' 2. Provide a pre-computed \code{thresh} object (from
#'    \code{CalculateThresholds} or \code{CalculateThresholdsFromProbs}),
#'    all plots are generated directly.  \code{model_obj} can be \code{NULL}.
#'
#' @param model_obj  A Train_Model or caret model.  If \code{thresh} is supplied,
#'   this is ignored.
#' @param thresh     Optional pre-computed threshold result.  When provided,
#'   \code{model_obj}, \code{newdata}, \code{model_name} and all \code{...}
#'   arguments are not used.
#' @param newdata    Data frame for prediction. Required when \code{model_obj}
#'   is provided and \code{thresh} is \code{NULL}. Defaults to the testing set
#'   from the \code{Train_Model} object if available.
#' @param model_name Character string specifying which model inside
#'   \code{model_obj} to use. Only applicable when \code{model_obj} is a
#'   \code{Train_Model} with multiple trained models. If \code{NULL}, the best
#'   model is used.
#' @param compare_model  Optional second \code{thresh_result} for ROC comparison and NRI.
#' @param compare_label  Label for the comparison model.
#' @param save_plot  Save plots?
#' @param save_dir   Output directory.
#' @param ...        Further arguments passed to \code{CalculateThresholds} when
#'   \code{thresh} is not supplied (e.g., \code{target_ppv}, \code{target_npv}).
#'
#' @return Invisible list with threshold results.
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # Example 1: Use a Train_Model object directly
#' result <- ClinicalThreshold(model_obj = model_obj,
#'                             newdata = test_data,
#'                             save_plot = FALSE)
#'
#' # Example 2: Use a pre-computed threshold object
#' thresh <- CalculateThresholds(model_obj, test_data)
#' result <- ClinicalThreshold(thresh = thresh)
#' }
ClinicalThreshold <- function(model_obj = NULL,
                              thresh = NULL,
                              newdata = NULL,
                              model_name = NULL,
                              compare_model = NULL,
                              compare_label = "Comparator",
                              save_plot = FALSE,
                              save_dir = NULL,
                              ...) {
  .check_clinical_pkgs()
  
  if (!is.null(thresh)) {
    if (!is.list(thresh) || is.null(thresh$thresholds))
      stop("'thresh' must be a list returned by CalculateThresholds or CalculateThresholdsFromProbs.")
    cat("Using pre-computed threshold object...\n")
  } else {
    if (is.null(model_obj))
      stop("Either 'model_obj' or 'thresh' must be provided.")
    if (is.null(newdata) && inherits(model_obj, "Train_Model"))
      newdata <- model_obj@filtered.set$testing
    thresh <- CalculateThresholds(model_obj, newdata, model_name, ...)
  }
  PlotThresholdAccuracy(thresh, save_plot = save_plot, save_dir = save_dir)
  PlotThresholdDensity(thresh, save_plot = save_plot, save_dir = save_dir)
  for (thr_name in names(thresh$thresholds)) {
    PlotThresholdWaterfall(thresh, thr_name, save_plot = save_plot, save_dir = save_dir)
    PlotThresholdConfusion(thresh, thr_name, save_plot = save_plot, save_dir = save_dir)
  }
  PlotThresholdROC(thresh, compare_model = compare_model,
                   compare_label = compare_label, save_plot = save_plot, save_dir = save_dir)
  if (!is.null(compare_model)) {
    CalculateNRI(thresh, compare_model,
                 label1 = "Model", label2 = compare_label,
                 save_plot = save_plot, save_dir = save_dir)
  }
  
  invisible(thresh)
}
#' Compare Two Models via Threshold ROC and NRI/IDI (Standalone)
#'
#' @param model_obj1 First model (Train_Model or caret train).
#' @param model_obj2 Second model (Train_Model or caret train).
#' @param newdata Data frame for prediction (default: testing set).
#' @param model_name1,model_name2 Optional model names within Train_Model objects.
#' @param label1,label2 Labels for the two models.
#' @param save_plot Logical. Save plots?
#' @param save_dir Output directory.
#' @param ... Additional arguments passed to `CalculateThresholds` for the first model
#'   (e.g., target_ppv, target_npv).
#' @return Invisible list with threshold results and comparison plots.
#' @export
CompareModelThresholds <- function(model_obj1,
                                   model_obj2,
                                   newdata = NULL,
                                   model_name1 = NULL,
                                   model_name2 = NULL,
                                   label1 = "Model 1",
                                   label2 = "Model 2",
                                   save_plot = FALSE,
                                   save_dir = NULL,
                                   ...) {
  .check_clinical_pkgs()
  thresh1 <- CalculateThresholds(model_obj1, newdata, model_name1, ...)
  thresh2 <- CalculateThresholds(model_obj2, newdata, model_name2)
  
  PlotThresholdROC(thresh1, compare_model = thresh2,
                   compare_label = label2,
                   save_plot = save_plot, save_dir = save_dir)
  
  CalculateNRI(thresh1, thresh2,
               label1 = label1, label2 = label2,
               save_plot = save_plot, save_dir = save_dir)
  
  invisible(list(thresh1 = thresh1, thresh2 = thresh2))
}

#' @keywords internal
.align_clinical_and_newdata <- function(model_obj, 
                                        clinical_data = NULL, 
                                        newdata       = NULL, 
                                        dataset_type  = c("testing", "training", "external")) {
  dataset_type <- match.arg(dataset_type)
  
  if (is.null(newdata) && inherits(model_obj, "Train_Model")) {
    if (dataset_type == "testing") {
      newdata <- model_obj@filtered.set$testing
      if (is.null(newdata)) newdata <- model_obj@split.data$testing
    } else if (dataset_type == "training") {
      newdata <- model_obj@filtered.set$training
      if (is.null(newdata)) newdata <- model_obj@split.data$training
    } else if (dataset_type == "external") {
      newdata <- model_obj@process.info$external_data
      if (is.null(newdata)) stop("No external_data found in model_obj@process.info$external_data.")
    }
  }
  if (is.null(newdata)) stop("newdata could not be determined.")
  
  if (is.null(clinical_data) && inherits(model_obj, "Train_Model")) {
    if (dataset_type == "external" && !is.null(model_obj@process.info$external_clinical)) {
      clinical_data <- model_obj@process.info$external_clinical
    } else {
      clinical_data <- model_obj@process.info$clinical_data
    }
  }
  if (is.null(clinical_data)) stop("clinical_data is NULL. Please attach clinical data first.")
  
  sample_ids <- rownames(newdata)
  if (is.null(sample_ids)) stop("newdata must have rownames matching clinical_data.")
  
  common_samples <- intersect(sample_ids, rownames(clinical_data))
  if (length(common_samples) == 0) {
    stop("No overlapping rownames between newdata and clinical_data.")
  }
  
  if (length(common_samples) < length(sample_ids)) {
    warning(sprintf("[%s] Matched %d out of %d samples in newdata with clinical_data.", 
                    dataset_type, length(common_samples), length(sample_ids)))
  }
  
  matched_newdata  <- newdata[common_samples, , drop = FALSE]
  matched_clinical <- clinical_data[common_samples, , drop = FALSE]
  
  gc <- model_obj@group_col
  if (!is.null(gc) && !gc %in% colnames(matched_newdata) && gc %in% colnames(matched_clinical)) {
    matched_newdata[[gc]] <- matched_clinical[[gc]]
  }
  
  return(list(
    newdata       = matched_newdata,
    clinical_data = matched_clinical,
    samples       = common_samples
  ))
}