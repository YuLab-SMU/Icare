# ==============================================================================
# prognosis_03_deploy.R
# PrognosiX Deployment Module
#
# Provides:
#   1. Stat_to_PrognosiX()      -- robust Stat -> PrognosiX conversion
#   2. run_prognosis_pipeline() -- end-to-end survival analysis pipeline
#   3. Prog_deploy_dispatcher() -- risk prediction for new samples
#   4. New_Prog_Manager()       -- lightweight deployment manager object
#   5. launch_prog_deploy_app() -- interactive Shiny prediction terminal
# ==============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b


# ==============================================================================
# 1.  Stat -> PrognosiX Conversion
# ==============================================================================

#' Convert a Stat Object to a PrognosiX Object
#'
#' This function extracts the survival time and status columns from a \code{Stat}
#' object and builds a \code{PrognosiX} object. It ensures that:
#' \itemize{
#'   \item The \code{time} and \code{status} columns are moved to \code{info.data}.
#'   \item \code{clean.data} contains only numeric feature columns.
#'   \item All other columns (clinical, IDs, etc.) are preserved in \code{info.data}.
#' }
#' Missing values are handled according to \code{na_action}.
#'
#' @param stat_obj A \code{Stat} S4 object.
#' @param time_col Character string naming the survival time column.
#' @param status_col Character string naming the event status column (1 = event, 0 = censored).
#' @param na_action Character; one of \code{"omit"} (remove rows with any NA in features),
#'   \code{"impute_median"} (replace NAs with median/mode), or \code{"allow"} (keep NAs).
#' @param min_events Integer; minimum number of events required; a warning is issued if fewer.
#' @param verbose Logical; print progress messages.
#' @return A \code{PrognosiX} S4 object.
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' stat <- CreateStatObject(raw.data = veteran, clean.data = veteran,
#'                          group_col = "status", na.action = "allow")
#' prog <- Stat_to_PrognosiX(stat, "time", "status", na_action = "omit",
#'                           min_events = 10, verbose = TRUE)
#' }
Stat_to_PrognosiX <- function(stat_obj,
                              time_col,
                              status_col,
                              na_action  = c("omit", "impute_median", "allow"),
                              min_events = 20,
                              verbose    = TRUE) {
  
  na_action <- match.arg(na_action)
  .log <- function(fmt, ...) if (verbose) message(sprintf(fmt, ...))
  
  # ---- 1. Validate and extract data ----
  if (!inherits(stat_obj, "Stat"))
    stop("[Stat_to_PrognosiX] 'stat_obj' must be a 'Stat' S4 object.")
  
  .log("[1/6] Extracting data from Stat object...")
  core_data <- if (nrow(stat_obj@clean.data) > 0) {
    stat_obj@clean.data
  } else if (nrow(stat_obj@raw.data) > 0) {
    warning("[Stat_to_PrognosiX] clean.data is empty; using raw.data.")
    stat_obj@raw.data
  } else {
    stop("[Stat_to_PrognosiX] Both clean.data and raw.data are empty.")
  }
  
  info_data <- stat_obj@info.data
  # Ensure row names match; if info.data is empty, create one
  if (nrow(info_data) == 0) {
    info_data <- data.frame(row.names = rownames(core_data))
  } else {
    common <- intersect(rownames(core_data), rownames(info_data))
    if (length(common) == 0) {
      warning("[Stat_to_PrognosiX] Row names do not match between clean.data and info.data. info.data ignored.")
      info_data <- data.frame(row.names = rownames(core_data))
    } else {
      core_data <- core_data[common, , drop = FALSE]
      info_data <- info_data[common, , drop = FALSE]
    }
  }
  
  # ---- 2. Ensure time/status columns exist ----
  .log("[2/6] Verifying time / status columns...")
  miss <- setdiff(c(time_col, status_col), colnames(core_data))
  if (length(miss) > 0) {
    # Try to find them in info.data
    if (all(miss %in% colnames(info_data))) {
      core_data[, miss] <- info_data[, miss]
      .log("   Copied missing columns from info.data: %s", paste(miss, collapse = ", "))
    } else {
      stop(sprintf("[Stat_to_PrognosiX] Column(s) not found: %s\n  Available in core: %s",
                   paste(miss, collapse = ", "),
                   paste(colnames(core_data), collapse = ", ")))
    }
  }
  
  # ---- 3. Coerce types and remove invalid rows ----
  .log("[3/6] Coercing column types...")
  core_data[[time_col]]   <- suppressWarnings(
    as.numeric(gsub("[^0-9.-]", "", as.character(core_data[[time_col]]))))
  core_data[[status_col]] <- suppressWarnings(
    as.numeric(as.character(core_data[[status_col]])))
  
  bad_t <- !is.finite(core_data[[time_col]]) | core_data[[time_col]] <= 0
  if (any(bad_t)) {
    warning(sprintf("[Stat_to_PrognosiX] Removing %d row(s) with time <= 0 or NA.", sum(bad_t)))
    core_data <- core_data[!bad_t, , drop = FALSE]
    info_data <- info_data[!bad_t, , drop = FALSE]
  }
  bad_s <- is.na(core_data[[status_col]]) | !(core_data[[status_col]] %in% c(0, 1))
  if (any(bad_s)) {
    warning(sprintf("[Stat_to_PrognosiX] Removing %d row(s) where status not in {0,1}.", sum(bad_s)))
    core_data <- core_data[!bad_s, , drop = FALSE]
    info_data <- info_data[!bad_s, , drop = FALSE]
  }
  
  # ---- 4. Move time/status to info.data and clean feature columns ----
  .log("[4/6] Separating features and metadata...")
  # Ensure time/status are in info.data
  info_data[[time_col]]   <- core_data[[time_col]]
  info_data[[status_col]] <- core_data[[status_col]]
  
  # All other columns: decide which are features (numeric) and which are metadata
  all_cols <- colnames(core_data)
  feature_candidates <- setdiff(all_cols, c(time_col, status_col))
  
  # Convert characters to factors (mlr3 requires factors, not characters)
  for (col in feature_candidates) {
    if (is.character(core_data[[col]])) {
      core_data[[col]] <- factor(core_data[[col]])
    }
  }
  
  # Separate numeric features from metadata
  numeric_features <- sapply(core_data[, feature_candidates, drop = FALSE], is.numeric)
  feature_cols <- feature_candidates[numeric_features]
  meta_cols    <- feature_candidates[!numeric_features]
  
  # Move non-numeric metadata to info.data (keep them for clinical analyses)
  if (length(meta_cols) > 0) {
    for (col in meta_cols) {
      info_data[[col]] <- core_data[[col]]
    }
    .log("   Moved %d non-numeric columns to info.data: %s",
         length(meta_cols), paste(meta_cols, collapse = ", "))
  }
  
  # Keep only numeric feature columns in clean.data
  clean_data <- core_data[, c(feature_cols), drop = FALSE]
  
  # ---- 5. Handle missing values ----
  .log("[5/6] Missing values (strategy: %s)...", na_action)
  if (na_action == "omit") {
    n_before  <- nrow(clean_data)
    keep_idx <- complete.cases(clean_data)
    clean_data <- clean_data[keep_idx, , drop = FALSE]
    info_data <- info_data[keep_idx, , drop = FALSE]
    removed   <- n_before - nrow(clean_data)
    if (removed > 0)
      warning(sprintf("[Stat_to_PrognosiX] Removed %d/%d row(s) with NA in features.",
                      removed, n_before))
  } else if (na_action == "impute_median") {
    for (col in colnames(clean_data)) {
      n_na <- sum(is.na(clean_data[[col]]))
      if (n_na > 0) {
        if (is.numeric(clean_data[[col]])) {
          fill <- median(clean_data[[col]], na.rm = TRUE)
        } else {
          fill <- names(sort(table(clean_data[[col]]), decreasing = TRUE))[1]
        }
        clean_data[[col]][is.na(clean_data[[col]])] <- fill
        .log("    [impute] %-14s  %d NA -> %s", col, n_na, fill)
      }
    }
  } else { 
    .log("   NA values are kept (na.action = 'allow').")
  }
  
  # ---- 6. Event count ----
  .log("[6/6] Checking event count...")
  n_ev <- sum(info_data[[status_col]] == 1)
  n_to <- nrow(info_data)
  .log("    N = %d | Events = %d | Censoring = %.1f%%",
       n_to, n_ev, (1 - n_ev / n_to) * 100)
  if (n_ev < min_events)
    warning(sprintf("[Stat_to_PrognosiX] Only %d events (< min_events=%d). Results may be unstable.",
                    n_ev, min_events))
  
  # ---- 7. Build PrognosiX object ----
  .log("[OK] Building PrognosiX object...")
  prog_obj <- CreatePrognosiXObject(
    clean.data = clean_data,
    info.data  = info_data,
    time_col   = time_col,
    status_col = status_col)
  
  .log("[OK] Done. Features: %d | Samples: %d", ncol(clean_data), nrow(clean_data))
  prog_obj
}


# ==============================================================================
# 2.  Complete Prognosis Analysis Pipeline
# ==============================================================================

#' Run the Complete Prognosis Analysis Pipeline
#'
#' Chains all prognosis analysis steps into a single call. By default, all
#' evaluation plots (KM, time-AUC, nomogram) are computed on the **training set**
#' (apparent performance). If you provide an external validation set via
#' \code{val_data}, the pipeline will also evaluate on that set and report
#' unbiased performance.
#'
#' @param object         \code{PrognosiX} or \code{Stat} S4 object.
#' @param time_col       Survival time column (required for \code{Stat} input).
#' @param status_col     Event status column (required for \code{Stat} input).
#' @param learner_ids    mlr3 learner IDs to benchmark.
#' @param p_threshold    Univariate Cox p-value cut-off.
#' @param tuning_budget  Number of hyperparameter evaluations per model.
#' @param cutoff_method  Risk stratification cut-point method.
#' @param time_points    Nomogram prediction horizons.
#' @param output_dir     Root folder for all outputs.
#' @param seed           Random seed.
#' @param run_shap       Run SHAP explanation (default \code{FALSE}).
#' @param run_nomogram   Draw nomogram (default \code{TRUE}).
#' @param val_data       Optional external validation data frame.
#' @param subgroup_vars  Optional subgroup variables for forest-plot analysis.
#'
#' @return Invisible named list with all results.
#' @export
run_prognosis_pipeline <- function(
    object,
    time_col      = NULL,
    status_col    = NULL,
    learner_ids   = c("surv.coxph", "surv.cv_glmnet", "surv.ranger"),
    p_threshold   = 0.05,
    tuning_budget = 30,
    cutoff_method = c("median", "tertile", "quartile", "p_optimize", "maxstat"),
    time_points   = c(1, 3, 5),
    output_dir    = NULL,
    seed          = 2025,
    run_shap      = FALSE,
    run_nomogram  = TRUE,
    val_data      = NULL,
    subgroup_vars = NULL) {
  .require_pkgs(c("mlr3", "mlr3proba", "mlr3tuning", "mlr3tuningspaces", "paradox"))
  
  cutoff_method <- match.arg(cutoff_method)
  set.seed(seed)
  
  if (is.null(output_dir)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path("./icare_output", "m4",
                            paste0("Prognosis_Pipeline_", ts))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Helper functions
  .sdir <- function(n, nm) {
    p <- file.path(output_dir, sprintf("Step%02d_%s", n, nm))
    dir.create(p, recursive = TRUE, showWarnings = FALSE); p
  }
  .sp <- function(plt, path, w = 8, h = 6)
    tryCatch(ggplot2::ggsave(path, plt, width = w, height = h),
             error = function(e) NULL)
  
  message("\n", strrep("=", 65))
  message(" PrognosiX Pipeline  |  seed=", seed,
          "  |  output: ", output_dir)
  message(strrep("=", 65))
  
  # Step 0: object conversion
  if (inherits(object, "Stat")) {
    if (is.null(time_col) || is.null(status_col))
      stop("[run_prognosis_pipeline] time_col and status_col required for Stat input.")
    message("\n[Step 0] Stat -> PrognosiX...")
    prog_obj <- Stat_to_PrognosiX(object, time_col, status_col, "omit")
  } else if (inherits(object, "PrognosiX")) {
    prog_obj <- object
  } else {
    stop("[run_prognosis_pipeline] object must be Stat or PrognosiX.")
  }
  saveRDS(prog_obj, file.path(output_dir, "prog_obj_initial.rds"))
  
  # Step 1: univariate Cox feature filtering
  d1 <- .sdir(1, "Feature_Selection")
  message("\n[Step 1] Univariate Cox feature filtering (p < ", p_threshold, ")...")
  filter_result <- surv_filter_features_clinical(prog_obj, p_threshold = p_threshold)
  write.csv(filter_result$table, file.path(d1, "Univariate_Cox_Results.csv"),
            row.names = FALSE)
  .sp(filter_result$plot, file.path(d1, "Univariate_Feature_Selection.pdf"))
  
  task_filtered <- filter_result$task
  selected_feats <- task_filtered$feature_names
  message(sprintf("  [OK] %d feature(s): %s",
                  length(selected_feats), paste(selected_feats, collapse = ", ")))
  
  if (length(selected_feats) == 0) {
    message("  [!] 0 features at p<", p_threshold, " -- relaxing to p<0.2 ...")
    filter_result <- surv_filter_features_clinical(prog_obj, p_threshold = 0.2)
    task_filtered <- filter_result$task
    selected_feats <- task_filtered$feature_names
    if (length(selected_feats) == 0)
      stop("Still 0 features at p<0.2. Check your data.")
  }
  prog_obj@univariate.analysis <- filter_result
  prog_obj@survival.var <- list(selected = selected_feats)
  
  # Step 2: algorithm benchmark
  d2 <- .sdir(2, "Algorithm_Benchmark")
  message("\n[Step 2] Algorithm mlr3::benchmark (", paste(learner_ids, collapse = ", "), ")...")
  lrn_list <- Filter(Negate(is.null), lapply(learner_ids, function(lid) {
    tryCatch(surv_get_learner(lid, task_filtered),
             error = function(e) { warning("Skipping ", lid, ": ", e$message); NULL })
  }))
  bmr_result <- surv_run_algorithm_benchmark(task_filtered, learners_list = lrn_list)
  perf_tab <- bmr_result$table
  best_id <- as.character(
    perf_tab[order(-perf_tab$surv.cindex), ][1, "learner_id"])
  
  write.csv(perf_tab, file.path(d2, "Benchmark_Leaderboard.csv"), row.names = FALSE)
  .sp(bmr_result$plot, file.path(d2, "Algorithm_Benchmark.pdf"), w = 9, h = 5)
  message(sprintf("  [OK] Winner: %s  (C-index = %.4f)",
                  best_id, perf_tab[perf_tab$learner_id == best_id, "surv.cindex"]))
  
  # Step 3: hyperparameter tuning
  d3 <- .sdir(3, "Final_Model")
  message("\n[Step 3] Tuning ", best_id, " (budget=", tuning_budget, " evals)...")
  tune_res <- surv_train_and_tune(task_filtered, best_id,
                                  tuning_budget = tuning_budget, seed = seed)
  best_learner <- tune_res$learner
  saveRDS(best_learner, file.path(d3, "best_learner.rds"))
  write.csv(
    data.frame(Algorithm = best_id,
               CV_CIndex = round(tune_res$cv_performance, 4),
               N_Features = length(selected_feats),
               Params = paste(names(tune_res$best_params),
                              unlist(tune_res$best_params),
                              sep = "=", collapse = "; ")),
    file.path(d3, "Best_Model_Summary.csv"), row.names = FALSE)
  message(sprintf("  [OK] CV C-index = %.4f", tune_res$cv_performance))
  
  if (is.null(val_data)) {
    message("\n", strrep("!", 65))
    message("!! WARNING: No validation data provided (val_data = NULL).")
    message("!! All subsequent evaluation plots (KM, time-AUC, nomogram)")
    message("!! are computed on the TRAINING set (apparent performance).")
    message("!! These estimates are optimistic and not suitable for")
    message("!! publication without independent validation.")
    message("!! If you have a validation set, pass it via val_data.")
    message(strrep("!", 65), "\n")
  }
  
  # Step 4: KM risk stratification
  d4 <- .sdir(4, "Risk_KM")
  message("\n[Step 4] KM risk stratification (", cutoff_method, ")...")
  km_plot <- tryCatch({
    p <- surv_plot_risk_km(best_learner, task_filtered,
                           cutoff_method = cutoff_method, risk_table = TRUE)
    cutoff_val <- attr(p, "cutoffs_used")
    if (!is.null(cutoff_val)) {
      tune_res$best_params$cutoff <- cutoff_val
      message(sprintf("  [OK] Cutoff extracted: %.4f", cutoff_val))
    }
    pdf(file.path(d4, paste0("KM_", cutoff_method, ".pdf")), width = 8, height = 7)
    print(p); dev.off(); p
  }, error = function(e) {
    if (grDevices::dev.cur() > 1) dev.off()
    message("  [!] KM skipped: ", e$message); NULL
  })
  
  # Step 5: time-dependent AUC
  d5 <- .sdir(5, "Time_AUC")
  message("\n[Step 5] Time-dependent AUC...")
  auc_data <- tryCatch({
    adf <- surv_plot_time_dependent_auc(best_learner, task_filtered)
    write.csv(adf, file.path(d5, "Time_AUC_Data.csv"), row.names = FALSE)
    .sp(last_plot(), file.path(d5, "Time_Dependent_AUC.pdf"))
    adf
  }, error = function(e) {
    message("  [!] AUC skipped (need risksetROC): ", e$message); NULL
  })
  
  # Step 6: nomogram (optional)
  if (run_nomogram && requireNamespace("rms", quietly = TRUE)) {
    d6 <- .sdir(6, "Nomogram")
    message("\n[Step 6] Nomogram...")
    tryCatch({
      pdf(file.path(d6, "Clinical_Nomogram.pdf"), width = 13, height = 7)
      surv_generate_nomogram(task_filtered,
                             selected_features = head(selected_feats, 6),
                             time_points = time_points)
      dev.off(); message("  [OK] Nomogram saved.")
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) dev.off()
      message("  [!] Nomogram failed: ", e$message)
    })
  }
  
  # Step 7: SHAP (optional)
  if (run_shap) {
    d7 <- .sdir(7, "SHAP")
    message("\n[Step 7] SHAP explanation...")
    tryCatch({
      sr <- surv_explain_shap(best_learner, task_filtered)
      .sp(sr$plot, file.path(d7, "SHAP_Importance.pdf"))
      message("  [OK] SHAP saved.")
    }, error = function(e)
      message("  [!] SHAP skipped (need survex): ", e$message))
  }
  
  # Step 7.5: Validation set evaluation (if provided)
  if (!is.null(val_data)) {
    message("\n[Step 7.5] Evaluating on validation set...")
    val_clean <- val_data[, c(selected_feats, time_col, status_col), drop = FALSE]
    val_prog <- CreatePrognosiXObject(
      clean.data = val_clean[, selected_feats, drop = FALSE],
      info.data = val_clean[, c(time_col, status_col)],
      time_col = time_col,
      status_col = status_col
    )
    val_task <- surv_extract_task(val_prog)
    
    val_cindex <- best_learner$predict(val_task)$score(mlr3::msr("surv.cindex"))
    message(sprintf("  Validation C-index = %.4f", val_cindex))
    
    med_time <- median(prog_obj@survival.data[[time_col]], na.rm = TRUE)
    cal_val <- tryCatch({
      surv_plot_calibration(best_learner, val_task, time_point = med_time, apparent = FALSE)
    }, error = function(e) NULL)
    
    km_val <- NULL
    if (!is.null(tune_res$best_params$cutoff)) {
      km_val <- tryCatch({
        surv_plot_risk_km(best_learner, val_task, cutoff_method = "custom",
                          custom_cutoffs = tune_res$best_params$cutoff)
      }, error = function(e) NULL)
    }
    
    prog_obj@subgroup.risk$validation <- list(
      cindex = val_cindex,
      calibration_plot = cal_val,
      km_plot = km_val,
      n = nrow(val_data)
    )
    
    if (!is.null(cal_val)) .sp(cal_val, file.path(d4, "Calibration_Validation.pdf"))
    if (!is.null(km_val)) .sp(km_val, file.path(d4, "KM_Validation.pdf"))
  } else {
    message("\n[Step 7.5] Skipped (no val_data provided).")
  }
  
  # Step 8: write back and save
  message("\n[Step 8] Saving results to PrognosiX object...")
  prog_obj@best.model <- list(
    learner_id = best_id,
    learner = best_learner,
    best_params = tune_res$best_params,
    cv_cindex = tune_res$cv_performance,
    features = selected_feats)
  prog_obj@subgroup.risk <- list(benchmark_table = perf_tab)
  saveRDS(prog_obj, file.path(output_dir, "prog_obj_final.rds"))
  
  message("\n", strrep("=", 65))
  message(" [OK] Pipeline complete -> ", output_dir)
  message(strrep("=", 65))
  
  invisible(list(
    prog_obj = prog_obj,
    best_learner = best_learner,
    filter_result = filter_result,
    benchmark_table = perf_tab,
    km_plot = km_plot,
    auc_data = auc_data,
    output_dir = output_dir))
}


# ==============================================================================
# PrognosiX Deployment Module
# ==============================================================================
# Functions for deploying trained PrognosiX models via Shiny app,
# with support for custom thresholds (binary, tertile, multi-group),
# full UI text and color theming, and a production-ready manager.
# ---- 1. Theme customization ----
#' Set custom theme for the Prognosis Terminal
#' @param primary_color Main accent color (borders, buttons, headers)
#' @param background_color Main background color
#' @param sidebar_color Sidebar and header background
#' @param box_background Box background color
#' @param text_color General text color
#' @param label_color Label text color
#' @param run_button_gradient_start Left side of button gradient
#' @param run_button_gradient_end Right side of button gradient
#' @param risk_high_color Text color for "High Risk"
#' @param risk_medium_color Text color for "Medium Risk"
#' @param risk_low_color Text color for "Low Risk"
#' @param table_header_color Table header text color
#' @param table_row_hover_color Hover background for table rows
#' @param font_family CSS font family (e.g., "Arial, sans-serif")
#' @param font_size_base Base font size in px
#' @return Invisibly stores the theme in options("prog_app_theme_css")
#' @export
#' 
#' @examples
#' \dontrun{
#' set_prog_app_theme(primary_color = "#2c7fb8", background_color = "#f8f9fa")
#' }
set_prog_app_theme <- function(
    primary_color = "#2c7fb8",
    background_color = "#f8f9fa",
    sidebar_color = "#e9ecef",
    box_background = "#ffffff",
    text_color = "#212529",
    label_color = "#2c7fb8",
    run_button_gradient_start = "#2c7fb8",
    run_button_gradient_end = "#1d4e6e",
    risk_high_color = "#d9534f",
    risk_medium_color = "#f0ad4e",
    risk_low_color = "#5cb85c",
    table_header_color = "#2c7fb8",
    table_row_hover_color = "#212529",
    font_family = NULL,
    font_size_base = 14
) {
  css <- sprintf("
    .content-wrapper, .right-side { background: %s !important; }
    .skin-black .main-header .logo, .skin-black .main-header .navbar,
    .skin-black .main-sidebar { background: %s !important; }
    .box { background: %s !important; color: %s !important;
           border-top: 3px solid %s !important; border-radius: 8px; }
    .box-header .box-title { color: %s !important; font-weight: bold; }
    label { color: %s !important; font-weight: bold; font-size: %dpx; }
    .form-control { background: %s !important; color: %s !important;
                    border: 1px solid #ced4da !important; }
    .btn-run { background: linear-gradient(90deg, %s, %s) !important;
               color: #fff !important; font-size: 18px !important; font-weight: bold !important;
               padding: 12px 40px !important; border-radius: 8px !important; border: none !important;
               margin: 16px 0; box-shadow: 0 2px 6px rgba(0,0,0,0.1); }
    .risk-high { font-size: 48px; font-weight: 900; color: %s; }
    .risk-med { font-size: 48px; font-weight: 900; color: %s; }
    .risk-low { font-size: 48px; font-weight: 900; color: %s; }
    .score-lbl { font-size: 14px; color: %s; font-family: monospace; }
    .badge-box { background: %s; border: 1px solid %s; border-radius: 6px;
                 padding: 8px 16px; margin: 4px; display: inline-block; color: %s; }
    table.dataTable { background: %s !important; color: %s !important; }
    table.dataTable thead th { border-bottom: 1px solid %s !important; color: %s !important; }
    table.dataTable tbody tr { background: %s !important; }
    table.dataTable tbody tr:hover { background: %s !important; }
    h4 { color: %s; font-weight: bold; border-left: 3px solid %s; padding-left: 10px; }
    p { color: %s; }
  ",
                 background_color, sidebar_color,
                 box_background, text_color, primary_color,
                 primary_color,
                 label_color, font_size_base,
                 box_background, text_color,
                 run_button_gradient_start, run_button_gradient_end,
                 risk_high_color, risk_medium_color, risk_low_color,
                 text_color, box_background, sidebar_color, primary_color,
                 box_background, text_color, sidebar_color, table_header_color,
                 box_background, table_row_hover_color,
                 primary_color, primary_color, text_color
  )
  if (!is.null(font_family)) {
    css <- paste0("body, .box, label, .form-control, .btn-run, .score-lbl, .badge-box, table, h4, p { font-family: ", font_family, " !important; }\n", css)
  }
  options("prog_app_theme_css" = css)
  invisible(css)
}

#' Predefined theme: advanced grey (light background, dark text)
#' @examples
#' \dontrun{
#' use_app_theme_grey()
#' }
#' @export
use_app_theme_grey <- function() {
  set_prog_app_theme()
}

#' Predefined theme: dark (original dark background)
#' @export
#' @examples
#' \dontrun{
#' use_app_theme_dark()
#' }
use_app_theme_dark <- function() {
  set_prog_app_theme(
    primary_color = "#58a6ff",
    background_color = "#0d1117",
    sidebar_color = "#161b22",
    box_background = "#1c2128",
    text_color = "#e6edf3",
    label_color = "#58a6ff",
    run_button_gradient_start = "#1f6feb",
    run_button_gradient_end = "#388bfd",
    risk_high_color = "#f85149",
    risk_medium_color = "#d29922",
    risk_low_color = "#3fb950",
    table_header_color = "#58a6ff",
    table_row_hover_color = "#21262d",
    font_family = NULL,
    font_size_base = 14
  )
}

#' Predefined theme: light (clean white, blue accents)
#' @examples
#' \dontrun{
#' use_app_theme_light()
#' }
#' @export
use_app_theme_light <- function() {
  set_prog_app_theme(
    primary_color = "#007bff",
    background_color = "#ffffff",
    sidebar_color = "#f8f9fa",
    box_background = "#ffffff",
    text_color = "#212529",
    label_color = "#007bff",
    run_button_gradient_start = "#007bff",
    run_button_gradient_end = "#0056b3",
    risk_high_color = "#dc3545",
    risk_medium_color = "#fd7e14",
    risk_low_color = "#28a745",
    table_header_color = "#007bff",
    table_row_hover_color = "#f1f3f5",
    font_family = NULL,
    font_size_base = 14
  )
}

# ---- 2. Text customization ----
#' Set custom text for the Prognosis Terminal
#' @param ... Named arguments for text keys (see default list)
#' @examples
#' \dontrun{
#' set_prog_app_text(title = "My Prognosis App", prediction_portal = "Risk Calculator")
#' }
#' @export
set_prog_app_text <- function(...) {
  default_text <- list(
    title = "Prognosis Terminal",
    prediction_portal = "Prediction Portal",
    model_info = "Model Info",
    documentation = "Documentation",
    overview_title = "Project Overview",
    abstract = "Abstract",
    reference = "Reference",
    abstract_text = "Prognostic risk stratification from a validated survival model.",
    citation_text = "Icare R package - PrognosiX framework",
    input_box_title = "1. Input Samples",
    risk_strat_label = "Risk Stratification",
    median_choice = "Median -- High/Low (2 groups)",
    tertile_choice = "Tertile -- Low/Med/High (3 groups)",
    custom_choice = "Custom thresholds (e.g., 30, 60)",
    custom_thresholds_label = "Thresholds (comma-separated)",
    show_scores_check = "Show raw risk scores",
    batch_help = "Batch CSV: first column = SampleID; remaining columns = features.",
    batch_tab = "Batch Upload (CSV)",
    single_tab = "Single Sample",
    sample_id_label = "Sample ID:",
    upload_button_label = "Upload .csv",
    download_template_label = "Download Template",
    calculate_button = "CALCULATE RISK",
    results_title = "2. Results",
    risk_group_heading = "Risk Group",
    sample_table_heading = "Sample Table",
    export_csv_button = "Export CSV",
    model_summary_title = "Model Summary",
    algorithm_label = "Algorithm",
    cv_cindex_label = "CV C-index",
    training_n_label = "Training N",
    events_label = "Events",
    selected_features_title = "Selected Features",
    variable_glossary_title = "Variable Glossary",
    feature_col = "Feature",
    description_col = "Description",
    units_col = "Units"
  )
  user_text <- list(...)
  final_text <- utils::modifyList(default_text, user_text)
  options("prog_app_text" = final_text)
  invisible(final_text)
}

#' Get Custom Application Text
#'
#' Retrieves a specific text string from the global application text options.
#' This is an internal helper used by the PrognosiX deployment Shiny app to
#' fetch user-customized labels and messages.
#'
#' @param key Character string specifying the text key to retrieve
#'   (e.g., `"title"`, `"calculate_button"`). If `NULL`, the entire
#'   text list is returned.
#'
#' @return If `key` is provided, returns the corresponding text string,
#'   or an empty string if the key does not exist. If `key` is `NULL`,
#'   returns the full named list of all application texts.
#'
#' @keywords internal
#' @noRd
get_prog_app_text <- function(key = NULL) {
  txt <- getOption("prog_app_text", list())
  if (is.null(key)) return(txt)
  return(txt[[key]] %||% "")
}

# ---- 3. Core prediction functions ----
#' Robust Prediction for PrognosiX Objects with Imputation
#'
#' Predict prognostic risk scores while automatically handling data type mismatches and missing values.
#'
#' @param prog_obj A \code{PrognosiX} object containing a trained survival model.
#' @param newdata A data frame (or data.table) with the same features as used for training.
#' @param impute Logical. Should missing values in \code{newdata} be imputed using the
#'   reference data (median for numeric, mode for factor/character)? Default is \code{TRUE}.
#'
#' @return A numeric vector of risk scores (crank) for each observation in \code{newdata}.
#'   Higher values indicate higher predicted risk.
#'
#' @details
#' The function first extracts the training task from the \code{PrognosiX} object to obtain
#' the exact feature types (integer, numeric, factor) and factor levels. It then coerces
#' the columns of \code{newdata} to these types. If \code{impute = TRUE}, any remaining
#' \code{NA} values are filled using the reference data stored in \code{prog_obj@clean.data}
#' (median for numeric, mode for categorical). Dummy \code{time} and \code{status} columns
#' are added to satisfy \code{TaskSurv} requirements, and the learner's \code{predict}
#' method is called to obtain risk scores.
#'
#' This function resolves the \code{Mlr3ErrorInput} issue that often occurs when
#' the original \code{predict_prognosix} is used with a learner like \code{surv.ranger}.
#'
#' @examples
#' \dontrun{
#' # Assume 'prog' is a trained PrognosiX object and 'new_data' contains the features
#' risk <- predict_prognosix_robust(prog, new_data, impute = TRUE)
#' head(risk)
#' }
#'
#' @importFrom data.table as.data.table
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' stat <- CreateStatObject(raw.data = veteran, clean.data = veteran,
#'                          group_col = "status", na.action = "allow")
#' prog <- Stat_to_PrognosiX(stat, "time", "status", na_action = "omit",
#'                           min_events = 10, verbose = FALSE)
#' task <- surv_extract_task(prog)
#' lrn <- surv_get_learner("surv.coxph", task)
#' lrn$train(task)
#'
#' prog@best.model <- list(
#'   learner_id = "surv.coxph", learner = lrn,
#'   features = task$feature_names, cutoff = 1.0,
#'   decision_type = "binary", train_cols = task$feature_names
#' )
#'
#' # Build new-patient data by reusing two rows from the training data itself,
#' # guaranteeing every required feature column is present.
#' new_patients <- as.data.frame(task$data())[1:2, task$feature_names]
#' rownames(new_patients) <- c("patient_1", "patient_2")
#'
#' risk_scores <- predict_prognosix(prog, new_patients, impute = TRUE)
#' risk_scores
#'
#' groups <- predict_risk_groups(prog, new_patients, cutoff_method = "median")
#' groups
#'
#' manager <- New_Prog_Manager(prog)
#' }
predict_prognosix <- function(prog_obj, newdata, impute = TRUE) {
  .require_pkgs(c("mlr3", "mlr3proba", "mlr3tuning", "mlr3tuningspaces", "paradox"))
  if (!inherits(prog_obj, "PrognosiX")) {
    stop("prog_obj must be a PrognosiX object.")
  }
  if (length(prog_obj@best.model) == 0) {
    stop("No best.model found in the PrognosiX object.")
  }
  
  # Extract training task and features
  train_task <- surv_extract_task(prog_obj)
  features <- prog_obj@best.model$features
  train_task$select(features)
  
  # Check that all required features are present in newdata
  newdata <- as.data.frame(newdata)
  missing_feats <- setdiff(features, colnames(newdata))
  if (length(missing_feats) > 0) {
    stop("Missing features: ", paste(missing_feats, collapse = ", "))
  }
  
  # Subset to training features
  new_clean <- newdata[, features, drop = FALSE]
  
  # Reference data for imputation (original cleaned data)
  ref <- prog_obj@clean.data[, features, drop = FALSE]
  
  # Coerce column types based on training task and optionally impute
  for (feat in features) {
    # Determine target type from training task
    target_type <- train_task$feature_types[train_task$feature_types$id == feat, "type"]
    
    # Convert to the exact storage type expected by the learner
    if (target_type == "integer") {
      new_clean[[feat]] <- as.integer(new_clean[[feat]])
    } else if (target_type == "numeric") {
      new_clean[[feat]] <- as.numeric(new_clean[[feat]])
    } else if (target_type == "factor") {
      # Obtain factor levels from the training task
      lev <- train_task$levels(feat)[[1]]
      new_clean[[feat]] <- factor(as.character(new_clean[[feat]]), levels = lev)
    } else {
      # fallback for other types (character, logical, etc.)
      new_clean[[feat]] <- new_clean[[feat]]
    }
    
    # Impute missing values if requested
    if (impute && anyNA(new_clean[[feat]])) {
      if (target_type %in% c("numeric", "integer")) {
        fill <- median(ref[[feat]], na.rm = TRUE)
        new_clean[[feat]][is.na(new_clean[[feat]])] <- fill
      } else if (target_type == "factor") {
        fill <- names(sort(table(ref[[feat]]), decreasing = TRUE))[1]
        new_clean[[feat]][is.na(new_clean[[feat]])] <- fill
        # Ensure factor levels remain correct after replacement
        new_clean[[feat]] <- factor(new_clean[[feat]], levels = lev)
      } else {
        # For character or other types, use mode
        fill <- names(sort(table(ref[[feat]]), decreasing = TRUE))[1]
        new_clean[[feat]][is.na(new_clean[[feat]])] <- fill
      }
    }
  }
  
  # Add dummy time and status columns (required by TaskSurv)
  new_clean$time <- 1L
  new_clean$status <- 0L
  
  # Create prediction task
  pred_task <- mlr3proba::TaskSurv$new(
    id = "pred_task",
    backend = data.table::as.data.table(new_clean),
    time = "time",
    event = "status"
  )
  
  # Predict risk scores
  learner <- prog_obj@best.model$learner
  if ("crank" %in% learner$predict_types) {
    learner$predict_type <- "crank"
  }
  risk_scores <- learner$predict(pred_task)$crank
  
  return(risk_scores)
}

#' Predict risk groups with arbitrary thresholds
#' @param prog_obj A `PrognosiX` object
#' @param newdata Data frame of new samples
#' @param cutoff_method `"median"`, `"tertile"`, or `"custom"`
#' @param custom_cutoffs Numeric vector (length >=1). For length >2, labels become "Low Risk", "Medium Risk 1", ..., "High Risk"
#' @param return_scores Logical, include raw risk scores
#' @return Data frame with SampleID, risk_group, and optionally risk_score
#' @examples
#' \dontrun{
#' # Requires a trained PrognosiX object
#' # risk_df <- predict_risk_groups(prog_obj, new_data, cutoff_method = "median")
#' }
#' @export
predict_risk_groups <- function(prog_obj, newdata, 
                                cutoff_method = c("median", "tertile", "custom"),
                                custom_cutoffs = NULL,
                                return_scores = TRUE) {
  cutoff_method <- match.arg(cutoff_method)
  risk_scores <- predict_prognosix(prog_obj, newdata, impute = TRUE)
  
  train_data <- prog_obj@survival.data
  train_feats <- prog_obj@best.model$features %||% colnames(prog_obj@clean.data)
  train_scores <- predict_prognosix(prog_obj, train_data[, train_feats, drop = FALSE], impute = FALSE)
  
  cuts <- switch(cutoff_method,
                 median = median(train_scores, na.rm = TRUE),
                 tertile = quantile(train_scores, probs = c(1/3, 2/3), na.rm = TRUE),
                 custom = {
                   if (is.null(custom_cutoffs))
                     stop("custom_cutoffs required")
                   sort(custom_cutoffs)
                 })
  n_groups <- length(cuts) + 1
  if (n_groups == 2) {
    risk_group <- ifelse(risk_scores > cuts, "High Risk", "Low Risk")
  } else {
    labels <- c("Low Risk", paste("Medium Risk", 1:(n_groups-2)), "High Risk")
    breaks <- c(-Inf, cuts, Inf)
    risk_group <- as.character(cut(risk_scores, breaks = breaks, labels = labels))
  }
  result <- data.frame(SampleID = rownames(newdata), risk_group = risk_group, stringsAsFactors = FALSE)
  if (return_scores) result$risk_score <- round(risk_scores, 6)
  return(result)
}

# ---- 4. Deployment dispatcher and manager ----
#' Deployment dispatcher (legacy compatibility)
#' @param prog_train_obj A `PrognosiX` object
#' @param newdata Data frame of new samples
#' @param cutoff_method One of `"median"`, `"tertile"`, `"custom"`
#' @param custom_cutoffs Numeric vector for custom thresholds
#' @param return_scores Logical
#' @examples
#' \dontrun{
#' # Legacy compatibility wrapper
#' # result <- Prog_deploy_dispatcher(prog_obj, new_data, cutoff_method = "tertile")
#' }
#' @export
Prog_deploy_dispatcher <- function(prog_train_obj,
                                   newdata,
                                   cutoff_method = c("median", "tertile", "custom"),
                                   custom_cutoffs = NULL,
                                   return_scores = TRUE) {
  cutoff_method <- match.arg(cutoff_method)
  predict_risk_groups(prog_train_obj, newdata, cutoff_method, custom_cutoffs, return_scores)
}

#' Create a deployment manager
#' @param prog_train_obj A `PrognosiX` object with a trained model
#' @return An S3 object of class `"Prog_Manager"`
#' @examples
#' \dontrun{
#' # Requires a trained PrognosiX object
#' # mgr <- New_Prog_Manager(prog_obj)
#' # result <- mgr$prog_predict(new_data)
#' }
#' @export
New_Prog_Manager <- function(prog_train_obj) {
  if (!inherits(prog_train_obj, "PrognosiX"))
    stop("[New_Prog_Manager] Input must be a PrognosiX object.")
  if (length(prog_train_obj@best.model) == 0)
    stop("[New_Prog_Manager] No best.model. Run a training pipeline first.")
  
  feats <- prog_train_obj@best.model$features %||% colnames(prog_train_obj@clean.data)
  info <- list(
    algorithm  = prog_train_obj@best.model$learner_id,
    cv_cindex  = prog_train_obj@best.model$cv_cindex,
    features   = feats,
    n_features = length(feats),
    n_train    = nrow(prog_train_obj@survival.data),
    n_events   = sum(prog_train_obj@survival.data[[prog_train_obj@status_col]] == 1)
  )
  mgr <- list(
    trained_obj = prog_train_obj,
    model_info  = info,
    prog_predict = function(newdata,
                            cutoff_method = c("median", "tertile", "custom"),
                            custom_cutoffs = NULL,
                            return_scores = TRUE) {
      cutoff_method <- match.arg(cutoff_method)
      predict_risk_groups(prog_train_obj, newdata, cutoff_method, custom_cutoffs, return_scores)
    }
  )
  class(mgr) <- "Prog_Manager"
  cat("\n-- PrognosiX Manager ---------------------------------\n")
  cat(sprintf("  Algorithm  : %s\n",   info$algorithm))
  cat(sprintf("  CV C-index : %.4f\n", info$cv_cindex %||% NA))
  cat(sprintf("  Features   : %d  (%s%s)\n", info$n_features,
              paste(head(info$features, 4), collapse = ", "),
              if (info$n_features > 4) ", ..." else ""))
  cat(sprintf("  Training N : %d (events: %d)\n", info$n_train, info$n_events))
  cat("------------------------------------------------------\n\n")
  return(mgr)
}

# ---- 5. Shiny application ----
#' Launch the interactive Prognosis Terminal
#' @param prog_manager A `Prog_Manager` object
#' @param title Browser title (overrides custom text)
#' @param var_dict Data frame with columns Feature, Description, Units
#' @param project_info List with elements `abstract` and `citation`
#' @return Runs the Shiny app
#' @examples
#' \dontrun{
#' # Requires a Prog_Manager object
#' # mgr <- New_Prog_Manager(prog_obj)
#' # launch_prog_deploy_app(mgr)
#' }
#' @export
launch_prog_deploy_app <- function(prog_manager,
                                   title = NULL,
                                   var_dict = NULL,
                                   project_info = NULL) {
  if (!inherits(prog_manager, "Prog_Manager"))
    stop("Input must be a Prog_Manager object.")
  for (pkg in c("shiny", "shinydashboard", "DT"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("'%s' required: install.packages('%s')", pkg, pkg))
  # Apply default theme if none set
  if (is.null(getOption("prog_app_theme_css"))) use_app_theme_grey()
  
  txt <- get_prog_app_text()
  if (is.null(title)) title <- txt$title %||% "Prognosis Terminal"
  .get_text <- function(key, def) { val <- txt[[key]]; if (is.null(val)) def else val }
  
  info <- prog_manager$model_info
  train_data <- prog_manager$trained_obj@clean.data
  req_vars <- info$features
  default_vals <- sapply(train_data[, req_vars, drop = FALSE],
                         function(x) round(median(as.numeric(x), na.rm = TRUE), 3))
  
  ui <- shinydashboard::dashboardPage(
    skin = "black",
    shinydashboard::dashboardHeader(title = shiny::span(shiny::icon("heartbeat"), " ", title)),
    shinydashboard::dashboardSidebar(shinydashboard::sidebarMenu(
      shinydashboard::menuItem(.get_text("prediction_portal", "Prediction Portal"), tabName = "portal", icon = shiny::icon("desktop")),
      shinydashboard::menuItem(.get_text("model_info", "Model Info"), tabName = "model", icon = shiny::icon("chart-bar")),
      shinydashboard::menuItem(.get_text("documentation", "Documentation"), tabName = "docs", icon = shiny::icon("info-circle"))
    )),
    shinydashboard::dashboardBody(
      shiny::tags$head(shiny::tags$style(shiny::HTML(getOption("prog_app_theme_css")))),
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "portal",
                shiny::fluidRow(shinydashboard::box(width = 12, title = .get_text("overview_title", "Project Overview"),
                             shiny::column(8, shiny::h4(.get_text("abstract", "Abstract")),
                                    shiny::p(project_info$abstract %||% .get_text("abstract_text", "Prognostic risk stratification from a validated survival model."))),
                             shiny::column(4, shiny::h4(.get_text("reference", "Reference")),
                                    shiny::p(project_info$citation %||% .get_text("citation_text", "Icare R package - PrognosiX framework")))
                )),
                shiny::fluidRow(shinydashboard::box(title = .get_text("input_box_title", "1. Input Samples"), width = 12,
                             shiny::column(3,
                                    shiny::selectInput("cutoff_m", .get_text("risk_strat_label", "Risk Stratification"),
                                                choices = {
                                                  lbl_med <- .get_text("median_choice", "Median -- High/Low (2 groups)")
                                                  lbl_ter <- .get_text("tertile_choice", "Tertile -- Low/Med/High (3 groups)")
                                                  lbl_cus <- .get_text("custom_choice", "Custom thresholds (e.g., 0.3, 0.6)")
                                                  setNames(c("median", "tertile", "custom"), c(lbl_med, lbl_ter, lbl_cus))
                                                }),
                                    shiny::conditionalPanel(condition = "input.cutoff_m == 'custom'",
                                                     shiny::textInput("custom_thresholds", .get_text("custom_thresholds_label", "Thresholds (comma-separated)"),
                                                               value = "0.3, 0.6", placeholder = "e.g., 0.2, 0.5, 0.8")),
                                    shiny::checkboxInput("show_sc", .get_text("show_scores_check", "Show raw risk scores"), TRUE),
                                    shiny::hr(),
                                    shiny::helpText(.get_text("batch_help", "Batch CSV: first column = SampleID; remaining columns = features."))
                             ),
                             shiny::column(9,
                                    shiny::tabsetPanel(id = "input_mode",
                                                shiny::tabPanel(.get_text("batch_tab", "Batch Upload (CSV)"), shiny::br(),
                                                         shiny::fileInput("up_file", .get_text("upload_button_label", "Upload .csv"), accept = ".csv"),
                                                         shiny::downloadButton("dl_tpl", .get_text("download_template_label", "Download Template"), class = "btn-xs btn-info")),
                                                shiny::tabPanel(.get_text("single_tab", "Single Sample"), shiny::br(),
                                                         shiny::fluidRow(shiny::column(4, shiny::textInput("sid", .get_text("sample_id_label", "Sample ID:"), "SAMPLE_001")),
                                                                  lapply(seq_along(req_vars), function(i)
                                                                    shiny::column(4, shiny::numericInput(paste0("f_", req_vars[i]), req_vars[i],
                                                                                           value = default_vals[[i]], step = 0.01))))
                                                )
                                    )
                             )
                )),
                shiny::fluidRow(shiny::column(12, align = "center",
                                shiny::actionButton("go", .get_text("calculate_button", "CALCULATE RISK"), icon = shiny::icon("play-circle"), class = "btn-run"))),
                shiny::fluidRow(shinydashboard::box(title = .get_text("results_title", "2. Results"), width = 12,
                             shiny::column(4, align = "center", shiny::h4(.get_text("risk_group_heading", "Risk Group")),
                                    shiny::uiOutput("risk_ui"), shiny::hr(), shiny::uiOutput("score_ui")),
                             shiny::column(8,
                                    shiny::div(style = "display:flex;justify-content:space-between;align-items:center",
                                        shiny::h4(.get_text("sample_table_heading", "Sample Table")),
                                        shiny::downloadButton("dl_res", .get_text("export_csv_button", "Export CSV"), class = "btn-success btn-xs")),
                                    DT::dataTableOutput("res_tbl"))
                ))
        ),
        shinydashboard::tabItem(tabName = "model",
                shiny::fluidRow(shinydashboard::box(title = .get_text("model_summary_title", "Model Summary"), width = 12,
                             shiny::fluidRow(
                               shiny::column(3, shiny::div(class="badge-box", shiny::icon("brain"), " ", .get_text("algorithm_label", "Algorithm"), shiny::br(), shiny::tags$b(info$algorithm))),
                               shiny::column(3, shiny::div(class="badge-box", shiny::icon("chart-line"), " ", .get_text("cv_cindex_label", "CV C-index"), shiny::br(), shiny::tags$b(round(info$cv_cindex %||% NA, 4)))),
                               shiny::column(3, shiny::div(class="badge-box", shiny::icon("users"), " ", .get_text("training_n_label", "Training N"), shiny::br(), shiny::tags$b(info$n_train))),
                               shiny::column(3, shiny::div(class="badge-box", shiny::icon("flag"), " ", .get_text("events_label", "Events"), shiny::br(), shiny::tags$b(info$n_events)))
                             ),
                             shiny::br(), shiny::h4(.get_text("selected_features_title", "Selected Features")),
                             DT::dataTableOutput("feat_tbl")
                ))
        ),
        shinydashboard::tabItem(tabName = "docs",
                shiny::fluidRow(shinydashboard::box(title = .get_text("variable_glossary_title", "Variable Glossary"), width = 12,
                             DT::dataTableOutput("doc_tbl")))
        )
      )
    )
  )
  
  server <- function(input, output, session) {
    parsed_data <- shiny::eventReactive(input$go, {
      if (input$input_mode == "Batch Upload (CSV)") {
        shiny::req(input$up_file)
        read.csv(input$up_file$datapath, row.names = 1, check.names = FALSE)
      } else {
        vals <- setNames(sapply(req_vars, function(v) input[[paste0("f_", v)]]), req_vars)
        df <- as.data.frame(t(vals)); rownames(df) <- input$sid; df
      }
    })
    preds <- shiny::reactive({
      shiny::req(parsed_data())
      custom_cutoffs <- NULL
      if (input$cutoff_m == "custom") {
        thr_str <- gsub(" ", "", input$custom_thresholds)
        custom_cutoffs <- as.numeric(strsplit(thr_str, ",")[[1]])
        if (any(is.na(custom_cutoffs))) shiny::showNotification("Invalid custom thresholds", type = "error")
      }
      prog_manager$prog_predict(parsed_data(),
                                cutoff_method = input$cutoff_m,
                                custom_cutoffs = custom_cutoffs,
                                return_scores = input$show_sc)
    })
    output$risk_ui <- shiny::renderUI({
      shiny::validate(shiny::need(preds(), "Awaiting input..."))
      df <- preds()
      if (nrow(df) == 1) {
        g <- df$risk_group[1]
        cl <- if (grepl("High", g)) "risk-high" else if (grepl("Med", g)) "risk-med" else "risk-low"
        shiny::tagList(shiny::div(class = "score-lbl", "Sample: ", df$SampleID[1]), shiny::div(class = cl, g))
      } else {
        tbl <- table(df$risk_group)
        shiny::tagList(shiny::div(class = "score-lbl", sprintf("Batch: %d samples processed", nrow(df))), shiny::br(),
                lapply(names(tbl), function(g) {
                  cl <- if (grepl("High", g)) "risk-high" else if (grepl("Med", g)) "risk-med" else "risk-low"
                  shiny::div(shiny::span(class = cl, style = "font-size:28px", g),
                      shiny::span(style = "color:#8b949e;font-size:13px;margin-left:8px", paste0("n=", tbl[[g]])))
                }))
      }
    })
    output$score_ui <- shiny::renderUI({
      shiny::validate(shiny::need(preds(), ""))
      df <- preds()
      if (nrow(df) == 1 && "risk_score" %in% names(df))
        shiny::div(class = "score-lbl", sprintf("Risk Score: %.4f", df$risk_score[1]))
    })
    output$res_tbl <- DT::renderDataTable({
      shiny::req(preds())
      DT::datatable(preds(), rownames = FALSE, options = list(dom = "tp", pageLength = 10)) %>%
        DT::formatStyle("risk_group",
                        color = DT::styleEqual(c("High Risk","Medium Risk","Low Risk"),
                                               c("#f85149","#d29922","#3fb950")), fontWeight = "bold")
    })
    output$dl_res <- shiny::downloadHandler(filename = function() paste0("Prognosis_Results_", Sys.Date(), ".csv"),
                                     content = function(f) write.csv(preds(), f, row.names = FALSE))
    output$dl_tpl <- shiny::downloadHandler(filename = "Input_Template.csv",
                                     content = function(f) write.csv(head(train_data[, req_vars, drop = FALSE], 5), f))
    output$feat_tbl <- DT::renderDataTable(
      DT::datatable(data.frame(`#` = seq_along(req_vars), Feature = req_vars,
                               Median_Training = round(default_vals, 3), check.names = FALSE),
                    rownames = FALSE, options = list(dom = "t", pageLength = 30)))
    output$doc_tbl <- DT::renderDataTable({
      df <- var_dict %||% data.frame(Feature = req_vars,
                                     Description = paste("Predictor:", req_vars), Units = "--", stringsAsFactors = FALSE)
      colnames(df) <- c(.get_text("feature_col", "Feature"),
                        .get_text("description_col", "Description"),
                        .get_text("units_col", "Units"))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 30))
    })
  }
  shiny::shinyApp(ui, server)
}
