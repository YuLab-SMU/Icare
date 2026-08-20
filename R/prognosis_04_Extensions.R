# =============================================================================
# icare Package — Module 4: Advanced Extensions for PrognosiX (Full Suite)
# =============================================================================
# This module provides 10 publication‑grade, robust extensions that integrate
# seamlessly with your main prognosis pipeline. Each function performs its own
# package checks using requireNamespace, avoiding helper function conflicts.
#
# Functions:
#   1. run_competing_risks()          – Fine‑Gray or cause‑specific Cox
#   2. plot_tdROC_comparison()        – Time‑dependent ROC with CIs and comparison
#   3. plot_multitime_calibration()   – Calibration at multiple time points
#   4. repeated_cv_evaluation()       – Repeated CV with parallel support
#   5. compare_models_delong()        – DeLong test for C‑index difference
#   6. tune_with_bayes()              – Bayesian hyperparameter optimisation
#   7. stepwise_variable_selection()  – AIC stepwise Cox variable selection
#   8. calibration_in_the_large()     – Overall calibration (mean predicted vs observed)
#   9. external_validation_test()     – Survival Hosmer‑Lemeshow test
#  10. plot_clinical_impact()         – Clinical impact curve with uncertainty
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Competing Risks: Fine‑Gray or Cause‑Specific Cox (Robust)
# ---------------------------------------------------------------------------
#' Run a competing risks survival model with fallback and checks
#'
#' @param object A `PrognosiX` or `TaskSurv` object.
#' @param model_type Either `"finegray"` (requires mlr3extralearners) or `"cox"`.
#' @param cause Event code for the event of interest (default 1).
#' @param learner_id Override default learner ID.
#' @param min_events Minimum events required (default 10).
#' @param ... Additional arguments passed to `surv_train_and_tune`.
#' @return A list with `learner`, `cv_cindex`, and `model_type`.
#' @export
run_competing_risks <- function(object,
                                model_type = c("finegray", "cox"),
                                cause = 1,
                                learner_id = NULL,
                                min_events = 10,
                                ...) {
  model_type <- match.arg(model_type)
  
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    stop("Package 'mlr3proba' is required. Install with: install.packages('mlr3proba')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  time_col <- task$target_names[1]
  status_col <- task$target_names[2]
  
  n_events <- sum(data[[status_col]] == cause, na.rm = TRUE)
  if (n_events < min_events) {
    stop(sprintf("Only %d events for cause %d (need >= %d).", n_events, cause, min_events))
  }
  
  if (model_type == "finegray") {
    if (requireNamespace("mlr3extralearners", quietly = TRUE)) {
      if (is.null(learner_id)) learner_id <- "surv.fgr"
    } else {
      warning("mlr3extralearners not installed. Falling back to cause‑specific Cox.")
      model_type <- "cox"
      learner_id <- "surv.coxph"
    }
  } else {
    data[[status_col]] <- ifelse(data[[status_col]] == cause, 1, 0)
    task <- surv_create_surv_task(data, time_col, status_col)
    if (is.null(learner_id)) learner_id <- "surv.coxph"
  }
  
  # Remove constant features
  feature_cols <- task$feature_names
  for (f in feature_cols) {
    if (length(unique(data[[f]])) < 2) {
      warning("Removing constant feature: ", f)
      task$select(setdiff(task$feature_names, f))
    }
  }
  
  tune_res <- tryCatch({
    surv_train_and_tune(object = task, learner_id = learner_id, ...)
  }, error = function(e) {
    warning("Fine‑Gray model failed: ", e$message, ". Falling back to cause‑specific Cox.")
    task <- surv_create_surv_task(data, time_col, status_col)
    surv_train_and_tune(object = task, learner_id = "surv.coxph", ...)
  })
  
  list(learner = tune_res$learner,
       cv_cindex = tune_res$cv_performance,
       model_type = model_type,
       n_events = n_events)
}

# ---------------------------------------------------------------------------
# 2. Time‑dependent ROC with Confidence Intervals and Model Comparison
# ---------------------------------------------------------------------------
#' Plot time‑dependent ROC curves with bootstrap CIs and compare two models
#'
#' @param learner1,learner2 Two trained learners (if comparing); learner2 can be NULL.
#' @param object A `TaskSurv` or `PrognosiX` object.
#' @param time_points Vector of time points; if NULL, uses 25th, 50th, 75th percentiles.
#' @param n_boot Bootstrap iterations for CIs (default 100).
#' @param seed Random seed.
#' @param save_plot Logical; save to `sub_dir("03_evaluation")`.
#' @return A list with AUC values, CIs, and a ggplot.
#' @export
plot_tdROC_comparison <- function(learner1, learner2 = NULL,
                                  object,
                                  time_points = NULL,
                                  n_boot = 100,
                                  seed = 123,
                                  save_plot = TRUE) {
  if (!requireNamespace("timeROC", quietly = TRUE)) {
    stop("Package 'timeROC' is required. Install with: install.packages('timeROC')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }
  
  set.seed(seed)
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  time <- data[[task$target_names[1]]]
  status <- data[[task$target_names[2]]]
  
  if (is.null(time_points)) {
    event_times <- time[status == 1]
    if (length(event_times) < 3) stop("Fewer than 3 event times.")
    time_points <- round(quantile(event_times, probs = c(0.25, 0.5, 0.75)), 0)
  }
  
  get_crank <- function(learner) {
    if (!("crank" %in% learner$predict_types)) {
      learner$predict_type <- "crank"
    }
    learner$predict(task)$crank
  }
  
  marker1 <- get_crank(learner1)
  if (!is.null(learner2)) marker2 <- get_crank(learner2)
  
  .roc_at_times <- function(marker, times) {
    roc <- timeROC::timeROC(T = time, delta = status,
                            marker = marker, cause = 1,
                            weighting = "marginal", times = times)
    auc <- roc$AUC
    auc_boot <- matrix(NA, nrow = n_boot, ncol = length(times))
    for (b in seq_len(n_boot)) {
      idx <- sample(length(marker), replace = TRUE)
      tryCatch({
        roc_b <- timeROC::timeROC(T = time[idx], delta = status[idx],
                                  marker = marker[idx], cause = 1,
                                  weighting = "marginal", times = times)
        auc_boot[b, ] <- roc_b$AUC
      }, error = function(e) NA)
    }
    ci_lower <- apply(auc_boot, 2, quantile, 0.025, na.rm = TRUE)
    ci_upper <- apply(auc_boot, 2, quantile, 0.975, na.rm = TRUE)
    data.frame(time = times, AUC = auc, lower = ci_lower, upper = ci_upper)
  }
  
  df1 <- .roc_at_times(marker1, time_points)
  df1$Model <- "Model 1"
  if (!is.null(learner2)) {
    df2 <- .roc_at_times(marker2, time_points)
    df2$Model <- "Model 2"
    roc_df <- rbind(df1, df2)
  } else {
    roc_df <- df1
  }
  
  p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = time, y = AUC, colour = Model)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = Model),
                         alpha = 0.2, colour = NA) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(x = "Time", y = "AUC", title = "Time‑dependent ROC with 95% CI") +
    ggprism::theme_prism()
  
  print(p)
  
  if (save_plot) {
    if (exists("sub_dir", mode = "function")) {
      save_path <- sub_dir("03_evaluation")
    } else {
      save_path <- "./03_evaluation"
      dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    }
    ggplot2::ggsave(file.path(save_path, "tdROC_comparison.pdf"),
                    p, width = 8, height = 5)
  }
  
  invisible(list(plot = p, data = roc_df))
}

# ---------------------------------------------------------------------------
# 3. Multi‑time Calibration Curves (with automatic time selection)
# ---------------------------------------------------------------------------
#' Plot calibration curves at multiple time points (auto‑selected if not provided)
#'
#' @param learner A trained mlr3 learner (supports "distr").
#' @param object A `TaskSurv` or `PrognosiX` object.
#' @param time_points Numeric vector; if NULL, uses 25th, 50th, 75th percentiles of event times.
#' @param n_bins Number of bins per curve (adaptive).
#' @param save_plot Logical; save to `sub_dir("03_evaluation")`.
#' @param ... Additional arguments passed to `surv_plot_calibration`.
#' @return A `patchwork` composite plot.
#' @export
plot_multitime_calibration <- function(learner, object,
                                       time_points = NULL,
                                       n_bins = 10,
                                       save_plot = TRUE, ...) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  time <- data[[task$target_names[1]]]
  status <- data[[task$target_names[2]]]
  
  if (is.null(time_points)) {
    event_times <- time[status == 1]
    if (length(event_times) < 3) {
      stop("Fewer than 3 event times; cannot select quantiles.")
    }
    time_points <- round(quantile(event_times, probs = c(0.25, 0.5, 0.75)), 0)
    message("Using automatically selected time points: ", paste(time_points, collapse = ", "))
  }
  
  valid_times <- time_points[sapply(time_points, function(tp) {
    any(time[status == 1] >= tp)
  })]
  if (length(valid_times) == 0) {
    stop("None of the specified time points have any events.")
  }
  time_points <- valid_times
  
  cal_plots <- lapply(time_points, function(tp) {
    p <- surv_plot_calibration(learner, object,
                               time_point = tp,
                               n_bins = n_bins,
                               ...)
    p + ggplot2::labs(subtitle = paste("t =", tp))
  })
  
  combined <- patchwork::wrap_plots(cal_plots, ncol = min(3, length(cal_plots)))
  print(combined)
  
  if (save_plot) {
    if (exists("sub_dir", mode = "function")) {
      save_path <- sub_dir("03_evaluation")
    } else {
      save_path <- "./03_evaluation"
      dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    }
    ggplot2::ggsave(file.path(save_path, "calibration_multitime.pdf"),
                    combined, width = 10, height = 6)
    message("Multi‑time calibration plot saved to: ", save_path)
  }
  
  invisible(combined)
}

# ---------------------------------------------------------------------------
# 4. Repeated Cross‑Validation (with parallel support)
# ---------------------------------------------------------------------------
#' Perform repeated k‑fold CV with optional parallelisation
#'
#' @param object A `TaskSurv` or `PrognosiX` object.
#' @param learner_id Character; learner to evaluate.
#' @param n_repeats Integer; number of CV repetitions (default 5).
#' @param folds Integer; number of folds per repeat (default 5).
#' @param tune Logical; whether to tune within each CV loop (time‑consuming).
#' @param parallel Logical; use parallel processing (requires `foreach`).
#' @param n_cores Number of cores for parallel (default detectCores() - 1).
#' @param seed Random seed for reproducibility.
#' @return A list with results and summary (including bootstrap CI).
#' @export
repeated_cv_evaluation <- function(object, learner_id,
                                   n_repeats = 5, folds = 5,
                                   tune = FALSE,
                                   parallel = FALSE,
                                   n_cores = NULL,
                                   seed = 123) {
  if (!requireNamespace("mlr3", quietly = TRUE)) {
    stop("Package 'mlr3' is required. Install with: install.packages('mlr3')")
  }
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    stop("Package 'mlr3proba' is required. Install with: install.packages('mlr3proba')")
  }
  if (parallel && !requireNamespace("foreach", quietly = TRUE)) {
    warning("Package 'foreach' not installed. Running sequentially.")
    parallel <- FALSE
  }
  
  set.seed(seed)
  task <- surv_extract_task(object)
  measure <- mlr3::msr("surv.cindex")
  
  if (parallel && requireNamespace("doParallel", quietly = TRUE)) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else if (parallel) {
    warning("doParallel not installed. Running sequentially.")
    parallel <- FALSE
  }
  
  .run_cv <- function(rep_idx) {
    set.seed(seed + rep_idx)
    resampling <- mlr3::rsmp("cv", folds = folds)$instantiate(task)
    fold_cindex <- numeric(folds)
    for (f in seq_len(folds)) {
      train_set <- resampling$train_set(f)
      test_set <- resampling$test_set(f)
      train_task <- task$clone()$filter(train_set)
      test_task <- task$clone()$filter(test_set)
      
      if (tune) {
        learner <- surv_train_and_tune(train_task, learner_id,
                                       tuning_budget = 10)$learner
      } else {
        learner <- surv_get_learner(learner_id, train_task)$train(train_task)
      }
      pred <- learner$predict(test_task)
      fold_cindex[f] <- pred$score(measure)
    }
    data.frame(Repeat = rep_idx, Fold = seq_len(folds), C_index = fold_cindex)
  }
  
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    results <- foreach::foreach(r = seq_len(n_repeats), .combine = rbind) %dopar% {
      .run_cv(r)
    }
  } else {
    results <- do.call(rbind, lapply(seq_len(n_repeats), .run_cv))
  }
  
  all_cindex <- results$C_index
  summary <- data.frame(
    Mean = mean(all_cindex),
    SD = sd(all_cindex),
    CI_lower = quantile(all_cindex, 0.025),
    CI_upper = quantile(all_cindex, 0.975),
    N = length(all_cindex)
  )
  
  list(results = results, summary = summary)
}

# ---------------------------------------------------------------------------
# 5. DeLong Test for Comparing Two Models (with bootstrap CIs)
# ---------------------------------------------------------------------------
#' Compare two survival models using DeLong's test with confidence intervals
#'
#' @param model1,model2 Trained mlr3 learners (must support "crank").
#' @param task A `TaskSurv` object (validation set recommended).
#' @param n_boot Bootstrap iterations (default 1000).
#' @param seed Random seed for reproducibility.
#' @return A list with p‑value, C‑indices, difference, and bootstrap CI.
#' @export
compare_models_delong <- function(model1, model2, task,
                                  n_boot = 1000, seed = 123) {
  if (!requireNamespace("survcomp", quietly = TRUE)) {
    stop("Package 'survcomp' is required. Install with: install.packages('survcomp')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  set.seed(seed)
  
  pred1 <- model1$predict(task)$crank
  pred2 <- model2$predict(task)$crank
  time <- task$data()[[task$target_names[1]]]
  status <- task$data()[[task$target_names[2]]]
  
  valid <- complete.cases(pred1, pred2, time, status)
  if (sum(valid) < 10) stop("Too few valid observations for comparison.")
  pred1 <- pred1[valid]; pred2 <- pred2[valid]
  time <- time[valid]; status <- status[valid]
  
  ci1 <- survcomp::concordance.index(pred1, time, status, method = "noether")
  ci2 <- survcomp::concordance.index(pred2, time, status, method = "noether")
  
  test <- tryCatch({
    survcomp::cindex.comp(pred1, pred2, time, status,
                          conf.int = 0.95, nboot = n_boot)
  }, error = function(e) {
    warning("cindex.comp failed: ", e$message, ". Using simpler test.")
    diff_boot <- replicate(n_boot, {
      idx <- sample(length(pred1), replace = TRUE)
      ci1b <- survcomp::concordance.index(pred1[idx], time[idx], status[idx],
                                          method = "noether")$c.index
      ci2b <- survcomp::concordance.index(pred2[idx], time[idx], status[idx],
                                          method = "noether")$c.index
      ci1b - ci2b
    })
    p_val <- 2 * min(mean(diff_boot < 0), mean(diff_boot > 0))
    list(p.value = p_val, diff_ci = quantile(diff_boot, c(0.025, 0.975)))
  })
  
  list(p.value = test$p.value,
       C_index_model1 = ci1$c.index,
       C_index_model2 = ci2$c.index,
       difference = ci1$c.index - ci2$c.index,
       ci_lower = if (!is.null(test$diff_ci)) test$diff_ci[1] else NA,
       ci_upper = if (!is.null(test$diff_ci)) test$diff_ci[2] else NA,
       test_result = test)
}

# ---------------------------------------------------------------------------
# 6. Bayesian Optimisation (with fallback to random search)
# ---------------------------------------------------------------------------
#' Tune a survival learner using Bayesian optimisation (robust fallback)
#'
#' @param object A `TaskSurv` or `PrognosiX` object.
#' @param learner_id Character; learner to tune.
#' @param tuning_budget Number of evaluations.
#' @param ... Additional arguments passed to `mlr3tuning::tnr("mbo")`.
#' @return Same structure as `surv_train_and_tune`.
#' @export
tune_with_bayes <- function(object, learner_id, tuning_budget = 30, ...) {
  if (!requireNamespace("mlr3tuning", quietly = TRUE)) {
    stop("Package 'mlr3tuning' is required. Install with: install.packages('mlr3tuning')")
  }
  if (!requireNamespace("mlr3mbo", quietly = TRUE)) {
    stop("Package 'mlr3mbo' is required. Install with: install.packages('mlr3mbo')")
  }
  if (!requireNamespace("bbotk", quietly = TRUE)) {
    stop("Package 'bbotk' is required. Install with: install.packages('bbotk')")
  }
  
  tryCatch({
    tuner <- mlr3tuning::tnr("mbo")
    tune_res <- surv_train_and_tune(object = object,
                                    learner_id = learner_id,
                                    tuning_budget = tuning_budget,
                                    tuner = tuner,
                                    ...)
    return(tune_res)
  }, error = function(e) {
    warning("Bayesian optimisation failed: ", e$message,
            ". Falling back to random search.")
    surv_train_and_tune(object = object,
                        learner_id = learner_id,
                        tuning_budget = tuning_budget,
                        ...)
  })
}

# ---------------------------------------------------------------------------
# 7. AIC Stepwise Variable Selection (Cox model)
# ---------------------------------------------------------------------------
#' Perform AIC stepwise variable selection on a Cox model
#'
#' @param object A `TaskSurv` or `PrognosiX` object.
#' @param direction One of `"both"`, `"backward"`, `"forward"` (default "both").
#' @param scope List of formulas for upper/lower scope (optional).
#' @param trace Integer; control tracing (default 0).
#' @param ... Additional arguments to `MASS::stepAIC`.
#' @return A list with selected features, the final model, and the stepwise trace.
#' @export
stepwise_variable_selection <- function(object,
                                        direction = c("both", "backward", "forward"),
                                        scope = NULL,
                                        trace = 0,
                                        ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install with: install.packages('MASS')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  direction <- match.arg(direction)
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  time_col <- task$target_names[1]
  status_col <- task$target_names[2]
  
  all_features <- task$feature_names
  formula_str <- paste("survival::Surv(", time_col, ", ", status_col, ") ~",
                       paste(all_features, collapse = " + "))
  full_model <- survival::coxph(as.formula(formula_str), data = data)
  
  step_model <- MASS::stepAIC(full_model,
                              direction = direction,
                              scope = scope,
                              trace = trace,
                              ...)
  
  selected_features <- names(coef(step_model))
  list(selected_features = selected_features,
       model = step_model,
       call = step_model$call,
       anova = step_model$anova)
}

# ---------------------------------------------------------------------------
# 8. Calibration‑in‑the‑Large (overall calibration)
# ---------------------------------------------------------------------------
#' Assess overall calibration: mean predicted vs observed survival probability
#'
#' @param learner A trained mlr3 learner (supports "distr").
#' @param object A `TaskSurv` or `PrognosiX` object.
#' @param time_point Time point for evaluation (default median event time).
#' @return A list with calibration slope, intercept, and a plot.
#' @export
calibration_in_the_large <- function(learner, object,
                                     time_point = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  time <- data[[task$target_names[1]]]
  status <- data[[task$target_names[2]]]
  
  if (is.null(time_point)) {
    event_times <- time[status == 1]
    if (length(event_times) == 0) stop("No events to determine time_point.")
    time_point <- median(event_times)
  }
  
  if (!("distr" %in% learner$predict_types)) {
    stop("Learner must support 'distr' predictions.")
  }
  learner$predict_type <- "distr"
  pred <- learner$predict(task)
  surv_prob <- as.numeric(pred$distr$survival(time_point))
  event_prob <- 1 - surv_prob
  
  km_fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = data)
  km_summary <- summary(km_fit, times = time_point, extend = TRUE)
  obs_prob <- 1 - km_summary$surv[1]
  
  # Compute calibration slope and intercept on logit scale
  logit_pred <- log(event_prob / (1 - event_prob))
  # Use observed status at time point (if time <= time_point & status == 1, else 0)
  obs_binary <- ifelse(time <= time_point & status == 1, 1, 0)
  idx <- time <= time_point
  if (sum(idx) < 10) stop("Too few patients with follow‑up up to time_point.")
  fit <- glm(obs_binary[idx] ~ logit_pred[idx], family = binomial())
  coefs <- coef(fit)
  
  list(intercept = coefs[1],
       slope = coefs[2],
       mean_predicted = mean(event_prob, na.rm = TRUE),
       observed = obs_prob,
       n_at_risk = sum(idx))
}

# ---------------------------------------------------------------------------
# 9. External Validation Test (Survival Hosmer‑Lemeshow)
# ---------------------------------------------------------------------------
#' Perform a survival version of the Hosmer‑Lemeshow test
#'
#' @param learner A trained mlr3 learner (supports "distr").
#' @param object A `TaskSurv` or `PrognosiX` object (validation set).
#' @param time_point Time point for evaluation (default median event time).
#' @param n_groups Number of risk groups for HL test (default 10).
#' @return A list with chi‑square statistic, p‑value, and observed/expected tables.
#' @export
external_validation_test <- function(learner, object,
                                     time_point = NULL,
                                     n_groups = 10) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  time <- data[[task$target_names[1]]]
  status <- data[[task$target_names[2]]]
  
  if (is.null(time_point)) {
    event_times <- time[status == 1]
    if (length(event_times) == 0) stop("No events to determine time_point.")
    time_point <- median(event_times)
  }
  
  if (!("distr" %in% learner$predict_types)) {
    stop("Learner must support 'distr' predictions.")
  }
  learner$predict_type <- "distr"
  pred <- learner$predict(task)
  surv_prob <- as.numeric(pred$distr$survival(time_point))
  event_prob <- 1 - surv_prob
  
  idx <- time >= time_point
  if (sum(idx) < n_groups) stop("Too few patients for grouping.")
  event_prob_sub <- event_prob[idx]
  status_sub <- status[idx]
  time_sub <- time[idx]
  
  groups <- cut(event_prob_sub,
                breaks = quantile(event_prob_sub, probs = seq(0, 1, length.out = n_groups + 1)),
                include.lowest = TRUE)
  obs_events <- tapply(status_sub, groups, sum)
  exp_events <- tapply(event_prob_sub, groups, sum)
  
  chi2 <- sum((obs_events - exp_events)^2 / exp_events, na.rm = TRUE)
  df <- n_groups - 1
  p_val <- pchisq(chi2, df, lower.tail = FALSE)
  
  list(chi_square = chi2,
       df = df,
       p_value = p_val,
       observed = obs_events,
       expected = exp_events,
       n_total = sum(idx))
}

# ---------------------------------------------------------------------------
# 10. Clinical Impact Curve (with bootstrapped CIs)
# ---------------------------------------------------------------------------
#' Plot clinical impact curve with uncertainty bands
#'
#' @param learner A trained mlr3 learner (supports "distr").
#' @param object A `TaskSurv` or `PrognosiX` object.
#' @param thresholds Numeric vector (default seq(0.05, 0.95, 0.05)).
#' @param time_point Time point for event probability (default median event time).
#' @param n_boot Bootstrap iterations for confidence bands (default 100).
#' @param save_plot Logical; save to `sub_dir("05_xai")`.
#' @return A `ggplot` object.
#' @export
plot_clinical_impact <- function(learner, object,
                                 thresholds = seq(0.05, 0.95, by = 0.05),
                                 time_point = NULL,
                                 n_boot = 100,
                                 save_plot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  time <- data[[task$target_names[1]]]
  status <- data[[task$target_names[2]]]
  
  if (is.null(time_point)) {
    event_times <- time[status == 1]
    if (length(event_times) == 0) stop("No events to determine time_point.")
    time_point <- median(event_times)
  }
  
  if (!("distr" %in% learner$predict_types)) {
    stop("Learner must support 'distr' predictions.")
  }
  learner$predict_type <- "distr"
  pred <- learner$predict(task)
  surv_prob <- as.numeric(pred$distr$survival(time_point))
  event_prob <- 1 - surv_prob
  
  n_high_risk_boot <- matrix(NA, nrow = length(thresholds), ncol = n_boot)
  n_events_boot <- matrix(NA, nrow = length(thresholds), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    idx <- sample(length(event_prob), replace = TRUE)
    ep_boot <- event_prob[idx]
    status_boot <- status[idx]
    for (i in seq_along(thresholds)) {
      high_idx <- ep_boot >= thresholds[i]
      n_high_risk_boot[i, b] <- sum(high_idx)
      n_events_boot[i, b] <- sum(status_boot[high_idx] == 1)
    }
  }
  
  impact_df <- data.frame(
    threshold = thresholds,
    n_high_risk = rowMeans(n_high_risk_boot, na.rm = TRUE),
    n_high_risk_lower = apply(n_high_risk_boot, 1, quantile, 0.025, na.rm = TRUE),
    n_high_risk_upper = apply(n_high_risk_boot, 1, quantile, 0.975, na.rm = TRUE),
    n_events_high = rowMeans(n_events_boot, na.rm = TRUE),
    n_events_lower = apply(n_events_boot, 1, quantile, 0.025, na.rm = TRUE),
    n_events_upper = apply(n_events_boot, 1, quantile, 0.975, na.rm = TRUE)
  )
  
  p <- ggplot2::ggplot(impact_df, ggplot2::aes(x = threshold)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = n_high_risk_lower,
                                      ymax = n_high_risk_upper,
                                      fill = "High‑risk patients"),
                         alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = n_high_risk, colour = "High‑risk patients"),
                       linewidth = 1.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = n_events_lower,
                                      ymax = n_events_upper,
                                      fill = "Events in high‑risk"),
                         alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = n_events_high, colour = "Events in high‑risk"),
                       linewidth = 1.2, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = c("High‑risk patients" = "#1f77b4",
                                          "Events in high‑risk" = "#d62728"),
                               name = NULL) +
    ggplot2::scale_colour_manual(values = c("High‑risk patients" = "#1f77b4",
                                            "Events in high‑risk" = "#d62728"),
                                 name = NULL) +
    ggplot2::labs(x = "Risk Threshold",
                  y = "Number of Patients",
                  title = "Clinical Impact Curve (with 95% CI)",
                  subtitle = paste("Evaluated at t =", time_point)) +
    ggprism::theme_prism()
  
  print(p)
  
  if (save_plot) {
    if (exists("sub_dir", mode = "function")) {
      save_path <- sub_dir("05_xai")
    } else {
      save_path <- "./05_xai"
      dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    }
    ggplot2::ggsave(file.path(save_path, "clinical_impact.pdf"),
                    p, width = 8, height = 5)
    message("Clinical impact plot saved to: ", save_path)
  }
  
  invisible(p)
}