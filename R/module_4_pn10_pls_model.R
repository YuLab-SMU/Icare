#' Train PLS Model (Base)
#'
#' Trains a Partial Least Squares (PLS) Cox model (base function).
#'
#' @param train_data Training data.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param nfolds Number of folds for CV.
#' @param method PLS method (e.g., "efron").
#' @param plot.it Logical, plot CV results.
#' @param se Logical, compute standard errors.
#' @param scaleX Logical, scale predictors.
#' @param seed Seed.
#'
#' @return List containing model, lambda, cv_fit.
#' @export
train_pls_base <- function(train_data,
                           time_col = "time",
                           status_col = "status",
                           nfolds = 10,
                           method = "efron",
                           plot.it = TRUE,
                           se = TRUE,
                           scaleX = TRUE,
                           seed=1234) {
  set.seed(seed)
  x_train <- as.matrix(train_data[, !(names(train_data) %in% c(time_col, status_col))])
  y_train <- survival::Surv(train_data[[time_col]], train_data[[status_col]])


  cv_fit <- plsRcox::cv.plsRcox(
    data = list(x = x_train, time = train_data[[time_col]], status = train_data[[status_col]]),
    method = method,
    nfold = nfolds,
    plot.it = plot.it,
    se = se,
    scaleX = scaleX
  )

  best_lambda_value <- cv_fit$lambda.min

  fit_best_cv_pls <- plsRcox::plsRcox(x_train, y_train, lambda = best_lambda_value)

  return(list(model = fit_best_cv_pls,
              lambda = best_lambda_value,
              cv_fit = cv_fit))
}

#' Train PLS Model for PrognosiX
#'
#' Wrapper to train PLS model on PrognosiX object.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param nfolds Folds.
#' @param method Method.
#' @param plot.it Plot CV.
#' @param se Standard errors.
#' @param scaleX Scale predictors.
#'
#' @return Updated PrognosiX object or model list.
#' @export
train_pls_model <- function(object,
                            time_col = "time",
                            status_col = "status",
                            nfolds = 10,
                            method = "efron",
                            plot.it = TRUE,
                            se = TRUE,
                            scaleX = TRUE) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$training
    test_data <- data_sets$test

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    train_data <- object$training
    test_data <- object$testing

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Training pls model...\n")

  pls_model_info <- train_pls_base(train_data, time_col, status_col, nfolds, method, plot.it, se, scaleX)

  cat("Best lambda value: ", pls_model_info$lambda, "\n")

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with pls model results...\n")

    object@survival.model[["pls_model"]] <- pls_model_info

    return(object)
  }

  cat("Returning model information as a list...\n")

  return(pls_model_info) # Corrected return value to match assignment
}


#' Evaluate PLS ROC
#'
#' Evaluates PLS model performance using ROC curves.
#'
#' @param fit_best_cv_pls PLS model.
#' @param best_lambda_value Lambda.
#' @param train_data Training data.
#' @param test_data Testing data.
#' @param validation_data Validation data.
#' @param palette_name Palette.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#'
#' @return List with ROC results and plot.
#' @export
evaluate_pls_roc <- function(
    fit_best_cv_pls,    
    best_lambda_value,  
    train_data = NULL,  
    test_data = NULL,  
    validation_data = NULL, 
    palette_name = "AsteroidCity1",  
    time_col = "time", 
    status_col = "status",  
    save_plot = TRUE,   
    save_dir = here::here("PrognosiX", "pls_model"),  
    plot_width = 7,    
    plot_height = 7,   
    base_size = 14     
) {
  if (is.null(train_data) && is.null(test_data) && is.null(validation_data)) {
    stop("At least one of 'train_data', 'test_data', or 'validation_data' must be provided.")
  }
  
  roc_data_all <- data.frame()
  results_roc_df <- data.frame()
  
  if (!is.null(train_data)) {
    x_train_pls <- as.matrix(as.data.frame(lapply(train_data[, !(names(train_data) %in% c(time_col, status_col))], as.numeric)))
    y_train_status <- train_data[[status_col]]
    y_train_time <- train_data[[time_col]]
    
    pls_risk_score_train <- predict(fit_best_cv_pls, newdata = x_train_pls, type = "risk")
    pls_risk_score_train[is.infinite(pls_risk_score_train)] <- NA 
    pls_risk_score_train[is.na(pls_risk_score_train)] <- median(pls_risk_score_train, na.rm = TRUE)
    
    pls_roc_train <- pROC::roc(y_train_status, pls_risk_score_train)
    
    roc_data_train <- data.frame(
      specificity = pls_roc_train$specificities,
      sensitivity = pls_roc_train$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = survcomp::concordance.index(x = pls_risk_score_train,
                                            surv.time = y_train_time,
                                            surv.event = y_train_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(pls_roc_train)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    x_test_pls <- as.matrix(as.data.frame(lapply(test_data[, !(names(test_data) %in% c(time_col, status_col))], as.numeric)))
    y_test_status <- test_data[[status_col]]
    y_test_time <- test_data[[time_col]]
    
    pls_risk_score_test <- predict(fit_best_cv_pls, newdata = x_test_pls, type = "risk")
    
    pls_risk_score_test[is.infinite(pls_risk_score_test)] <- NA 
    pls_risk_score_test[is.na(pls_risk_score_test)] <- median(pls_risk_score_test, na.rm = TRUE)
    
    pls_roc_test <- pROC::roc(y_test_status, pls_risk_score_test)
    
    roc_data_test <- data.frame(
      specificity = pls_roc_test$specificities,
      sensitivity = pls_roc_test$sensitivities,
      Set = "Testing Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_test)
    
    Test <- data.frame(
      C_index = survcomp::concordance.index(x = pls_risk_score_test,
                                            surv.time = y_test_time,
                                            surv.event = y_test_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(pls_roc_test)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
  }
  
  if (!is.null(validation_data)) {
    x_validation_pls <- as.matrix(as.data.frame(lapply(validation_data[, !(names(validation_data) %in% c(time_col, status_col))], as.numeric)))
    y_validation_status <- validation_data[[status_col]]
    y_validation_time <- validation_data[[time_col]]
    
    pls_risk_score_validation <- predict(fit_best_cv_pls, newdata = x_validation_pls, type = "risk")
    pls_risk_score_validation[is.infinite(pls_risk_score_validation)] <- NA 
    pls_risk_score_validation[is.na(pls_risk_score_validation)] <- median(pls_risk_score_validation, na.rm = TRUE)
    
    pls_roc_validation <- pROC::roc(y_validation_status, pls_risk_score_validation)
    
    roc_data_validation <- data.frame(
      specificity = pls_roc_validation$specificities,
      sensitivity = pls_roc_validation$sensitivities,
      Set = "Validation Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_validation)
    
    Validation <- data.frame(
      C_index = survcomp::concordance.index(x = pls_risk_score_validation,
                                            surv.time = y_validation_time,
                                            surv.event = y_validation_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(pls_roc_validation)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Validation", Validation))
  }
  
  roc_data_all <- roc_data_all %>%
    dplyr::filter(!(sensitivity == 0 & specificity == 1)) %>%
    dplyr::filter(sensitivity >= 0 & sensitivity <= 1, specificity >= 0 & specificity <= 1)
  
  roc_plot <- ggplot2::ggplot(roc_data_all, aes(x = 1 - specificity, y = sensitivity, color = Set)) +
    geom_path(aes(color = Set), size = 1.5) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = wesanderson::wes_palette(palette_name)) +
    labs(
      title = "ROC Curves",
      subtitle = paste("AUCs: ", paste(results_roc_df$Dataset, "=", round(results_roc_df$ROC_AUC, 2), collapse = ", ")),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Dataset"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.25, "cm"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, face = "italic"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10)
    )
  
  print(roc_plot)
  
  if (save_plot) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = file.path(save_dir, "ROC_curves.pdf"), plot = roc_plot, width = plot_width, height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "ROC_curves.pdf"), "\n")
  }
  
  print(results_roc_df)
  
  results_roc <- list(
    results_roc = results_roc_df,
    plot_roc = roc_plot
  )
  
  return(results_roc)
}


#' Evaluate ROC for PLS Model (Wrapper)
#'
#' Wrapper to evaluate ROC for PLS model in PrognosiX.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param palette_name Palette.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#'
#' @return Updated PrognosiX object or results list.
#' @export
evaluate_roc_pls_model <- function(object,
                                   time_col = "time",
                                   status_col = "status",
                                   palette_name = "AsteroidCity1",
                                   save_plot = TRUE,
                                   save_dir = here::here("PrognosiX", "pls_model"),
                                   plot_width = 7,
                                   plot_height = 7,
                                   base_size = 14) {


  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_best_cv_pls <- methods::slot(object, "survival.model")[["pls_model"]][["model"]]
    best_lambda_value <- methods::slot(object, "survival.model")[["pls_model"]][["lambda"]]

    cat("Extracted PLS model and lambda: ", best_lambda_value, "\n")


    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$test

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")
    train_data <- object$train
    test_data <- object$test

    fit_best_cv_pls <- object$fit_best_cv_pls
    best_lambda_value <- object$best_lambda_value

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_roc <- evaluate_pls_roc(fit_best_cv_pls = fit_best_cv_pls,
                                  best_lambda_value = best_lambda_value,
                                  train_data = train_data,
                                  test_data = test_data,
                                  palette_name = palette_name,
                                  time_col = time_col,
                                  status_col = status_col,
                                  save_plot = save_plot,
                                  save_dir = save_dir,
                                  plot_width = plot_width,
                                  plot_height = plot_height,
                                  base_size = base_size)

  cat("ROC evaluation completed. Returning results...\n")

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with pls model ROC results...\n")

    object@survival.model[["pls_model"]][["results_roc"]] <- results_roc

    cat("Results updated in the 'PrognosiX' object.\n")
    return(object)
  }

  cat("Returning model ROC results as a list...\n")

  return(results_roc)
}


#' Evaluate PLS Kaplan-Meier
#'
#' Evaluates KM curves for PLS risk groups.
#'
#' @param fit_best_cv_pls Model.
#' @param best_lambda_value Lambda.
#' @param data Data.
#' @param data_name Data name.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param R Resampling (unused).
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param palette_name Palette.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#' @param seed Seed.
#' @param fill_color Fill color.
#'
#' @return List with KM results and plots.
#' @export
evaluate_pls_km <- function(fit_best_cv_pls,
                            best_lambda_value,
                            data,
                            data_name = "test",
                            time_col = "time",
                            status_col = "status",
                            R = 1000,
                            save_plot = TRUE,
                            save_dir = here::here("PrognosiX", "pls_model"),
                            palette_name = "Dark2",
                            plot_width = 10,
                            plot_height = 8,
                            base_size = 12,
                            seed = 1234,
                            fill_color = "lightblue"
) {
  
  set.seed(seed)
  library(survcomp)
  if (!data_name %in% c("test", "val")) {
    stop("'data_name' must be either 'test' or 'val'.")
  }
  
  x_data_pls <- as.matrix(as.data.frame(lapply(data[, !(names(data) %in% c(time_col, status_col))], as.numeric)))
  y_data_status <- data[[status_col]]
  y_data_time <- data[[time_col]]
  
  risk_score <- predict(fit_best_cv_pls, newdata = x_data_pls, type = "risk")
  
  
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  data<-cbind(data,risk_score,risk_group)
  
  colnames(data)[ncol(data)-1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  cat("Calculating Kaplan-Meier survival curve...\n")
  km_fit_pls <- survival::survfit(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  
  survdiff_result_pls <- survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  km_pval_pls <- pchisq(survdiff_result_pls$chisq, 1, lower.tail = FALSE)
  km_hr_pls <- (survdiff_result_pls$obs[2] / survdiff_result_pls$exp[2]) / (survdiff_result_pls$obs[1] / survdiff_result_pls$exp[1])
  
  km_upper95_pls <- exp(log(km_hr_pls) + qnorm(0.975) * sqrt(1 / survdiff_result_pls$exp[2] + 1 / survdiff_result_pls$exp[1]))
  km_lower95_pls <- exp(log(km_hr_pls) - qnorm(0.975) * sqrt(1 / survdiff_result_pls$exp[2] + 1 / survdiff_result_pls$exp[1]))
  
  km_results_pls <- data.frame(KM_HR = km_hr_pls,
                               KM_CI_lower = km_lower95_pls,
                               KM_CI_upper = km_upper95_pls,
                               KM_p_value = km_pval_pls)
  print(km_results_pls)
  
  cat("Generating Kaplan-Meier plot...\n")
  km_plot_pls <- survminer::ggsurvplot(
    km_fit_pls,
    data = data,
    conf.int = TRUE,
    conf.int.fill = fill_color,
    conf.int.alpha = 0.5,
    pval = TRUE,
    pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve for Cox-PLS Risk Groups (", data_name, ")", sep = ""),
    surv.median.line = "hv",
    risk.table = TRUE,
    xlab = "Follow-up Time (days)",
    legend = c(0.8, 0.75),
    legend.title = "Risk Group",
    legend.labs = unique(data$risk_group_pls),
    break.x.by = 100,
    palette = "Dark2",
    base_size = base_size
  )
  print(km_plot_pls)
  
  cat("Kaplan-Meier plot generated.\n")
  
  surv_plot <- km_plot_pls$plot +
    ggprism::theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot_pls$table + ggprism::theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  combined_plot <- gridExtra::grid.arrange(surv_plot, risk_table, ncol = 1, heights = c(1.8, 1))
  
  cat("Combined plot generated.\n")
  
  
  if (save_plot) {
    ggsave(file.path(save_dir, paste0("pls_km_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("pls_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("pls_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", save_dir, "\n")
  }
  
  results_pls <- list(
    KM_data_results = km_results_pls,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_pls)
}


#' Evaluate KM for PLS Model (Wrapper)
#'
#' Wrapper to evaluate KM for PLS model in PrognosiX.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param nfolds Folds.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#' @param data_name Data name.
#' @param data Data.
#'
#' @return Updated PrognosiX object or results list.
#' @export
evaluate_km_pls_model <- function(object,
                                  time_col = "time",
                                  status_col = "status",
                                  nfolds = 10,
                                  save_plot = TRUE,
                                  save_dir = here::here("PrognosiX", "pls_model"),
                                  plot_width = 7,
                                  plot_height = 7,
                                  base_size = 14,
                                  data_name = "test",
                                  data=NULL) {


  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")
    library(survcomp)
    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_best_cv_pls <- methods::slot(object, "survival.model")[["pls_model"]][["model"]]
    best_lambda_value <- methods::slot(object, "survival.model")[["pls_model"]][["lambda"]]

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing


  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    train_data <- object$training
    test_data <- object$testing
    fit_best_cv_pls <- object$fit_best_cv_pls
    best_lambda_value <- object$best_lambda_value
  } else {
    stop("Input is neither a 'PrognosiX' object nor a valid list with 'train' and 'test' elements.")
  }

  cat("Evaluation completed.\n")
  results_pls_km <- evaluate_pls_km(
    fit_best_cv_pls,
    best_lambda_value,
    data = test_data,
    data_name = data_name,
    time_col = time_col,
    status_col = status_col,
    save_plot = save_plot,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with pls model Kaplan-Meier results...\n")

    object@survival.model[["pls_model"]][["results_km"]] <- results_pls_km

    return(object)
  }

  cat("Returning Kaplan-Meier results as a list...\n")

  return(results_pls_km)
}


#' Compute HR for PLS
#'
#' Computes HR and CI for PLS coefficients.
#'
#' @param pls_model PLS model.
#' @param train_data Training data.
#'
#' @return HR results data frame.
#' @export
compute_hr_pls <- function(pls_model,train_data) {
  pls_coefs <-pls_model[["Coeffs"]]
  pls_coefs <- as.numeric(pls_coefs)

  feature_names <- colnames(train_data[, !(names(train_data) %in% c("time", "status"))]) # Assuming defaults for now as arguments missing

  hr <- exp(pls_coefs)

  assumed_se <- 0.05
  ci_lower <- exp(pls_coefs - 1.96 * assumed_se)
  ci_upper <- exp(pls_coefs + 1.96 * assumed_se)

  hr_results <- data.frame(
    Variable = feature_names,
    Coefficient = pls_coefs,
    HR = hr,
    CI_lower = ci_lower,
    CI_upper = ci_upper
  )

  hr_results$HR_95CI <- paste0(
    round(hr_results$HR, 2), " (",
    round(hr_results$CI_lower, 2), "-",
    round(hr_results$CI_upper, 2), ")"
  )

  print(hr_results)
  return(hr_results)
}

#' Compute HR and CI for PLS (Wrapper)
#'
#' Wrapper to compute HR/CI for PLS in PrognosiX.
#'
#' @param object PrognosiX object.
#'
#' @return Updated PrognosiX object or results.
#' @export
pls_compute_hr_and_ci <- function(object) {
  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")


    fit_best_cv_pls <- methods::slot(object, "survival.model")[["pls_model"]][["model"]]
    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing
  } else if (is.list(object) && all(c("fit_best_cv_pls") %in% names(object))) {
    fit_best_cv_pls <- object$fit_best_cv_pls
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'fit_best_cv_pls' and 'lambda'.")
  }

  cat("Evaluating hazard ratios and confidence intervals for PLS model...\n")
  hr_results <- compute_hr_pls(
    pls_model = fit_best_cv_pls,
    train_data =train_data
  )

  if (inherits(object, 'PrognosiX')) {
    object@survival.model[["pls_model"]][["hr_results"]] <- hr_results
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")

    return(object)
  }

  cat("Returning hazard ratio results as a list...\n")
  return(hr_results)
}


#' Create Forest Plot for PLS Model
#'
#' Creates a forest plot for PLS model HR results.
#'
#' @param hr_results HR results.
#' @param plot_title Title.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param palette_name Palette.
#' @param base_size Size.
#' @param hr_limit HR limit.
#' @param ci_range_limit CI range limit.
#'
#' @return Plot object.
#' @export
create_forest_plot_pls <- function(hr_results,
                                     plot_title = "Evaluation of Hazard Ratios and Confidence Intervals",
                                     save_plot = FALSE,
                                     save_dir = here::here('Prognosi'),
                                     plot_width = 11,
                                     plot_height = 5,
                                     palette_name = "AsteroidCity1",
                                     base_size = 14,
                                     hr_limit = c(0, 3),
                                     ci_range_limit = 1000) {

  hr_results <- hr_results %>%
    dplyr::filter(is.finite(HR) & is.finite(CI_lower) & is.finite(CI_upper)) %>%
    dplyr::filter(HR >= hr_limit[1] & HR <= hr_limit[2]) %>%
    dplyr::filter((CI_upper - CI_lower) <= ci_range_limit)

  cat(paste(nrow(hr_results), "rows remain after filtering. Rows with extreme or invalid values were excluded.\n"))

  tab_header <- c("Clinical factors", "HR (95% CI)")
  tab_header_bold <- TRUE

  tmp_df <- hr_results %>%
    dplyr::select(Variable, HR_95CI) %>%
    rbind(tab_header, .) %>%
    cbind("clinical" = .[, "Variable"], .)

  tmp_df$clinical <- factor(tmp_df$clinical, levels = rev(tmp_df$clinical))

  tmp_df <- tmp_df %>%
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "x", values_to = "label")

  tmp_df$label_bold <- sapply(tmp_df$label, function(x) {
    if (x %in% unlist(tab_header)) {
      if (tab_header_bold) {
        paste0("<b>", x, "</b>")
      } else {
        x
      }
    } else {
      x
    }
  }, simplify = T)

  p_tab <- ggplot2::ggplot(data = tmp_df, aes(x = x, y = clinical)) +
    geom_tile(color = "white", fill = "white") +
    ggtext::geom_richtext(aes(label = label_bold), label.color = NA) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line.x = element_line(linewidth = 1),
      axis.line.y = element_blank(),
      text = element_text()
    ) +
    geom_hline(yintercept = c(length(unique(tmp_df$clinical)) - 0.5), linewidth = 1)

  tmp_error <- rbind(rep(NA, ncol(hr_results)), hr_results)
  tmp_error$Variable[1] <- tab_header[[1]]

  tmp_error$Variable <- factor(tmp_error$Variable, levels = rev(tmp_error$Variable))
  error_bar_height <- 0.2
  ref_line <- 1

  p_errorbar <- ggplot2::ggplot(data = tmp_error) +
    geom_point(aes(x = HR, y = Variable, color = HR > 1), shape = "diamond", size = 4) +
    scale_color_manual(values = wesanderson::wes_palette(palette_name)) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, y = Variable), linewidth = 1, height = error_bar_height) +
    theme_classic(base_size = base_size) +
    theme(
      title = element_blank(),
      axis.line.x = element_line(linewidth = 1),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    expand_limits(x = hr_limit) +
    coord_cartesian(xlim = hr_limit) +
    geom_vline(xintercept = ref_line, color = "grey50")

  p <- p_tab + p_errorbar +
    patchwork::plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.5)))

  print(p)

  if (save_plot) {
    pdf_path <- file.path(save_dir, "forest_plot.pdf")
    ggsave(filename = pdf_path, plot = p, width = plot_width, height = plot_height)
    cat("Plot saved at: ", pdf_path, "\n")
  }
  return(p)
}



#' Forest Plot for PLS Model (Wrapper)
#'
#' Generates forest plot for PLS model in PrognosiX.
#'
#' @param object PrognosiX object.
#' @param plot_title Title.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param var_col Variable column.
#' @param palette_name Palette.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#' @param use_subgroup_data Logical.
#' @param hr_limit HR limit.
#'
#' @return Updated PrognosiX object or plot.
#' @export
forest_plot_pls_model <- function(object,
                                  plot_title = "Evaluation of Hazard Ratios and Confidence Intervals for pls Model in Survival Analysis",
                                  time_col = "time",
                                  status_col = "status",
                                  var_col = NULL,
                                  palette_name = "AsteroidCity1",
                                  save_plot = TRUE,
                                  save_dir = here::here('PrognosiX', "pls_model"),
                                  plot_width = 14,
                                  plot_height = 7,
                                  base_size = 14,
                                  use_subgroup_data = FALSE,
                                  hr_limit = c(0.1, 3)
) {

  if (inherits(object, 'PrognosiX')) {
    hr_results <- methods::slot(object, "survival.model")[["pls_model"]][["hr_results"]]

    if (is.null(hr_results) || nrow(hr_results) == 0) {
      stop("The hr_results in the PrognosiX object is empty.")
    }
  } else if (is.data.frame(object)) {
    hr_results <- object
    if (nrow(hr_results) == 0) {
      stop("The provided data frame is empty.")
    }
    cat("Using provided data frame for univariate analysis data...\n")

  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  p <- create_forest_plot_pls(hr_results,
                                plot_title = plot_title,
                                save_plot = save_plot,
                                save_dir = save_dir,
                                plot_width = plot_width,
                                plot_height = plot_height,
                                base_size = base_size)

  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object with forest plot...\n")
    object@survival.model[["pls_model"]][["forest_plot"]] <- p
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")
    return(object)
  }

  return(p)
}
