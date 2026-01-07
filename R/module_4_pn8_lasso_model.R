#' Train LASSO Model (Base)
#'
#' Trains a LASSO Cox model (base function).
#'
#' @param train_data Training data.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param nfolds Number of folds for CV.
#' @param palette_name Palette name.
#' @param save_plot Logical.
#' @param save_dir Save directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base font size.
#' @param seed Seed.
#'
#' @return List containing model, lambda, cv_fit, and plot.
#' @export
train_lasso_base <- function(train_data,
                             time_col = "time",
                             status_col = "status",
                             nfolds = 10,
                             palette_name = "AsteroidCity1",
                             save_plot = TRUE,
                             save_dir = here::here("PrognosiX", "lasso_model"),
                             plot_width = 5,
                             plot_height = 5,
                             base_size = 14,
                             seed = 123) {

  set.seed(seed)

  x_train <- as.matrix(train_data[, !(names(train_data) %in% c(time_col, status_col))])
  y_train <- survival::Surv(train_data[[time_col]], train_data[[status_col]])

  cv_fit <- glmnet::cv.glmnet(x_train, y_train, family = "cox", alpha = 1, nfolds = nfolds)
  best_lambda_value <- cv_fit$lambda.min

  fit_best_cv_lasso <- glmnet::glmnet(x_train, y_train, family = "cox", alpha = 1, lambda = best_lambda_value)

  palette <-  wesanderson::wes_palette(palette_name, type = "discrete")

  cv_data <- data.frame(
    lambda = cv_fit$lambda,
    mean_error = cv_fit$cvm,
    se_error = cv_fit$cvsd
  )

  p <- ggplot2::ggplot(cv_data, aes(x = lambda, y = mean_error)) +
    geom_line(color = palette[1]) +
    geom_ribbon(aes(ymin = mean_error - se_error, ymax = mean_error + se_error),
                fill = palette[2], alpha = 0.3) +
    scale_x_log10() +
    labs(x = "Lambda", y = "Mean Cross-Validation Error",
         title = "Cross-Validation for Lasso model") +
    geom_vline(xintercept = cv_fit$lambda.min, color = palette[3], linetype = "dashed") +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  print(p)

  if (save_plot) {
    plot_filename <- file.path(save_dir, paste0("cv_lasso_plot", ".pdf"))
    ggsave(plot_filename, plot = p, width = plot_width, height = plot_height,
           device = "pdf")
    cat("Plot saved to:", plot_filename, "\n")
  }

  return(list(model = fit_best_cv_lasso,
              lambda = best_lambda_value,
              cv_fit = cv_fit,
              cv_plot = p))
}


#' Train LASSO Model for PrognosiX
#'
#' Wrapper to train LASSO model on PrognosiX object.
#'
#' @param object PrognosiX object or list.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param nfolds Folds.
#' @param palette_name Palette.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#' @param seed Seed.
#'
#' @return Updated PrognosiX object or model list.
#' @export
train_lasso_model <- function(object,
                              time_col = "time",
                              status_col = "status",
                              nfolds = 10,
                              palette_name = "AsteroidCity1",
                              save_plot = TRUE,
                              save_dir = here::here("PrognosiX", "lasso_model"),
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 14,
                              seed = 123) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    train_data <- object$training
    test_data <- object$testing

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Training lasso model...\n")

  lasso_model_info <- train_lasso_base(train_data,
                                       time_col = time_col,
                                       status_col = status_col,
                                       nfolds = nfolds,
                                       palette_name = palette_name,
                                       save_plot = save_plot,
                                       save_dir = save_dir,
                                       plot_width = plot_width,
                                       plot_height = plot_height,
                                       base_size = base_size,
                                       seed = seed)

  cat("Best lambda value: ", lasso_model_info$lambda, "\n")

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with lasso model results...\n")

    object@survival.model[["lasso_model"]] <- lasso_model_info

    return(object)
  }

  cat("Returning model information as a list...\n")

  return(lasso_model_info)
}

#' Evaluate LASSO ROC
#'
#' Evaluates LASSO model performance using ROC curves.
#'
#' @param fit_best_cv_lasso Fitted LASSO model.
#' @param best_lambda_value Best lambda.
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
evaluate_lasso_roc <- function(
    fit_best_cv_lasso,
    best_lambda_value,
    train_data = NULL,
    test_data = NULL,
    validation_data = NULL,
    palette_name = "AsteroidCity1",
    time_col = "time",
    status_col = "status",
    save_plot = TRUE,
    save_dir = here::here("PrognosiX", "lasso_model"),
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
    x_train <- as.matrix(train_data[, !(names(train_data) %in% c(time_col, status_col))])
    lasso_risk_score_train <- predict(fit_best_cv_lasso,
                                      newx = x_train,
                                      s = best_lambda_value,
                                      type = "link")
    lasso_risk_score_train[is.infinite(lasso_risk_score_train)] <- NA 
    lasso_risk_score_train[is.na(lasso_risk_score_train)] <- median(lasso_risk_score_train, na.rm = TRUE)
    
    lasso_roc_train <- pROC::roc(train_data[[status_col]], as.numeric(lasso_risk_score_train))
 
    
    roc_data_train <- data.frame(
      specificity = lasso_roc_train$specificities,
      sensitivity = lasso_roc_train$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = survcomp::concordance.index(x = lasso_risk_score_train,
                                  surv.time = train_data[[time_col]],
                                  surv.event = train_data[[status_col]],
                                  method = "noether")$c.index,
      ROC_AUC = pROC::auc(lasso_roc_train)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    x_test <- as.matrix(test_data[, !(names(test_data) %in% c(time_col, status_col))])
    lasso_risk_score_test <- predict(fit_best_cv_lasso,
                          newx = x_test,
                          s = best_lambda_value,
                          type = "link")
    lasso_risk_score_test[is.infinite(lasso_risk_score_test)] <- NA 
    lasso_risk_score_test[is.na(lasso_risk_score_test)] <- median(lasso_risk_score_test, na.rm = TRUE)
    
    lasso_roc_test <- pROC::roc(test_data[[status_col]], as.numeric(lasso_risk_score_test))
    
    roc_data_test <- data.frame(
      specificity = lasso_roc_test$specificities,
      sensitivity = lasso_roc_test$sensitivities,
      Set = "Testing Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_test)
    
    Test <- data.frame(
      C_index = survcomp::concordance.index(x = lasso_risk_score_test,
                                  surv.time = test_data[[time_col]],
                                  surv.event = test_data[[status_col]],
                                  method = "noether")$c.index,
      ROC_AUC = pROC::auc(lasso_roc_test)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
  }
  
  if (!is.null(validation_data)) {
    x_validation <- as.matrix(validation_data[, !(names(validation_data) %in% c(time_col, status_col))])
    lasso_risk_score_validation <- predict(fit_best_cv_lasso,
                                           newx = x_validation,
                                           s = best_lambda_value,
                                           type = "link")
    
    lasso_risk_score_validation[is.infinite(lasso_risk_score_validation)] <- NA 
    lasso_risk_score_validation[is.na(lasso_risk_score_validation)] <- median(lasso_risk_score_validation, na.rm = TRUE)
    
    lasso_roc_validation <- pROC::roc(validation_data[[status_col]], as.numeric(lasso_risk_score_validation))
    
    roc_data_validation <- data.frame(
      specificity = lasso_roc_validation$specificities,
      sensitivity = lasso_roc_validation$sensitivities,
      Set = "Validation Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_validation)
    
    Validation <- data.frame(
      C_index = survcomp::concordance.index(x = lasso_risk_score_validation,
                                  surv.time = validation_data[[time_col]],
                                  surv.event = validation_data[[status_col]],
                                  method = "noether")$c.index,
      ROC_AUC = pROC::auc(lasso_roc_validation)
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


#' Evaluate ROC for LASSO Model (Wrapper)
#'
#' Wrapper to evaluate ROC for LASSO model in PrognosiX object.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param nfolds Folds (unused).
#' @param palette_name Palette.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#'
#' @return Updated PrognosiX object or results list.
#' @export
evaluate_roc_lasso_model <- function(object,
                                     time_col = "time",
                                     status_col = "status",
                                     nfolds = 10,
                                     palette_name = "AsteroidCity1",
                                     save_plot = TRUE,
                                     save_dir = here::here("PrognosiX", "lasso_model"),
                                     plot_width = 7,
                                     plot_height = 7,
                                     base_size = 14) {

  if (inherits(object, 'PrognosiX')) {
    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_best_cv_lasso <- methods::slot(object, "survival.model")[["lasso_model"]][["model"]]
    best_lambda_value <- methods::slot(object, "survival.model")[["lasso_model"]][["lambda"]]

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    train_data <- object$training
    test_data <- object$testing

    fit_best_cv_lasso <- object$fit_best_cv_lasso
    best_lambda_value <- object$best_lambda_value

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_roc <- evaluate_lasso_roc(fit_best_cv_lasso, best_lambda_value,
                                    train_data, test_data,
                                    palette_name = palette_name,
                                    time_col = time_col, status_col = status_col,
                                    save_plot = save_plot, save_dir = save_dir,
                                    plot_width = plot_width, plot_height = plot_height,
                                    base_size = base_size)

  if (inherits(object, 'PrognosiX')) {
    object@survival.model[["lasso_model"]][["results_roc"]] <- results_roc
    return(object)
  }

  return(results_roc)
}


#' Evaluate LASSO Kaplan-Meier
#'
#' Evaluates KM curves for LASSO risk groups.
#'
#' @param fit_best_cv_lasso Model.
#' @param best_lambda_value Lambda.
#' @param data Data.
#' @param data_name Data name ("test", "val").
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
evaluate_lasso_km <- function(
    fit_best_cv_lasso,
    best_lambda_value,
    data,
    data_name = "test",
    time_col = "time",
    status_col = "status",
    R = 1000,
    save_plot = TRUE,
    save_dir = here::here("PrognosiX", "lasso_model"),
    palette_name = "Dark2",
    plot_width = 10,
    plot_height = 8,
    base_size = 14,
    seed = 1234,
    fill_color = "lightblue"
) {
  set.seed(seed)
  
  if (!data_name %in% c("test", "val")) {
    stop("'data_name' must be either 'test' or 'val'.")
  }
  
  cat("Predicting lasso risk scores for the", data_name, "data...\n")
  risk_score <- predict(fit_best_cv_lasso,
                        newx = as.matrix(data[, !(names(data) %in% c(time_col, status_col))]),
                        s = best_lambda_value,
                        type = "response")
  
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  cat("Risk groups assigned (High/Low) based on lasso risk scores.\n")
  
  data <- cbind(data, risk_score, risk_group)
  colnames(data)[ncol(data) - 1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  fit_km_lasso <- survival::survfit(survival::Surv(data[[time_col]], data[[status_col]]) ~ risk_group, data = data)
  cat("Kaplan-Meier fit for lasso risk groups completed.\n")
  
  km_pval_lasso <- pchisq(survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                                   data = data)$chisq, 1, lower.tail = FALSE)
  
  km_hr_lasso <- (survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                           data = data)$obs[2] / survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                                                          data = data)$exp[2]) /
    (survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
              data = data)$obs[1] / survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                                             data = data)$exp[1])
  
  km_upper95_lasso <- exp(log(km_hr_lasso) + qnorm(0.975) * sqrt(1 / survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                                                                              data = data)$exp[2] + 1 / survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                                                                                                                 data = data)$exp[1]))
  km_lower95_lasso <- exp(log(km_hr_lasso) - qnorm(0.975) * sqrt(1 / survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                                                                              data = data)$exp[2] + 1 / survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group,
                                                                                                                 data = data)$exp[1]))
  
  km_results_lasso <- data.frame(
    KM_HR = km_hr_lasso,
    KM_CI_lower = km_lower95_lasso,
    KM_CI_upper = km_upper95_lasso,
    KM_p_value = km_pval_lasso
  )
  cat("Kaplan-Meier results computed: HR, CI, and p-value.\n")
  
  # Plot Kaplan-Meier Curve
  km_plot_lasso <- survminer::ggsurvplot(
    fit_km_lasso,
    data = data,
    conf.int = TRUE,
    conf.int.fill = fill_color,
    conf.int.alpha = 0.5,
    pval = TRUE,
    pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve for Cox-lasso Risk Groups (", data_name, ")", sep = ""),
    surv.median.line = "hv",
    risk.table = TRUE,
    xlab = "Follow-up Time (days)",
    legend = c(0.8, 0.75),
    legend.title = "Risk Group",
    legend.labs = unique(risk_group),
    break.x.by = 100,
    palette = palette_name,
    base_size = base_size
  )
  
  surv_plot <- km_plot_lasso$plot + ggprism::theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot_lasso$table + ggprism::theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  combined_plot <- gridExtra::grid.arrange(surv_plot, risk_table, ncol = 1, heights = c(1.8, 1))
  
  if (save_plot) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    cat("Saving plots to directory:", save_dir, "\n")
    ggsave(file.path(save_dir, paste0("lasso_km_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("lasso_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("lasso_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", save_dir, "\n")
    }
  
  results_lasso <- list(
    KM_test_results = km_results_lasso,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_lasso)
}


#' Evaluate KM for LASSO Model (Wrapper)
#'
#' Wrapper to evaluate KM for LASSO model in PrognosiX object.
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
#' @return Updated PrognosiX object or results.
#' @export
evaluate_km_lasso_model <- function(object,
                                    time_col = "time",
                                    status_col = "status",
                                    nfolds = 10,
                                    save_plot = TRUE,
                                    save_dir = here::here("PrognosiX", "lasso_model"),
                                    plot_width = 7,
                                    plot_height = 7,
                                    base_size = 14,
                                    data_name = "test",
                                    data=NULL) {

  cat("Evaluating Kaplan-Meier for lasso model...\n")

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_best_cv_lasso <- methods::slot(object, "survival.model")[["lasso_model"]][["model"]]
    best_lambda_value <- methods::slot(object, "survival.model")[["lasso_model"]][["lambda"]]

    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$training
    test_data <- data_sets$testing
    cat("Data extracted for training and testing.\n")

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    train_data <- object$training
    test_data <- object$testing
    fit_best_cv_lasso <- object$fit_best_cv_lasso
    best_lambda_value <- object$best_lambda_value

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_lasso_km <- evaluate_lasso_km(
    fit_best_cv_lasso,
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

  cat("Kaplan-Meier evaluation completed. Returning results...\n")

  if (inherits(object, 'PrognosiX')) {
    object@survival.model[["lasso_model"]][["results_km"]] <- results_lasso_km
    return(object)
  }

  return(results_lasso_km)
}


#' Compute HR for LASSO
#'
#' Computes HR and CI for LASSO coefficients.
#'
#' @param lasso_model LASSO model.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param lambda Lambda.
#'
#' @return HR results data frame.
#' @export
compute_hr_lasso <- function(lasso_model,
                             time_col = "time",
                             status_col = "status",
                             lambda =NULL
) {

  lasso_coefs <- coef(lasso_model, s = lambda)
  lasso_coefs <- as.numeric(lasso_coefs)

  # Note: Variable names need to be retrieved or passed correctly.
  # The original code accessed 'train_data' from global scope or similar which is risky.
  # Assuming coefficients order matches features.
  # Better implementation would extract names from model or coefficient object.
  # For now, following structure but beware of missing feature names.
  
  hr <- exp(lasso_coefs)

  assumed_se <- 0.05
  ci_lower <- exp(lasso_coefs - 1.96 * assumed_se)
  ci_upper <- exp(lasso_coefs + 1.96 * assumed_se)

  # Attempt to get names
  feature_names <- rownames(coef(lasso_model, s=lambda))
  
  hr_results <- data.frame(
    Variable = feature_names,
    Coefficient = lasso_coefs,
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


#' Compute HR and CI for LASSO (Wrapper)
#'
#' Wrapper to compute HR/CI for LASSO in PrognosiX.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param lambda Lambda.
#'
#' @return Updated PrognosiX or results.
#' @export
lasso_compute_hr_and_ci <- function(object,
                                    time_col = "time",
                                    status_col = "status",
                                    lambda =NULL
) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_best_cv_lasso <- methods::slot(object, "survival.model")[["lasso_model"]][["model"]]
    lambda<-methods::slot(object, "survival.model")[["lasso_model"]][["lambda"]]
  } else if (is.list(object) && all(c("training", "testing") %in% names(object))) {

    fit_best_cv_lasso <- object$fit_best_cv_lasso

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Evaluating hazard ratios and confidence intervals for lasso model...\n")

  hr_results <- compute_hr_lasso(lasso_model  = fit_best_cv_lasso,
                                 time_col = time_col,
                                 status_col = status_col,
                                 lambda = lambda
  )

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with lasso model results...\n")

    object@survival.model[["lasso_model"]][["hr_results"]] <- hr_results
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")

    return(object)
  }

  cat("Returning Kaplan-Meier results as a list...\n")

  return(hr_results)
}


#' Create Forest Plot for LASSO Model
#'
#' Creates a forest plot for LASSO model HR results.
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
create_forest_plot_lasso <- function(hr_results,
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


#' Forest Plot for LASSO Model (Wrapper)
#'
#' Generates forest plot for LASSO model in PrognosiX.
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
forest_plot_lasso_model <- function(object,
                                    plot_title = "Evaluation of Hazard Ratios and Confidence Intervals for Lasso Model in Survival Analysis",
                                    time_col = "time",
                                    status_col = "status",
                                    var_col = NULL,
                                    palette_name = "AsteroidCity1",
                                    save_plot = TRUE,
                                    save_dir = here::here('PrognosiX', "lasso_model"),
                                    plot_width = 14,
                                    plot_height = 7,
                                    base_size = 14,
                                    use_subgroup_data = FALSE,
                                    hr_limit = c(0.1, 3)
) {

  if (inherits(object, 'PrognosiX')) {
    hr_results <- methods::slot(object, "survival.model")[["lasso_model"]][["hr_results"]]

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

  p <- create_forest_plot_lasso(hr_results,
                                plot_title = plot_title,
                                save_plot = save_plot,
                                save_dir = save_dir,
                                plot_width = plot_width,
                                plot_height = plot_height,
                                base_size = base_size)

  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object with forest plot...\n")
    object@survival.model[["lasso_model"]][["forest_plot"]] <- p
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")
    return(object)
  }

  return(p)
}
