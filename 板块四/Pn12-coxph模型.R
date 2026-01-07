train_coxph_base <- function(train_data,
                             time_col = "time",
                             status_col = "status",
                             nfolds = 5,
                             palette_name = "AsteroidCity1",
                             save_plot = TRUE,
                             save_dir = here("PrognosiX", "coxph_model"),
                             plot_width = 5,
                             plot_height = 5,
                             base_size = 14,
                             seed = 123,
                             kernel_func = NULL) {
  set.seed(seed)

  if (!(time_col %in% colnames(train_data)) || !(status_col %in% colnames(train_data))) {
    stop("Error: time or status column does not exist in the data.")
  }

  if (!is.null(kernel_func)) {
    train_data <- kernel_func(train_data)
  }

  train_data <- train_data[complete.cases(train_data), ]
  if (nrow(train_data) < 2) stop("Insufficient valid data after removing missing values.")

  folds <- caret::createFolds(train_data[[status_col]], k = nfolds, list = FALSE)

  c_statistics <- numeric(nfolds)
  best_c_index <- -Inf
  best_model <- NULL
  skipped_folds <- 0

  for (i in 1:nfolds) {
    training_set <- train_data[folds != i, ]
    test_set <- train_data[folds == i, ]

    if (length(unique(test_set[[status_col]])) < 2) {
      cat("Skipping fold", i, "due to insufficient event diversity in the test set\n")
      skipped_folds <- skipped_folds + 1
      next
    }

    cox_model <- tryCatch({
      coxph(Surv(training_set[[time_col]], training_set[[status_col]]) ~ ., data = training_set)
    }, error = function(e) {
      cat("Skipping fold", i, "due to model fitting error:", e$message, "\n")
      skipped_folds <- skipped_folds + 1
      return(NULL)
    })

    if (is.null(cox_model)) next

    predicted_risk <- tryCatch({
      predict(cox_model, newdata = test_set, type = "risk")
    }, error = function(e) {
      cat("Skipping fold", i, "due to prediction error:", e$message, "\n")
      skipped_folds <- skipped_folds + 1
      return(NULL)
    })

    if (is.null(predicted_risk)) next

    predicted_risk <- pmax(pmin(predicted_risk, 1e6), -1e6)

    c_index_result <- tryCatch({
      survcomp::concordance.index(predicted_risk, test_set[[time_col]], test_set[[status_col]])
    }, error = function(e) {
      cat("Skipping fold", i, "due to C-index calculation error:", e$message, "\n")
      skipped_folds <- skipped_folds + 1
      return(NULL)
    })

    if (is.null(c_index_result) || is.na(c_index_result$c.index) || c_index_result$c.index == 0) {
      cat("Skipping fold", i, "due to invalid C-index\n")
      skipped_folds <- skipped_folds + 1
      next
    }

    c_statistics[i] <- c_index_result$c.index
    cat("C-index for fold", i, ":", c_statistics[i], "\n")

    if (c_statistics[i] > best_c_index) {
      best_c_index <- c_statistics[i]
      best_model <- cox_model
    }
  }

  if (skipped_folds > nfolds / 2) warning("More than 50% of folds were skipped.")
  if (skipped_folds == nfolds) stop("All folds were skipped due to missing or invalid data.")

  palette <-  wes_palette(palette_name, type = "discrete")
  cv_data <- data.frame(fold = 1:nfolds, c_index = c_statistics)

  p <- ggplot(cv_data, aes(x = fold, y = c_index)) +
    geom_line(color = palette[1]) +
    geom_ribbon(aes(ymin = c_index - 0.01, ymax = c_index + 0.01), fill = palette[2], alpha = 0.5) +
    scale_x_continuous() +
    labs(x = "Fold", y = "Concordance Index (C-index)", title = "Cross-Validation C-index for CoxPH model") +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

  print(p)

  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, "cv_coxph_plot.pdf"), plot = p, width = plot_width, height = plot_height, device = "pdf")
    
    cat("Plot saved to:", file.path(save_dir, "cv_coxph_plot.pdf"), "\n")
    }

  return(list(model = best_model, mean_cv_error = mean(c_statistics, na.rm = TRUE),
              cv_errors = c_statistics, skipped_folds = skipped_folds, cv_plot = p))
}


train_coxph_model <- function(object,
                              time_col ,
                              status_col ,
                              nfolds = 5,
                              palette_name = "AsteroidCity1",
                              save_plot = TRUE,
                              save_dir = here("PrognosiX", "coxph_model"),
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 14,
                              seed = 123,
                              kernel_func = NULL) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    train_data <- object$training
    test_data <- object$testing

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Training CoxPH model...\n")

  coxph_model_info <- train_coxph_base(
    train_data = train_data,
    time_col = time_col,
    status_col = status_col,
    nfolds = nfolds,
    palette_name = palette_name,
    save_plot = save_plot,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    seed = seed,
    kernel_func = kernel_func
  )

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with CoxPH model results...\n")

    object@survival.model[["coxph_model"]] <- coxph_model_info

    return(object)
  }

  cat("Returning model information as a list...\n")

  return(coxph_model_info)
}


evaluate_coxph_roc <- function(
    fit_coxph_model,
    train_data = NULL,
    test_data = NULL,
    validation_data = NULL,
    palette_name = "AsteroidCity1",
    time_col = "time",
    status_col = "status",
    save_plot = TRUE,
    save_dir = here("PrognosiX", "coxph_model"),
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
    train_data_clean <- train_data[complete.cases(train_data[, c(time_col, status_col)]), ]
    risk_score_train <- predict(fit_coxph_model, newdata = train_data_clean, type = "risk")
    
    risk_score_train[is.infinite(risk_score_train)] <- NA 
    risk_score_train[is.na(risk_score_train)] <- median(risk_score_train, na.rm = TRUE)
    risk_score_train <- as.vector(risk_score_train)
    train_data_clean[[status_col]] <- as.numeric(as.factor(train_data_clean[[status_col]]))
    roc_train_curve <- roc(train_data_clean[[status_col]], risk_score_train)
    roc_train_auc <- auc(roc_train_curve)
    
    roc_data_train <- data.frame(
      specificity = roc_train_curve$specificities,
      sensitivity = roc_train_curve$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = survcomp::concordance.index(
        x = risk_score_train,
        surv.time = train_data_clean[[time_col]],
        surv.event = train_data_clean[[status_col]],
        method = "noether"
      )$c.index,
      ROC_AUC = roc_train_auc
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    test_data_clean <- test_data[complete.cases(test_data[, c(time_col, status_col)]), ]
    risk_score_test <- predict(fit_coxph_model, newdata = test_data_clean, type = "risk")
   
    risk_score_test[is.infinite(risk_score_test)] <- NA 
    risk_score_test[is.na(risk_score_test)] <- median(risk_score_test, na.rm = TRUE)
     risk_score_test <- as.vector(risk_score_test)
    test_data_clean[[status_col]] <- as.numeric(as.factor(test_data_clean[[status_col]]))
    roc_test_curve <- roc(test_data_clean[[status_col]], risk_score_test)
    roc_test_auc <- auc(roc_test_curve)
    
    roc_data_test <- data.frame(
      specificity = roc_test_curve$specificities,
      sensitivity = roc_test_curve$sensitivities,
      Set = "Testing Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_test)
    
    Test <- data.frame(
      C_index = survcomp::concordance.index(
        x = risk_score_test,
        surv.time = test_data_clean[[time_col]],
        surv.event = test_data_clean[[status_col]],
        method = "noether"
      )$c.index,
      ROC_AUC = roc_test_auc
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
  }
  
  if (!is.null(validation_data)) {
    validation_data_clean <- validation_data[complete.cases(validation_data[, c(time_col, status_col)]), ]
    risk_score_validation <- predict(fit_coxph_model, newdata = validation_data_clean, type = "risk")
    
    risk_score_validation[is.infinite(risk_score_validation)] <- NA 
    risk_score_validation[is.na(risk_score_validation)] <- median(risk_score_validation, na.rm = TRUE)
    risk_score_validation <- as.vector(risk_score_validation)
    
    validation_data_clean[[status_col]] <- as.numeric(as.factor(validation_data_clean[[status_col]]))
    
    roc_validation_curve <- roc(validation_data_clean[[status_col]], risk_score_validation)
    roc_validation_auc <- auc(roc_validation_curve)
    
    roc_data_validation <- data.frame(
      specificity = roc_validation_curve$specificities,
      sensitivity = roc_validation_curve$sensitivities,
      Set = "Validation Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_validation)
    
    Validation <- data.frame(
      C_index = survcomp::concordance.index(
        x = risk_score_validation,
        surv.time = validation_data_clean[[time_col]],
        surv.event = validation_data_clean[[status_col]],
        method = "noether"
      )$c.index,
      ROC_AUC = roc_validation_auc
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Validation", Validation))
  }
  
  roc_data_all <- roc_data_all %>%
    filter(!(sensitivity == 0 & specificity == 1)) %>%
    filter(sensitivity >= 0 & sensitivity <= 1, specificity >= 0 & specificity <= 1)
  
  roc_plot <- ggplot(roc_data_all, aes(x = 1 - specificity, y = sensitivity, color = Set)) +
    geom_path(aes(color = Set), size = 1.5) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = wes_palette(palette_name)) +
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
    ggsave(filename = file.path(save_dir, "ROC_curves.pdf"), plot = roc_plot, width = plot_width, height = plot_height,device = "pdf")
    
    cat("Plot saved to:", file.path(save_dir, "ROC_curves.pdf"), "\n")
  }
  
  print(results_roc_df)
  
  results_roc <- list(
    results_roc = results_roc_df,
    plot_roc = roc_plot
  )
  
  return(results_roc)
}



evaluate_roc_coxph_model <- function(object,
                                     time_col = "time",
                                     status_col = "status",
                                     save_plot = TRUE,
                                     save_dir = here("PrognosiX", "coxph_model"),
                                     plot_width = 7,
                                     plot_height = 7,
                                     base_size = 14,
                                     palette_name = "AsteroidCity1") {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_coxph_model <- slot(object, "survival.model")[["coxph_model"]][["model"]]

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    train_data <- object$training
    test_data <- object$testing

    fit_coxph_model <- object$fit_coxph_model
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_roc <- evaluate_coxph_roc(fit_coxph_model,
                                    train_data, test_data,
                                    palette_name = palette_name,
                                    time_col = time_col, status_col = status_col,
                                    save_plot = save_plot, save_dir = save_dir,
                                    plot_width = plot_width, plot_height = plot_height,
                                    base_size = base_size)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is updated with new ROC results.\n")
    object@survival.model[["coxph_model"]][["results_roc"]] <- results_roc
    return(object)
  }

  return(results_roc)
}


evaluate_coxph_km <- function(
    fit_best_cv_coxph, data, data_name = "test",
    time_col = "time", status_col = "status",
    save_plot = TRUE, save_dir = here("PrognosiX", "coxph_model"),
    palette_name = "Dark2", plot_width = 10, plot_height = 8,
    base_size = 12, seed = 1234, fill_color = "lightblue"
) {
  set.seed(seed)
  
  if (!data_name %in% c("train", "test", "val")) {
    stop("'data_name' must be one of 'train', 'test', or 'val'.")
  }
  
  data_x <- as.data.frame(data)
  
  risk_score <- predict(fit_best_cv_coxph, newdata = data_x, type = "risk")
  cat("Risk scores calculated for", data_name, "data.\n")
  
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  data <- cbind(data, risk_score, risk_group)
  
  colnames(data)[ncol(data) - 1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  fit_km_coxph <- survfit(Surv(data[[time_col]], data[[status_col]]) ~ risk_group, data = data)
  cat("Kaplan-Meier model fitted for", data_name, "data.\n")
  
  km_pval_coxph <- pchisq(survdiff(Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)$chisq, 1, lower.tail = FALSE)
  
  coxph_fit <- coxph(Surv(data[[time_col]], data[[status_col]]) ~ risk_group, data = data)
  km_hr_coxph <- exp(coef(coxph_fit))
  km_ci_coxph <- confint(coxph_fit)
  
  km_results_coxph <- data.frame(
    KM_HR = km_hr_coxph,
    KM_CI_lower = km_ci_coxph[1],
    KM_CI_upper = km_ci_coxph[2],
    KM_p_value = km_pval_coxph
  )
  cat("Kaplan-Meier results stored in a dataframe for", data_name, "data.\n")
  
  km_plot_coxph <- ggsurvplot(
    fit_km_coxph, data = data, conf.int = TRUE,
    conf.int.fill = fill_color, conf.int.alpha = 0.5,
    pval = TRUE, pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve for Coxph Risk Groups (", data_name, ")", sep = ""),
    surv.median.line = "hv", risk.table = TRUE,
    xlab = "Follow-up Time (days)", legend = c(0.8, 0.75),
    legend.title = "Risk Group", legend.labs = unique(data$risk_group),
    break.x.by = 100, palette = palette_name, base_size = base_size
  )
  print(km_plot_coxph)
  
  cat("Kaplan-Meier plot generated for", data_name, "data.\n")
  
  surv_plot <- km_plot_coxph$plot + theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot_coxph$table + theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  combined_plot <- grid.arrange(surv_plot, risk_table, ncol = 1, heights = c(1.8, 1))
  
  cat("Combined plot generated for", data_name, "data.\n")
  
  if (save_plot) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat(paste("Created directory:", save_dir, "\n"))
    }
    
    ggsave(file.path(save_dir, paste0("coxph_time_distribution_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("coxph_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("coxph_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    
    cat("Plot saved to:", save_dir ,"\n")
  }
  
  results_coxph <- list(
    KM_test_results = km_results_coxph,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_coxph)
}

evaluate_km_coxph_model <- function(object,
                                    time_col = "time",
                                    status_col = "status",
                                    save_plot = TRUE,
                                    save_dir = here("PrognosiX", "coxph_model"),
                                    palette_name = "Dark2",
                                    plot_width = 10,
                                    plot_height = 8,
                                    base_size = 12,
                                    seed = 1234,
                                    fill_color = "lightblue",
                                    data_name = "test",
                                    data=NULL) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_best_cv_coxph <- slot(object, "survival.model")[["coxph_model"]][["model"]]

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing
  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    train_data <- object$training
    test_data <- object$testing

    fit_best_cv_coxph <- object$fit_best_cv_coxph

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_coxph_km <- evaluate_coxph_km(
    fit_best_cv_coxph = fit_best_cv_coxph,
    data = test_data,
    data_name = data_name,
    time_col = time_col,
    status_col = status_col,
    save_plot = save_plot,
    save_dir = save_dir,
    palette_name = palette_name,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    seed = seed,
    fill_color = fill_color)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with coxph model Kaplan-Meier results...\n")

    object@survival.model[["coxph_model"]][["results_km"]] <- results_coxph_km
    return(object)
  }

  return(results_coxph_km)
}





compute_hr_coxph <- function(cox_model,
                             time_col = "time",
                             status_col = "status"
) {

  coxph_coefs <- coef(cox_model)
  coxph_se <- sqrt(diag(cox_model$var))

  hr <- exp(coxph_coefs)
  ci_lower <- exp(coxph_coefs - 1.96 * coxph_se)
  ci_upper <- exp(coxph_coefs + 1.96 * coxph_se)

  p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]

  hr_results <- data.frame(
    Variable = names(coxph_coefs),
    HR = hr,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    P_value = p_values
  )

  hr_results$HR_95CI <- paste0(
    round(hr_results$HR, 2),
    " (",
    round(hr_results$CI_lower, 2),
    "-",
    round(hr_results$CI_upper, 2),
    ")"
  )
  hr_results$se <- (log(hr_results$CI_upper) - log(hr_results$CI_lower)) / (2 * qnorm(0.975))

  hr_results$P_value <- format.pval(hr_results$P_value, digits = 3)

  print(hr_results)

  return(hr_results)
}


coxph_compute_hr_and_ci <- function(object,
                                    time_col = "time",
                                    status_col = "status"
) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_best_cv_coxph <- slot(object, "survival.model")[["coxph_model"]][["model"]]

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {

    fit_best_cv_coxph <- object$fit_best_cv_coxph

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Evaluating hazard ratios and confidence intervals for Cox model...\n")

  hr_results <- compute_hr_coxph(cox_model = fit_best_cv_coxph,
                                 time_col = time_col,
                                 status_col = status_col
  )

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with Cox model results...\n")

    object@survival.model[["coxph_model"]][["hr_results"]] <- hr_results

    return(object)
  }

  cat("Returning Kaplan-Meier results as a list...\n")

  return(hr_results)
}



create_forest_model_plot <- function(hr_results,
                                     plot_title = "Evaluation of Hazard Ratios and Confidence Intervals",
                                     save_plot = FALSE,
                                     save_dir = here('Prognosi'),
                                     plot_width = 11,
                                     plot_height = 5,
                                     palette_name = "AsteroidCity1",
                                     base_size = 14,
                                     hr_limit = c(0, 3),
                                     ci_range_limit = 1000) {

  hr_results <- hr_results %>%
    filter(is.finite(HR) & is.finite(CI_lower) & is.finite(CI_upper)) %>%
    filter(HR >= hr_limit[1] & HR <= hr_limit[2]) %>%
    filter((CI_upper - CI_lower) <= ci_range_limit)

  cat(paste(nrow(hr_results), "rows remain after filtering. Rows with extreme or invalid values were excluded.\n"))

  tab_header <- c("Clinical factors", "HR (95% CI)")
  tab_header_bold <- TRUE

  tmp_df <- hr_results %>%
    dplyr::select(Variable, HR_95CI) %>%
    rbind(tab_header, .) %>%
    cbind("clinical" = .[, "Variable"], .)

  tmp_df$clinical <- factor(tmp_df$clinical, levels = rev(tmp_df$clinical))

  tmp_df <- tmp_df %>%
    pivot_longer(cols = 2:ncol(.), names_to = "x", values_to = "label")

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

  p_tab <- ggplot(data = tmp_df, aes(x = x, y = clinical)) +
    geom_tile(color = "white", fill = "white") +
    geom_richtext(aes(label = label_bold), label.color = NA) +
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

  p_errorbar <- ggplot(data = tmp_error) +
    geom_point(aes(x = HR, y = Variable, color = HR > 1), shape = "diamond", size = 4) +
    scale_color_manual(values = wes_palette(palette_name)) +
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
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.5)))

  print(p)

  if (save_plot) {
    pdf_path <- file.path(save_dir, "forest_plot.pdf")
    ggsave(filename = pdf_path, plot = p, width = plot_width, height = plot_height)
    cat("Plot saved at: ", pdf_path, "\n")
  }
}



forest_plot_coxph_model <- function(object,
                                    plot_title = "Evaluation of Hazard Ratios and Confidence Intervals for coxph Model in Survival Analysis",
                                    time_col = "time",
                                    status_col = "status",
                                    var_col = NULL,
                                    palette_name = "AsteroidCity1",
                                    save_plot = TRUE,
                                    save_dir = here('PrognosiX', "coxph_model"),
                                    plot_width = 14,
                                    plot_height = 7,
                                    base_size = 14,
                                    use_subgroup_data = FALSE,
                                    hr_limit = c(0.1, 3)
) {

  if (inherits(object, 'PrognosiX')) {
    hr_results <- object@survival.model[["coxph_model"]][["hr_results"]]

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

  p <- create_forest_model_plot(hr_results,
                                plot_title = plot_title,
                                save_plot = save_plot,
                                save_dir = save_dir,
                                plot_width = plot_width,
                                plot_height = plot_height,
                                base_size = base_size)

  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object with forest plot...\n")
    object@survival.model[["coxph_model"]][["forest_plot"]] <- p
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")
    return(object)
  }

  return(p)
}


