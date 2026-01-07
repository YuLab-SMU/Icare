train_superpc_base <- function(train_data,
                               time_col = "time",
                               status_col = "status",
                               nfolds = 10,
                               s0_perc = 0.5,
                               n_threshold = 20,
                               n_components = 3,
                               min_features = 3,
                               max_features = NULL,
                               compute_fullcv = TRUE,
                               compute_preval = TRUE,
                               base_size = 14,
                               seed = 1234) {

  set.seed(seed)
  x_train <- t(as.matrix(train_data[, !(names(train_data) %in% c(time_col, status_col))]))
  y_train <- train_data[[time_col]]
  x_train <- apply(x_train, 2, as.numeric)

  train <- list(
    x = x_train,
    y = y_train,
    censoring.status = train_data[[status_col]],
    featurenames = colnames(train_data)[!(colnames(train_data) %in% c(time_col, status_col))]
  )


  if (is.null(max_features)) {
    max_features <- nrow(train$x)
  }

  cat("Training SuperPC model...\n")

  fit_train_superpc <- superpc.train(data = train, type = 'survival', s0.perc = s0_perc)
  summary(fit_train_superpc)


  cat("Running cross-validation...\n")
  cv_superpc_fit <- superpc.cv(fit_train_superpc,
                               train,
                               n.threshold = n_threshold,
                               n.fold = nfolds,
                               n.components = n_components,
                               min.features = min_features,
                               max.features = max_features,
                               compute.fullcv = compute_fullcv,
                               compute.preval = compute_preval)

  cat("Extracting cross-validation scores and thresholds...\n")

  cv_plot <- superpc.plotcv(cv_superpc_fit)

  print(cv_plot)

  best_threshold <- cv_superpc_fit$thresholds[which.max(cv_superpc_fit[["scor"]][1,])]
  cat("Best Threshold: ", best_threshold, "\n")
  cat("Returning model and cross-validation results...\n")
  return(list(model = fit_train_superpc,
              cv_fit = cv_superpc_fit,
              threshold = best_threshold,
              cv_plot = cv_plot))
}

train_superpc_model <- function(object,
                                time_col = "time",
                                status_col = "status",
                                nfolds = 10,
                                s0_perc = 0.5,
                                n_threshold = 20,
                                n_components = 3,
                                min_features = 3,
                                max_features = NULL,
                                compute_fullcv = TRUE,
                                compute_preval = TRUE,
                                base_size = 14,
                                seed = 1234) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")

    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("training", "testing") %in% names(object))) {
    train_data <- object$training
    test_data <- object$testing

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  superpc_model_info <- train_superpc_base(train_data = train_data,
                                           time_col = time_col,
                                           status_col = status_col,
                                           nfolds = nfolds,
                                           s0_perc = s0_perc,
                                           n_threshold = n_threshold,
                                           n_components = n_components,
                                           min_features = min_features,
                                           max_features = max_features,
                                           compute_fullcv = compute_fullcv,
                                           compute_preval = compute_preval,
                                           base_size = base_size,
                                           seed = seed)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with SuperPC model results...\n")

    object@survival.model[["superpc_model"]] <- superpc_model_info

    return(object)
  }

  cat("Returning model information as a list...\n")

  return(superpc_model_info)
}


evaluate_superpc_roc <- function(
    fit_best_cv_superpc,
    train_data = NULL,
    test_data = NULL,
    validation_data = NULL,
    threshold,
    cv_superpc_fit,
    palette_name = "AsteroidCity1",
    time_col = "time",
    status_col = "status",
    save_plot = TRUE,
    save_dir = here("PrognosiX", "superpc_model"),
    plot_width = 7,
    plot_height = 7,
    base_size = 14,
    n_components = 3
) {
  if (is.null(train_data) && is.null(test_data) && is.null(validation_data)) {
    stop("At least one of 'train_data', 'test_data', or 'validation_data' must be provided.")
  }
  
  roc_data_all <- data.frame()
  results_roc_df <- data.frame()
  
  if (!is.null(train_data)) {
    x_train_superpc <- t(as.matrix(train_data[, !(names(train_data) %in% c(time_col, status_col))]))
    y_train_superpc <- train_data[[time_col]]
    x_train_superpc <- apply(x_train_superpc, 2, as.numeric)
    
    train <- list(
      x = x_train_superpc,
      y = y_train_superpc,
      censoring.status = train_data[[status_col]],
      featurenames = colnames(train_data)[!(colnames(train_data) %in% c(time_col, status_col))]
    )
    
    pred_train_superpc <- superpc.predict(fit_best_cv_superpc,
                                          data = train,
                                          newdata = train,
                                          threshold = threshold,
                                          n.components = n_components)
    
    pred_train_superpc <- pred_train_superpc[["v.pred.1df"]]
    pred_train_superpc <- as.numeric(pred_train_superpc)
    pred_train_superpc[is.infinite(pred_train_superpc)] <- NA 
    pred_train_superpc[is.na(pred_train_superpc)] <- median(pred_train_superpc, na.rm = TRUE)
    
    
    roc_train_curve <- roc(train_data[[status_col]], pred_train_superpc)
    roc_train_auc <- auc(roc_train_curve)
    cat("ROC for training set calculated. AUC: ", round(roc_train_auc, 2), "\n")
    
    roc_data_train <- data.frame(
      specificity = roc_train_curve$specificities,
      sensitivity = roc_train_curve$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = concordance.index(
        x = pred_train_superpc,
        surv.time = train_data[[time_col]],
        surv.event = train_data[[status_col]],
        method = "noether"
      )$c.index,
      ROC_AUC = roc_train_auc
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    x_test_superpc <- t(as.matrix(test_data[, !(names(test_data) %in% c(time_col, status_col))]))
    y_test_superpc <- test_data[[time_col]]
    x_test_superpc <- apply(x_test_superpc, 2, as.numeric)
    
    test <- list(
      x = x_test_superpc,
      y = y_test_superpc,
      censoring.status = test_data[[status_col]],
      featurenames = colnames(test_data)[!(colnames(test_data) %in% c(time_col, status_col))]
    )
    
    pred_test_superpc <- superpc.predict(fit_best_cv_superpc,
                                         data = train,
                                         newdata = test,
                                         threshold = threshold,
                                         n.components = n_components)
    
    pred_test_superpc <- pred_test_superpc[["v.pred.1df"]]
    pred_test_superpc <- as.numeric(pred_test_superpc)
    pred_test_superpc[is.infinite(pred_test_superpc)] <- NA 
    pred_test_superpc[is.na(pred_test_superpc)] <- median(pred_test_superpc, na.rm = TRUE)
    
    if (length(test_data[[status_col]]) == length(pred_test_superpc)) {
      cat("Prediction vector length matches the test set.\n")
      
      roc_test_curve <- roc(test_data[[status_col]], pred_test_superpc)
      roc_test_auc <- auc(roc_test_curve)
      cat("ROC for test set calculated. AUC: ", round(roc_test_auc, 2), "\n")
      
      roc_data_test <- data.frame(
        specificity = roc_test_curve$specificities,
        sensitivity = roc_test_curve$sensitivities,
        Set = "Testing Set"
      )
      roc_data_all <- rbind(roc_data_all, roc_data_test)
      
      Test <- data.frame(
        C_index = concordance.index(
          x = pred_test_superpc,
          surv.time = test_data[[time_col]],
          surv.event = test_data[[status_col]],
          method = "noether"
        )$c.index,
        ROC_AUC = roc_test_auc
      )
      results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
    } else {
      cat("Error: Length mismatch between response and prediction vectors.\n")
    }
  }
  
  # Process validation data
  if (!is.null(validation_data)) {
    x_validation_superpc <- t(as.matrix(validation_data[, !(names(validation_data) %in% c(time_col, status_col))]))
    y_validation_superpc <- validation_data[[time_col]]
    x_validation_superpc <- apply(x_validation_superpc, 2, as.numeric)
    
    validation <- list(
      x = x_validation_superpc,
      y = y_validation_superpc,
      censoring.status = validation_data[[status_col]],
      featurenames = colnames(validation_data)[!(colnames(validation_data) %in% c(time_col, status_col))]
    )
    
    pred_validation_superpc <- superpc.predict(fit_best_cv_superpc,
                                               data = train,
                                               newdata = validation,
                                               threshold = threshold,
                                               n.components = n_components)
    
    pred_validation_superpc <- pred_validation_superpc[["v.pred.1df"]]
    pred_validation_superpc <- as.numeric(pred_validation_superpc)
    pred_validation_superpc[is.infinite(pred_validation_superpc)] <- NA 
    pred_validation_superpc[is.na(pred_validation_superpc)] <- median(pred_validation_superpc, na.rm = TRUE)
    
    if (length(validation_data[[status_col]]) == length(pred_validation_superpc)) {
      cat("Prediction vector length matches the validation set.\n")
      
      roc_validation_curve <- roc(validation_data[[status_col]], pred_validation_superpc)
      roc_validation_auc <- auc(roc_validation_curve)
      cat("ROC for validation set calculated. AUC: ", round(roc_validation_auc, 2), "\n")
      
      roc_data_validation <- data.frame(
        specificity = roc_validation_curve$specificities,
        sensitivity = roc_validation_curve$sensitivities,
        Set = "Validation Set"
      )
      roc_data_all <- rbind(roc_data_all, roc_data_validation)
      
      Validation <- data.frame(
        C_index = concordance.index(
          x = pred_validation_superpc,
          surv.time = validation_data[[time_col]],
          surv.event = validation_data[[status_col]],
          method = "noether"
        )$c.index,
        ROC_AUC = roc_validation_auc
      )
      results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Validation", Validation))
    } else {
      cat("Error: Length mismatch between response and prediction vectors.\n")
    }
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
    cat("Plot saved to:", file.path(save_dir, "ROC_curves.pdf"), "\n")  }
  
  print(results_roc_df)
  
  results_roc <- list(
    results_roc = results_roc_df,
    plot_roc = roc_plot
  )
  
  return(results_roc)
}





evaluate_roc_superpc_model <- function(object,
                                       time_col = "time",
                                       status_col = "status",
                                       threshold = 0.5,
                                       n_components = 3,
                                       save_plot = TRUE,
                                       save_dir = here("PrognosiX", "superpc_model"),
                                       plot_width = 7,
                                       plot_height = 7,
                                       base_size = 14,
                                       palette_name = "AsteroidCity1") {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_superpc_model <- slot(object, "survival.model")[["superpc_model"]][["model"]]
    threshold <- slot(object, "survival.model")[["superpc_model"]][["threshold"]]



    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$training
    test_data <- data_sets$testing

    cat("Training data and test data extracted successfully.\n")

  } else if (is.list(object) && all(c("training", "testing") %in% names(object))) {
    cat("Input is a list with 'training' and 'testing' elements.\n")

    train_data <- object$training
    test_data <- object$testing

    fit_superpc_model <- object$fit_superpc_model
    threshold <- object$threshold

    cat("model and data extracted from the list successfully.\n")

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Calling evaluate_superpc_roc function to calculate ROC and AUC...\n")

  results_roc <- evaluate_superpc_roc(fit_best_cv_superpc = fit_superpc_model,
                                      train_data = train_data,
                                      test_data = test_data,
                                      threshold = threshold,
                                      cv_superpc_fit = threshold,
                                      palette_name = palette_name,
                                      time_col = time_col,
                                      status_col = status_col,
                                      save_plot = save_plot,
                                      save_dir = save_dir,
                                      plot_width = plot_width,
                                      plot_height = plot_height,
                                      base_size = base_size)


  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is updated with new ROC results.\n")
    object@survival.model[["superpc_model"]][["results_roc"]] <- results_roc
    return(object)
  }

  cat("Returning ROC results...\n")
  return(results_roc)
}

evaluate_superpc_km <- function(
    fit_best_cv_superpc, data, data_name = "test",
    time_col = "time", status_col = "status",
    threshold = 0.5, save_plot = TRUE,
    save_dir = here("PrognosiX", "superpc_model"),
    palette_name = "Dark2", plot_width = 10, plot_height = 8,
    base_size = 12, seed = 1234, n.components = 1,
    fill_color = "lightblue"
) {
  set.seed(seed)
  
  if (!data_name %in% c("train", "test", "val")) {
    stop("'data_name' must be one of 'train', 'test', or 'val'.")
  }
  
  x_data_superpc <- t(as.matrix(data[, !(names(data) %in% c(time_col, status_col))]))
  y_data_superpc <- data[[time_col]]
  x_data_superpc <- apply(x_data_superpc, 2, as.numeric)
  
  data_list <- list(
    x = x_data_superpc,
    y = y_data_superpc,
    censoring.status = data[[status_col]],
    featurenames = colnames(data)[!(colnames(data) %in% c(time_col, status_col))]
  )
  
  risk_score <- superpc.predict(fit_best_cv_superpc,
                                data = data_list,
                                newdata = data_list,
                                threshold = threshold,
                                n.components = n.components)
  
  risk_score <- as.numeric(risk_score$v.pred)
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  
  cat("Risk groups assigned based on predicted scores for", data_name, "data.\n")
  
  data <- cbind(data, risk_score, risk_group)
  colnames(data)[ncol(data) - 1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  fit_km_superpc <- survfit(Surv(data[[time_col]], data[[status_col]]) ~ risk_group, data = data)
  sdf_superpc <- survdiff(Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  
  km_pval_superpc <- pchisq(sdf_superpc$chisq, length(sdf_superpc$n) - 1, lower.tail = FALSE)
  km_hr_superpc <- (sdf_superpc$obs[2] / sdf_superpc$exp[2]) / (sdf_superpc$obs[1] / sdf_superpc$exp[1])
  
  km_upper95_superpc <- exp(log(km_hr_superpc) + qnorm(0.975) * sqrt(1 / sdf_superpc$exp[2] + 1 / sdf_superpc$exp[1]))
  km_lower95_superpc <- exp(log(km_hr_superpc) - qnorm(0.975) * sqrt(1 / sdf_superpc$exp[2] + 1 / sdf_superpc$exp[1]))
  
  km_results_superpc <- data.frame(
    KM_HR = km_hr_superpc,
    KM_CI_lower = km_lower95_superpc,
    KM_CI_upper = km_upper95_superpc,
    KM_p_value = km_pval_superpc
  )
  cat("Kaplan-Meier results calculated for", data_name, "data.\n")
  
  km_plot_superpc <- ggsurvplot(
    fit_km_superpc, data = data, conf.int = TRUE,
    conf.int.fill = fill_color, conf.int.alpha = 0.5,
    pval = TRUE, pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve by Risk Group (SuperPC) (", data_name, ")", sep = ""),
    surv.median.line = "hv", risk.table = TRUE,
    xlab = "Follow up time (days)", legend = c(0.8, 0.75),
    legend.title = "Risk Group (SuperPC)", legend.labs = unique(risk_group),
    break.x.by = 100, palette = palette_name, base_size = base_size
  )
  print(km_plot_superpc)
  
  cat("Kaplan-Meier plot generated for", data_name, "data.\n")
  
  surv_plot <- km_plot_superpc$plot + theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot_superpc$table + theme_prism(base_size = base_size) +
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
    
    ggsave(file.path(save_dir, paste0("superpc_time_distribution_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("superpc_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("superpc_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", save_dir, "\n")
  }
  
  results_superpc <- list(
    KM_test_results = km_results_superpc,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_superpc)
}


evaluate_km_superpc_model <- function(object,
                                      time_col = "time",
                                      status_col = "status",
                                      threshold = 0.5,
                                      save_plot = TRUE,
                                      save_dir = here("PrognosiX", "superpc_model"),
                                      palette_name = "Dark2",
                                      plot_width = 10,
                                      plot_height = 8,
                                      base_size = 12,
                                      seed = 1234,
                                      n.components = 1,
                                      data_name = "test",
                                      data=NULL) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_best_cv_superpc <- slot(object, "survival.model")[["superpc_model"]][["model"]]
    threshold <- slot(object, "survival.model")[["superpc_model"]][["threshold"]]
    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    train_data <- object$train
    test_data <- object$test
    fit_best_cv_superpc <- object$fit_best_cv_superpc
    threshold <- object$threshold

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_superpc_km <- evaluate_superpc_km(
    fit_best_cv_superpc,
    data = test_data,
    data_name = data_name,
    time_col = time_col,
    status_col = status_col,
    save_plot = save_plot,  # Corrected variable name here
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    n.components = n.components,
    threshold = threshold
  )

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with SuperPC model Kaplan-Meier results...\n")

    object@survival.model[["superpc_model"]][["results_km"]] <- results_superpc_km

    return(object)
  }

  return(results_superpc_km)
}


compute_hr_surperpc <- function(surperpc_model, train_data,
                                time_col = "time",
                                status_col = "status",
                                assumed_se = 0.05) {

  # Extract feature scores (coefficients)
  surperpc_coefs <- surperpc_model$feature.scores
  surperpc_coefs <- as.numeric(surperpc_coefs)

  # Feature names (excluding time and status columns)
  feature_names <- colnames(train_data[, !(names(train_data) %in% c(time_col, status_col))])

  # Compute Hazard Ratios (HR)
  hr <- exp(surperpc_coefs)

  # Assuming a fixed standard error (or use bootstrapping if you prefer)
  ci_lower <- exp(surperpc_coefs - 1.96 * assumed_se)
  ci_upper <- exp(surperpc_coefs + 1.96 * assumed_se)

  # Prepare the results data frame
  hr_results <- data.frame(
    Variable = feature_names,
    Coefficient = surperpc_coefs,
    HR = hr,
    CI_lower = ci_lower,
    CI_upper = ci_upper
  )

  # Create a column for HR with 95% CI
  hr_results$HR_95CI <- paste0(
    round(hr_results$HR, 2), " (",
    round(hr_results$CI_lower, 2), "-",
    round(hr_results$CI_upper, 2), ")"
  )

  # Print the results
  print(hr_results)

  return(hr_results)
}

surperpc_compute_hr_and_ci <- function(object,
                                       time_col = "time",
                                       status_col = "status",
                                       assumed_se = 0.05) {
  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    fit_best_cv_surperpc_model <- object@survival.model[["superpc_model"]][["model"]]
    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing
    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
  } else if (is.list(object) && all(c("fit_best_cv_surperpc_model") %in% names(object))) {
    fit_best_cv_surperpc_model <- object$fit_best_cv_surperpc_model
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'fit_best_cv_surperpc_model'.")
  }

  cat("Evaluating hazard ratios and confidence intervals for SuperPC model...\n")
  hr_results <- compute_hr_surperpc(
    surperpc_model = fit_best_cv_surperpc_model,
    train_data = train_data,
    status_col = status_col,
    time_col = time_col,
    assumed_se=assumed_se

  )

  if (inherits(object, 'PrognosiX')) {
    object@survival.model[["superpc_model"]][["hr_results"]] <- hr_results
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")

    return(object)
  }

  cat("Returning hazard ratio results as a list...\n")
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



forest_plot_superpc_model <- function(object,
                                      plot_title = "Evaluation of Hazard Ratios and Confidence Intervals for superpc Model in Survival Analysis",
                                      time_col = "time",
                                      status_col = "status",
                                      var_col = NULL,
                                      palette_name = "AsteroidCity1",
                                      save_plot = TRUE,
                                      save_dir = here('PrognosiX', "superpc_model"),
                                      plot_width = 14,
                                      plot_height = 7,
                                      base_size = 14,
                                      use_subgroup_data = FALSE,
                                      hr_limit = c(0.1, 3)
) {

  if (inherits(object, 'PrognosiX')) {
    hr_results <- object@survival.model[["superpc_model"]][["hr_results"]]

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
    object@survival.model[["superpc_model"]][["forest_plot"]] <- p
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")
    return(object)
  }

  return(p)
}


