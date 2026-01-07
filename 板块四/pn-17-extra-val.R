ExtractPrognosiXValidata <- function(
    data = NULL,
    object_stats = NULL,
    object_prognosix = NULL,
    group_col = NULL,
    time_col = "time",  
    status_col = "status",  
    ...
) {
  if (is.null(object_prognosix)) {
    stop("The 'object_prognosix' parameter must be provided.")
  }
  
  if (!is.null(data) && !is.null(object_stats)) {
    stop("Only one of 'data' and 'object_stats' should be provided.")
  }
    if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("The 'data' parameter must be a data frame.")
    }
    
    if (any(is.na(data))) {
      stop("The data frame contains missing values. Please handle them before proceeding.")
    }
    
    if (!time_col %in% colnames(data)) {
      stop(paste("Column '", time_col, "' not found in the data frame.", sep = ""))
    }
    
    if (!status_col %in% colnames(data)) {
      stop(paste("Column '", status_col, "' not found in the data frame.", sep = ""))
    }
    colnames(data)[which(names(data) == time_col)] <- "time"
    colnames(data)[which(names(data) == status_col)] <- "status"
    data.df <- data
  } 
  else if (!is.null(object_stats)) {
    if (!inherits(object_stats, "Stat")) {
      stop("The 'object_stats' parameter must be an instance of class 'Stat'.")
    }
    
    clean.data <- ExtractCleanData(object_stats)
    info.data <- ExtractInfoData(object_stats)
    
    if (is.null(clean.data) || nrow(clean.data) == 0) {
      stop("Failed to extract valid clean data from the provided 'Stat' object.")
    }
    
    if (is.null(info.data) || nrow(info.data) == 0) {
      info.data <- data.frame(row.names = rownames(clean.data))
      cat("info.data was created from clean.data with", nrow(info.data), "rows.\n")
    }
    rownames_clean_data <- rownames(clean.data)
    
    info_filtered <- info.data[rownames(info.data) %in% rownames_clean_data, ]
    
    
    survival.data <- cbind(clean.data, info_filtered[, c(time_col, status_col)])
    colnames(survival.data)[which(names(survival.data) == time_col)] <- "time"
    colnames(survival.data)[which(names(survival.data) == status_col)] <- "status"
    data.df <- survival.data
    
    
  } 
  else {
    stop("At least one of 'data' and 'object_stats' must be provided.")
  }
  
  object_prognosix@filtered.set[["validation"]] <- data.df
  
  cat("The independent validation set has been added.\n")
  cat("Updating 'PrognosiX' object...\n")
  cat("The 'PrognosiX' object has been updated with the following slots:\n")
  cat("- 'filtered.set' slot updated.\n")
  
  return(object_prognosix)
}

evaluate_best_model_val <- function(object, 
                                    validation_data,
                                    save_plot = TRUE,
                                    save_dir = here("PrognosiX", "best_model_val"),
                                    plot_width = 5,
                                    plot_height = 5,
                                    base_size = 14) {
  model_name <- object@best.model[["model_name"]]
  cat("Selected model name: ", model_name, "\n") 
  validation_data<-object@filtered.set$validation
  training_data<-object@filtered.set$training
  
  model<-object@best.model[["model"]][["model"]]
  
  feature_names <- switch(class(model)[1],
                          "coxnet" = rownames(as.matrix(model$beta)),
                          "coxph" = names(model$coefficients),
                          "plsRcoxmodel" = rownames(as.matrix(model$Coeffs)),
                          "superpc" = names(model$feature.scores),
                          "rfsrc" = model$xvar.names,
                          "CoxBoost" = model[["xnames"]],
                          stop("Unsupported model type!")
  )
  
  important_vars<-feature_names
  missing_vars <- setdiff(important_vars, names(validation_data))
  prediction_data <- validation_data[, intersect(feature_names, colnames(validation_data)), drop = FALSE]
  for (var in missing_vars) {
    if (is.numeric(training_data[[var]])) {
      prediction_data[[var]] <- median(training_data[[var]], na.rm = TRUE)
    }
  }
  for (var in missing_vars) {
    if (is.factor(training_data[[var]])) {
      most_common <- names(which.max(table(training_data[[var]])))
      prediction_data[[var]] <- factor(most_common, 
                                       levels = levels(training_data[[var]]))
    }
  }
  
  validation_data<-prediction_data
  if (model_name == "lasso_model") {
    cat("Evaluating Lasso model...\n")  
    fit_best_cv_lasso <- slot(object, "survival.model")[["lasso_model"]][["model"]]
    best_lambda_value <- slot(object, "survival.model")[["lasso_model"]][["lambda"]]
    
    results_roc <- evaluate_lasso_roc(fit_best_cv_lasso, 
                                      best_lambda_value,
                                      validation_data =validation_data,
                                      save_plot = save_plot,
                                      save_dir = save_dir,
                                      plot_width = plot_width,
                                      plot_height = plot_height,
                                      base_size = base_size)
    
    results_km <- evaluate_lasso_km(
      fit_best_cv_lasso = fit_best_cv_lasso,
      best_lambda_value = best_lambda_value,
      data = validation_data,
      data_name = "val",
      save_plot = save_plot,
      save_dir = save_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size
    )
    
  } else if (model_name == "ridge_model") {
    cat("Evaluating Ridge model...\n") 
    fit_best_cv_ridge <- slot(object, "survival.model")[["ridge_model"]][["model"]]
    best_lambda_value <- slot(object, "survival.model")[["ridge_model"]][["lambda"]]
    
    results_roc <- evaluate_ridge_roc(fit_best_cv_ridge, 
                                      best_lambda_value,
                                      validation_data=validation_data,
                                      save_plot = save_plot,
                                      save_dir = save_dir,
                                      plot_width = plot_width,
                                      plot_height = plot_height,
                                      base_size = base_size)
    
    results_km<- evaluate_ridge_km(
      fit_best_cv_ridge = fit_best_cv_ridge,
      best_lambda_value = best_lambda_value,
      data = validation_data,
      data_name = "val",
      save_plot = save_plot,
      save_dir = save_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size
    )
    
  } else if (model_name == "pls_model") {
    cat("Evaluating PLS model...\n")  # 调试信息
    fit_best_cv_pls <- slot(object, "survival.model")[["pls_model"]][["model"]]
    best_lambda_value <- slot(object, "survival.model")[["pls_model"]][["lambda"]]
    
    results_roc <- evaluate_pls_roc(fit_best_cv_pls, 
                                    best_lambda_value,
                                    validation_data=validation_data,
                                    save_plot = save_plot,
                                    save_dir = save_dir,
                                    plot_width = plot_width,
                                    plot_height = plot_height,
                                    base_size = base_size)
    
    results_km <- evaluate_pls_km(
      fit_best_cv_pls = fit_best_cv_pls,
      data = validation_data,
      data_name = "val",
      save_plot = save_plot,
      save_dir = save_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size
    )
    
  } else if (model_name == "coxboost_model") {
    cat("Evaluating Coxboost model...\n")  # 调试信息
    fit_best_cv_coxboost <- slot(object, "survival.model")[["coxboost_model"]][["model"]]
    
    results_roc <- evaluate_coxboost_roc(fit_best_cv_coxboost, 
                                         validation_data=validation_data,
                                         save_plot = save_plot,
                                         save_dir = save_dir,
                                         plot_width = plot_width,
                                         plot_height = plot_height,
                                         base_size = base_size)
    
    results_km <- evaluate_coxboost_km(
      fit_best_cv_coxboost = fit_best_cv_coxboost,
      data = validation_data,
      data_name = "val",
      save_plot = save_plot,
      save_dir = save_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size
    )
    
  } else if (model_name == "coxph_model") {
    
    cat("Evaluating CoxPH model...\n") 
    fit_coxph_model <- slot(object, "survival.model")[["coxph_model"]][["model"]]
    
    results_roc <- evaluate_coxph_roc(fit_coxph_model,
                                      validation_data=validation_data,
                                      save_plot = save_plot,
                                      save_dir = save_dir,
                                      plot_width = plot_width,
                                      plot_height = plot_height,
                                      base_size = base_size)
    
    results_km <- evaluate_coxph_km(
      fit_best_cv_coxph = fit_coxph_model,
      data = validation_data,
      data_name = "val",
      save_plot = save_plot,
      save_dir = save_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size
    )
    
  } else if (model_name == "superpc_model") {
    cat("Evaluating SuperPC model...\n")  
    fit_superpc_model <- slot(object, "survival.model")[["superpc_model"]][["model"]]
    threshold <- slot(object, "survival.model")[["superpc_model"]][["threshold"]]
    
    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$training
    
    results_roc <- evaluate_superpc_roc(fit_best_cv_superpc = fit_superpc_model, 
                                        threshold = threshold,
                                        train_data = train_data,
                                        validation_data=validation_data,
                                        save_plot = save_plot,
                                        save_dir = save_dir,
                                        plot_width = plot_width,
                                        plot_height = plot_height,
                                        base_size = base_size)
    
    results_km <- evaluate_superpc_km(
      fit_best_cv_superpc = fit_superpc_model,
      data = validation_data,
      data_name = "val",
      save_plot = save_plot,
      save_dir = save_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size
    )
    
  } else if (model_name == "rsf_model") {
    cat("Evaluating RSF model...\n")  
    fit_rsf_model <- slot(object, "survival.model")[["rsf_model"]][["model"]]
    
    results_roc <- evaluate_rsf_roc(fit_rsf_model, 
                                    validation_data=validation_data,
                                    save_plot = save_plot,
                                    save_dir = save_dir,
                                    plot_width = plot_width,
                                    plot_height = plot_height,
                                    base_size = base_size)
    
    results_km <- evaluate_rsf_km(
      fit_best_rsf = fit_rsf_model,
      data = validation_data,
      data_name = "val",
      save_plot = save_plot,
      save_dir = save_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      base_size = base_size
    )
    
  } else {
    stop("Unsupported model name!")
  }
  if (inherits(object, 'PrognosiX')) {
    
    object@best.model[["val_results"]] <- list(results_roc=results_roc,
                                               results_km=results_km)
    cat("Model evaluation completed successfully.\n") 
    return(object)
  }
  return(object)
}




PrognosiXClinicalPrediction <- function(
    object,  
    new_data, 
    group_col = "group", 
    palette_name = "Royal1",  
    save_dir = here("PrognosiX", "clinical_predictions"),  
    plot_width = 6,  
    plot_height = 6,  
    alpha = 1,  
    base_size = 14 
) {
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  if (is.null(new_data)) {
    stop("New clinical data is missing.")
  }
  
  best_model <- object@best.model[["model"]][["model"]]
  training_data<-object@filtered.set$training
  cat("Extracting feature names based on model type...\n")
  feature_names <- switch(class(best_model)[1],
                          "coxnet" = rownames(as.matrix(best_model$beta)),
                          "coxph" = names(best_model$coefficients),
                          "plsRcoxmodel" = rownames(as.matrix(best_model$Coeffs)),
                          "superpc" = names(best_model$feature.scores),
                          "rfsrc" = best_model$xvar.names,
                          "CoxBoost" = best_model[["xnames"]],
                          stop("Unsupported model type!")
  )
  
  important_vars<-feature_names
  prediction_data <- new_data[, intersect(feature_names, colnames(new_data)), drop = FALSE]
  missing_vars <- setdiff(important_vars, names(prediction_data))

  for (var in missing_vars) {
    if (is.numeric(training_data[[var]])) {
      prediction_data[[var]] <- median(training_data[[var]], na.rm = TRUE)
    }
  }
  for (var in missing_vars) {
    if (is.factor(training_data[[var]])) {
      most_common <- names(which.max(table(training_data[[var]])))
      prediction_data[[var]] <- factor(most_common, 
                                       levels = levels(training_data[[var]]))
    }
  }
  
  cat("Calculating risk scores based on model type...\n")
  risk_scores <- switch(class(best_model)[1],
                        "coxnet" = predict(best_model, newx = as.matrix(prediction_data), type = "response"),
                        "coxph" = predict(best_model, newdata = prediction_data, type = "risk"),
                        "superpc" = {
                          x_new_superpc <- t(as.matrix(prediction_data))
                          x_new_superpc <- apply(x_new_superpc, 2, as.numeric)
                          
                          test <- list(
                            x = x_new_superpc,
                            y = rep(0, nrow(prediction_data)),  # 伪时间
                            censoring.status = rep(1, nrow(prediction_data)),  # 伪状态
                            featurenames = colnames(prediction_data)
                          )
                          
                          risk_score <- superpc.predict(best_model,
                                                        data = test,
                                                        newdata = test,
                                                        threshold = 0.5,
                                                        n.components = 1)
                          as.numeric(risk_score$v.pred)
                        },
                        "rfsrc" = predict(best_model, newdata = prediction_data, type = "risk")$predicted,
                        "plsRcoxmodel" = predict(best_model, newdata = as.matrix(prediction_data), type = "risk"),
                        "CoxBoost" = predict(best_model, newdata = as.matrix(prediction_data), type = "lp"),
                        stop("Unsupported model type!")
  )
  
  cat("Calculating risk groups (High/Low) based on median risk score...\n")
  median_risk <- median(risk_scores, na.rm = TRUE)
  risk_group <- ifelse(risk_scores > median_risk, "High", "Low")
  
  cat("Applying log transformation to risk scores...\n")
  risk_scores_log <- log(abs(risk_scores) + 1) * sign(risk_scores)
  
  final_result <- data.frame(
    Sample = rownames(prediction_data),
    RiskGroup = as.factor(risk_group),
    RiskScore = round(risk_scores_log, 3)
  )
  risk_groups <- unique(final_result$RiskGroup)
  
  colors <- wes_palette(n = length(risk_groups), name = palette_name, type = "discrete")
  names(colors) <- risk_groups 
  
  
  p <- ggplot(final_result, aes(x = RiskGroup, y = RiskScore, fill = RiskGroup)) +
    geom_boxplot(outlier.shape = 19, outlier.colour = "black", outlier.size = 1) +  
    geom_jitter(width = 0.2, size = 2, aes(color = RiskGroup), alpha = 0.6) +  
    scale_fill_manual(values = colors) +  
    scale_color_manual(values = colors) +  
    labs(
      title = "Visualization of Predicted Group and Probabilities",
      x = "Group",
      y = "Predicted Probability"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 2),
      axis.text.y = element_text(size = base_size - 2),
      axis.title.x = element_text(size = base_size),
      axis.title.y = element_text(size = base_size),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 2)
    ) +
    theme_prism(base_size = base_size)  
  
  print(p)
  
  
  ggsave(file.path(save_dir, "risk_group_prediction.pdf"), p, width = plot_width, height = plot_height,device = "pdf")
  cat("Risk group prediction visualization saved to:", file.path(save_dir, "risk_group_prediction.pdf"), "\n")
  
  return(list(
    predictions = final_result,
    median_risk = median_risk
  ))
}
