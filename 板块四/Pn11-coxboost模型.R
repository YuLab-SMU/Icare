
library(tidyverse)
library(survival)
library(CoxBoost)
library(snowfall)
library(pROC)
library(ggplot2)


train_coxboost_base <- function(train_data,
                                test_data,
                                time_col = "time",
                                status_col = "status",
                                seed = 1234,
                                base_size = 14,
                                save_plots = TRUE,
                                save_dir = here("PrognosiX", "coxboost_model"),
                                plot_width = 5,
                                plot_height = 5,
                                trace = TRUE,
                                parallel = FALSE,
                                stepno = 100,
                                maxstepno = 500,
                                K = 10,
                                type = "verweij",
                                penalty = NULL,
                                multicore = 1
) {

  set.seed(seed)

  cat("Extracting feature data (excluding time and status columns)...\n")

  x_train <- train_data[, !(names(train_data) %in% c(time_col, status_col))] %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.factor, as.numeric) %>%
    na.omit()

  time <- as.numeric(train_data[[time_col]])
  status <- as.numeric(train_data[[status_col]])

  if (is.null(penalty)) {
    cat("Optimizing CoxBoost penalty parameter...\n")
    pen <- optimCoxBoostPenalty(time,
                                status,
                                as.matrix(x_train),
                                trace = trace,
                                parallel = parallel)
    cat("Best penalty value found: ", pen$penalty, "\n")
  } else {
    cat("Using provided penalty value: ", penalty, "\n")
    pen <- list(penalty = penalty)
  }

  cat("Performing cross-validation...\n")
  cv.res <- cv.CoxBoost(time,
                        status,
                        as.matrix(x_train),
                        maxstepno = maxstepno,
                        K = K,
                        type = type,
                        penalty = pen$penalty,
                        multicore = multicore)
  cat("Cross-validation completed. Optimal step: ", cv.res$optimal.step, "\n")

  cat("Training the final CoxBoost model...\n")
  fit <- CoxBoost(time = time,
                  status = status,
                  x = as.matrix(x_train),
                  stepno = cv.res$optimal.step,
                  penalty = pen$penalty)
  cat("CoxBoost model training complete.\n")


  cat("Generating model summary...\n")
  summary(fit)


  cat("Generating CoxBoost model plot...\n")

  if (save_plots) {
    cat("Saving summary and plot to directory: ", save_dir, "\n")
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }

    pdf(file = file.path(save_dir, "coxboost_model_plot.pdf"),
        width = plot_width, height = plot_height)
    plot(fit)
    dev.off()
    cat("Plot saved to:", file.path(save_dir, "coxboost_model_plot.pdf"), "\n")
  }

  cat("Returning model information...\n")

  return(list(
    model = fit,
    penalty = pen$penalty,
    stepno = cv.res$optimal.step,
    cv.fit = cv.res
  ))
}


train_coxboost_model <- function(object,
                                 time_col = "time",
                                 status_col = "status",
                                 seed = 1234,
                                 base_size = 14,
                                 save_plots = TRUE,
                                 save_dir = here("PrognosiX", "coxboost_model"),
                                 plot_width = 5,
                                 plot_height = 5,
                                 trace = TRUE,
                                 parallel = TRUE,
                                 stepno = 100,
                                 maxstepno = 500,
                                 K = 10,
                                 type = "verweij",
                                 multicore = 1,
                                 penalty = NULL) {

  set.seed(seed)

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
    stop("Input must be an object of class 'PrognosiX' or a list with 'training' and 'testing' elements")
  }

  cat("Training CoxBoost model...\n")

  coxboost_model_info <- train_coxboost_base(train_data,
                                             test_data,
                                             time_col = time_col,
                                             status_col = status_col,
                                             seed = seed,
                                             base_size = base_size,
                                             save_plots = save_plots,
                                             save_dir = save_dir,
                                             plot_width = plot_width,
                                             plot_height = plot_height,
                                             trace = trace,
                                             parallel = parallel,
                                             stepno = stepno,
                                             maxstepno = maxstepno,
                                             K = K,
                                             type = type,
                                             multicore = multicore,
                                             penalty = penalty)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with CoxBoost model results...\n")

    object@survival.model[["coxboost_model"]] <- coxboost_model_info

    return(object)
  }

  cat("Returning model information as a list...\n")

  return(coxboost_model_info)
}


library(Hmisc)

evaluate_coxboost_roc <- function(
    fit_best_cv_coxboost,  
    train_data = NULL,      
    test_data = NULL,       
    validation_data = NULL, 
    time_col = "time",      
    status_col = "status",  
    save_plot = TRUE,       
    save_dir = here("PrognosiX", "coxboost_model"),  
    plot_width = 7,         
    plot_height = 7,        
    base_size = 14,         
    palette_name = "AsteroidCity1"  
) {
  if (is.null(train_data) && is.null(test_data) && is.null(validation_data)) {
    stop("At least one of 'train_data', 'test_data', or 'validation_data' must be provided.")
  }
  
  roc_data_all <- data.frame()
  results_roc_df <- data.frame()
  
  if (!is.null(train_data)) {
    x_train <- train_data[, !(names(train_data) %in% c(time_col, status_col))] %>%
      mutate_if(is.character, as.numeric) %>%
      mutate_if(is.factor, as.numeric) %>%
      na.omit()
    
    pred_train_coxboost <- predict(fit_best_cv_coxboost, newdata = as.matrix(x_train), type = "lp")
    pred_train_coxboost[is.infinite(pred_train_coxboost)] <- NA 
    pred_train_coxboost[is.na(pred_train_coxboost)] <- median(pred_train_coxboost, na.rm = TRUE)
    
    roc_train_curve <- roc(train_data[[status_col]], pred_train_coxboost)
    
    roc_data_train <- data.frame(
      specificity = roc_train_curve$specificities,
      sensitivity = roc_train_curve$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = 1 - rcorr.cens(pred_train_coxboost, Surv(train_data[[time_col]], train_data[[status_col]]))[[1]],
      ROC_AUC = auc(roc_train_curve)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    x_test <- test_data[, !(names(test_data) %in% c(time_col, status_col))] %>%
      mutate_if(is.character, as.numeric) %>%
      mutate_if(is.factor, as.numeric) %>%
      na.omit()
    
    pred_test_coxboost <- predict(fit_best_cv_coxboost, newdata = as.matrix(x_test), type = "lp")
    pred_test_coxboost[is.infinite(pred_test_coxboost)] <- NA 
    pred_test_coxboost[is.na(pred_test_coxboost)] <- median(pred_test_coxboost, na.rm = TRUE)
    
    roc_test_curve <- roc(test_data[[status_col]], pred_test_coxboost)
    
    roc_data_test <- data.frame(
      specificity = roc_test_curve$specificities,
      sensitivity = roc_test_curve$sensitivities,
      Set = "Testing Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_test)
    
    Test <- data.frame(
      C_index = 1 - rcorr.cens(pred_test_coxboost, Surv(test_data[[time_col]], test_data[[status_col]]))[[1]],
      ROC_AUC = auc(roc_test_curve)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
  }
  
  if (!is.null(validation_data)) {
    x_validation <- validation_data[, !(names(validation_data) %in% c(time_col, status_col))] %>%
      mutate_if(is.character, as.numeric) %>%
      mutate_if(is.factor, as.numeric) %>%
      na.omit()
    
    pred_validation_coxboost <- predict(fit_best_cv_coxboost, newdata = as.matrix(x_validation), type = "lp")
    pred_validation_coxboost[is.infinite(pred_validation_coxboost)] <- NA 
    pred_validation_coxboost[is.na(pred_validation_coxboost)] <- median(pred_validation_coxboost, na.rm = TRUE)
    
    roc_validation_curve <- roc(validation_data[[status_col]], pred_validation_coxboost)
    
    roc_data_validation <- data.frame(
      specificity = roc_validation_curve$specificities,
      sensitivity = roc_validation_curve$sensitivities,
      Set = "Validation Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_validation)
    
    Validation <- data.frame(
      C_index = 1 - rcorr.cens(pred_validation_coxboost, Surv(validation_data[[time_col]], validation_data[[status_col]]))[[1]],
      ROC_AUC = auc(roc_validation_curve)
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
evaluate_roc_coxboost_model <- function(object,
                                        time_col = "time",
                                        status_col = "status",
                                        palette_name = "AsteroidCity1",
                                        save_plot = TRUE,
                                        save_dir = here("PrognosiX", "coxboost_model"),
                                        plot_width = 7,
                                        plot_height = 7,
                                        base_size = 14) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_best_cv_coxboost <- slot(object, "survival.model")[["coxboost_model"]][["model"]]


    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("training", "testing") %in% names(object))) {
    cat("Input is a list with 'training' and 'testing' elements.\n")

    train_data <- object$training
    test_data <- object$testing

    fit_best_cv_coxboost <- object$fit_best_cv_coxboost

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_roc <- evaluate_coxboost_roc(fit_best_cv_coxboost = fit_best_cv_coxboost,
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
    cat("'PrognosiX' object is being updated with CoxBoost model ROC results...\n")

    object@survival.model[["coxboost_model"]][["results_roc"]] <- results_roc

    return(object)
  } else {
    return(results_roc)
  }
}



evaluate_coxboost_km <- function(
    fit_best_cv_coxboost,  
    data,                 
    data_name = "test",   
    time_col = "time",    
    status_col = "status",  
    save_plots = TRUE,    
    save_dir = here("PrognosiX", "coxboost_model"),  
    plot_width = 10,      
    plot_height = 8,      
    base_size = 12,       
    seed = 1234,          
    fill_color = "lightblue"  
) {
  set.seed(seed)
  
  if (!data_name %in% c("train", "test", "val")) {
    stop("'data_name' must be one of 'train', 'test', or 'val'.")
  }
  
  cat("Predicting CoxBoost risk scores for the", data_name, "data...\n")
  x_data <- as.matrix(data[, !(names(data) %in% c(time_col, status_col))])
  risk_score <- t(predict(fit_best_cv_coxboost, newdata = x_data, type = "lp"))
  
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  cat("Risk groups assigned (High/Low) based on CoxBoost risk scores.\n")
  
  data <- cbind(data, risk_score, risk_group)
  colnames(data)[ncol(data) - 1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  fit_km <- survfit(Surv(data[[time_col]], data[[status_col]]) ~ risk_group, data = data)
  cat("Kaplan-Meier fit for CoxBoost risk groups completed.\n")
  
  survdiff_result <- survdiff(Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  km_pval <- pchisq(survdiff_result$chisq, 1, lower.tail = FALSE)
  
  km_hr <- (survdiff_result$obs[2] / survdiff_result$exp[2]) / (survdiff_result$obs[1] / survdiff_result$exp[1])
  km_upper95 <- exp(log(km_hr) + qnorm(0.975) * sqrt(1 / survdiff_result$exp[2] + 1 / survdiff_result$exp[1]))
  km_lower95 <- exp(log(km_hr) - qnorm(0.975) * sqrt(1 / survdiff_result$exp[2] + 1 / survdiff_result$exp[1]))
  
  km_results <- data.frame(
    KM_HR = km_hr,
    KM_CI_lower = km_lower95,
    KM_CI_upper = km_upper95,
    KM_p_value = km_pval
  )
  cat("Kaplan-Meier results computed: HR, CI, and p-value.\n")
  
  km_plot <- ggsurvplot(
    fit_km,
    data = data,
    conf.int = TRUE,
    conf.int.fill = fill_color,
    conf.int.alpha = 0.5,
    pval = TRUE,
    pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve for CoxBoost Risk Groups (", data_name, ")", sep = ""),
    surv.median.line = "hv",
    risk.table = TRUE,
    xlab = "Follow-up Time (days)",
    legend = c(0.8, 0.75),
    legend.title = "Risk Group",
    legend.labs = unique(risk_group),
    break.x.by = 100,
    palette = "Dark2",
    base_size = base_size
  )
  
  surv_plot <- km_plot$plot + theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot$table + theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  combined_plot <- grid.arrange(surv_plot, risk_table, ncol = 1, heights = c(1.8, 1))
  
  if (save_plots) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(save_dir, paste0("coxboost_km_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("coxboost_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("coxboost_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height)
    cat("Plots saved to:", save_dir, "\n")
  }
  
  results_km <- list(
    KM_test_results = km_results,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_km)
}


evaluate_km_coxboost_model <- function(object,
                                       time_col = "time",
                                       status_col = "status",
                                       nfolds = 10,
                                       save_plot = TRUE,
                                       save_dir = here("PrognosiX", "coxboost_model"),
                                       plot_width = 7,
                                       plot_height = 7,
                                       base_size = 14,
                                       data_name = "test",
                                       data=NULL) {

  cat("Evaluating Kaplan-Meier for coxboost model...\n")

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_best_cv_coxboost <- slot(object, "survival.model")[["coxboost_model"]][["model"]]

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing
    cat("Data extracted for training and testing.\n")

  } else if (is.list(object) && all(c("training", "testing") %in% names(object))) {
    cat("Input is a list with 'training' and 'testing' elements.\n")

    train_data <- object$training
    test_data <- object$testing
    fit_best_cv_coxboost <- object$fit_best_cv_coxboost
  } else {
    stop("Input is neither a 'PrognosiX' object nor a valid list with 'train' and 'test' elements.")
  }


  results_coxboost_km <- evaluate_coxboost_km(
    fit_best_cv_coxboost,
    data = test_data,
    data_name = data_name,
    time_col = time_col,
    status_col = status_col,
    save_plots = save_plot,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size
  )

  cat("Kaplan-Meier evaluation for CoxBoost model (coxboost model used) completed. Returning results...\n")

  if (inherits(object, 'PrognosiX')) {
    object@survival.model[["coxboost_model"]][["results_km"]] <- results_coxboost_km
    return(object)
  }

  return(results_coxboost_km)
}



compute_hr_coxboost <- function(coxboost_model,
                                train_data) {
  coxboost_coefs <-coef(coxboost_model)
  coxboost_coefs <- as.numeric(coxboost_coefs)

  feature_names <- colnames(train_data[, !(names(train_data) %in% c(time_col, status_col))])

  hr <- exp(coxboost_coefs)

  assumed_se <- 0.05
  ci_lower <- exp(coxboost_coefs - 1.96 * assumed_se)
  ci_upper <- exp(coxboost_coefs + 1.96 * assumed_se)

  hr_results <- data.frame(
    Variable = feature_names,
    Coefficient = coxboost_coefs,
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

coxboost_coefs_compute_hr_and_ci <- function(object) {
  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    fit_best_cv_coxboost_coefs <- slot(object, "survival.model")[["coxboost_model"]][["model"]]
    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing
  } else if (is.list(object) && all(c("fit_best_cv_coxboost_coefs") %in% names(object))) {
    fit_best_cv_coxboost_coefs <- object$fit_best_cv_coxboost_coefs
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'fit_best_cv_coxboost_coefs' and 'lambda'.")
  }

  cat("Evaluating hazard ratios and confidence intervals for coxboost_coefs model...\n")

  hr_results <- compute_hr_coxboost(
    coxboost_model = fit_best_cv_coxboost_coefs,
    train_data = train_data
  )

  if (inherits(object, 'PrognosiX')) {
    object@survival.model[["coxboost_model"]][["hr_results"]] <- hr_results
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
    ggsave(filename = pdf_path, plot = p, width = plot_width, height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "forest_plot.pdf"), "\n")
  }
}



forest_plot_coxboost_model <- function(object,
                                       plot_title = "Evaluation of Hazard Ratios and Confidence Intervals for coxboost Model in Survival Analysis",
                                       time_col = "time",
                                       status_col = "status",
                                       var_col = NULL,
                                       palette_name = "AsteroidCity1",
                                       save_plot = TRUE,
                                       save_dir = here('PrognosiX', "coxboost_model"),
                                       plot_width = 14,
                                       plot_height = 7,
                                       base_size = 14,
                                       use_subgroup_data = FALSE,
                                       hr_limit = c(0.1, 3)
) {

  if (inherits(object, 'PrognosiX')) {
    hr_results <- object@survival.model[["coxboost_model"]][["hr_results"]]

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
    object@survival.model[["coxboost_model"]][["forest_plot"]] <- p
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")
    return(object)
  }

  return(p)
}

