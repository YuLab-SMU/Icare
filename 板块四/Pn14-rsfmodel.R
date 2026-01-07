library(randomForestSRC)
library(survival)
library(randomForestSRC)
library(survival)

train_rsf_base <- function(train_data,
                           time_col = "time",
                           status_col = "status",
                           nfolds = 10,
                           ntree = 1000,
                           nodesize = 10,
                           splitrule = "logrank",
                           seed = 1234) {

  cat("Starting train_rsf_base function...\n")
  cat("Parameters: time_col =", time_col, ", status_col =", status_col, "\n")

  if (!(time_col %in% names(train_data)) || !(status_col %in% names(train_data))) {
    stop("Specified time or status column does not exist in the data.")
  }

  cat("Setting up cross-validation with", nfolds, "folds...\n")
  folds <- sample(rep(1:nfolds, length.out = nrow(train_data)))

  rsf_models <- list()
  for (i in 1:nfolds) {
    cat("Training fold", i, "...\n")

    train_fold <- train_data[folds != i, ]

    formula <- as.formula(paste("Surv(", time_col, ",", status_col, ") ~ ."))
    cat("Using formula:", deparse(formula), "\n")

    rsf_models[[i]] <- rfsrc(formula,
                             data = train_fold,
                             ntree = ntree,
                             nodesize = nodesize,
                             splitrule = splitrule,
                             importance = TRUE,
                             proximity = TRUE,
                             forest = TRUE,
                             seed = seed)

    cat("Fold", i, "completed.\n")
  }

  cat("Selecting the best model based on OOB errors...\n")
  oob_errors <- sapply(rsf_models, function(model) model$err.rate[length(model$err.rate)])
  best_model <- rsf_models[[which.min(oob_errors)]]

  cat("Best model selected from fold:", which.min(oob_errors), "\n")

  return(list(model = best_model,
              oob_errors = oob_errors,
              best_fold = which.min(oob_errors)))
}


train_rsf_model <- function(object,
                            time_col = "time",
                            status_col = "status",
                            nfolds = 10,
                            ntree = 1000,
                            nodesize = 10,
                            splitrule = "logrank",
                            seed = 1234) {

  cat("Starting train_rsf_model function...\n")

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")

    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$train
    test_data <- data_sets$test

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with train and test datasets...\n")
    train_data <- object$train
    test_data <- object$test

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Training Random Survival Forest model...\n")

  rsf_model_info <- train_rsf_base(train_data, time_col, status_col, nfolds, ntree, nodesize, splitrule, seed)

  best_model <- rsf_model_info$model
  best_fold <- rsf_model_info$best_fold
  oob_errors <- rsf_model_info$oob_errors

  cat("Best model selected from fold:", best_fold, "\n")

  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object with Random Survival Forest model results...\n")

    object@survival.model[["rsf_model"]] <- list(model = best_model,
                                                 best_fold = best_fold,
                                                 oob_errors = oob_errors)

    cat("'PrognosiX' object updated successfully.\n")
    return(object)
  }

  cat("Returning model information as a list...\n")

  return(list(model = best_model,
              best_fold = best_fold,
              oob_errors = oob_errors))
}


evaluate_rsf_roc <- function(
    fit_rsf_model,          
    train_data = NULL,      
    test_data = NULL,       
    validation_data = NULL, 
    time_col = "time",      
    status_col = "status",  
    save_plot = TRUE,       
    save_dir = here("PrognosiX", "rsf_model"),  
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
    rsf_risk_score_train <- predict(fit_rsf_model, newdata = train_data)$predicted
    rsf_risk_score_train[is.infinite(rsf_risk_score_train)] <- NA 
    rsf_risk_score_train[is.na(rsf_risk_score_train)] <- median(rsf_risk_score_train, na.rm = TRUE)
    
    rsf_roc_train <- roc(train_data[[status_col]], rsf_risk_score_train)
    
    roc_data_train <- data.frame(
      specificity = rsf_roc_train$specificities,
      sensitivity = rsf_roc_train$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = concordance.index(
        x = rsf_risk_score_train,
        surv.time = train_data[[time_col]],
        surv.event = train_data[[status_col]],
        method = "noether"
      )$c.index,
      ROC_AUC = auc(rsf_roc_train)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    rsf_risk_score_test <- predict(fit_rsf_model, newdata = test_data)$predicted
    rsf_risk_score_test[is.infinite(rsf_risk_score_test)] <- NA 
    rsf_risk_score_test[is.na(rsf_risk_score_test)] <- median(rsf_risk_score_test, na.rm = TRUE)
    rsf_roc_test <- roc(test_data[[status_col]], rsf_risk_score_test)
    
    roc_data_test <- data.frame(
      specificity = rsf_roc_test$specificities,
      sensitivity = rsf_roc_test$sensitivities,
      Set = "Testing Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_test)
    
    Test <- data.frame(
      C_index = concordance.index(
        x = rsf_risk_score_test,
        surv.time = test_data[[time_col]],
        surv.event = test_data[[status_col]],
        method = "noether")$c.index,
      ROC_AUC = auc(rsf_roc_test)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
  }
  
  if (!is.null(validation_data)) {
    rsf_risk_score_validation <- predict(fit_rsf_model, newdata = validation_data)$predicted
    rsf_risk_score_validation[is.infinite(rsf_risk_score_validation)] <- NA 
    rsf_risk_score_validation[is.na(rsf_risk_score_validation)] <- median(rsf_risk_score_validation, na.rm = TRUE)
    
    rsf_roc_validation <- roc(validation_data[[status_col]], rsf_risk_score_validation)
    
    roc_data_validation <- data.frame(
      specificity = rsf_roc_validation$specificities,
      sensitivity = rsf_roc_validation$sensitivities,
      Set = "Validation Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_validation)
    
    Validation <- data.frame(
      C_index = concordance.index(
        x = rsf_risk_score_validation,
        surv.time = validation_data[[time_col]],
        surv.event = validation_data[[status_col]],
        method = "noether"
      )$c.index,
      ROC_AUC = auc(rsf_roc_validation)
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



evaluate_roc_rsf_model <- function(object,
                                   time_col = "time",
                                   status_col = "status",
                                   save_plot = TRUE,
                                   save_dir = here("PrognosiX", "rsf_model"),
                                   plot_width = 7,
                                   plot_height = 7,
                                   base_size = 14,
                                   palette_name = "AsteroidCity1") {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_rsf_model <- slot(object, "survival.model")[["rsf_model"]][["model"]]

    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$train
    test_data <- data_sets$test

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    train_data <- object$train
    test_data <- object$test

    fit_rsf_model <- object$fit_rsf_model
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_roc <- evaluate_rsf_roc(fit_rsf_model,
                                  train_data, test_data,
                                  time_col = time_col, status_col = status_col,
                                  save_plot = save_plot, save_dir = save_dir,
                                  plot_width = plot_width, plot_height = plot_height,
                                  base_size = base_size, palette_name = palette_name)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is updated with new ROC results.\n")
    object@survival.model[["rsf_model"]][["results_roc"]] <- results_roc
    return(object)
  }

  return(results_roc)
}


evaluate_rsf_km <- function(
    fit_best_rsf, data, data_name = "test",
    time_col = "time", status_col = "status",
    save_plot = TRUE, save_dir = here("PrognosiX", "rsf_model"),
    palette_name = "Dark2", plot_width = 10, plot_height = 8,
    base_size = 12, seed = 1234, fill_color = "lightblue"
) {
  set.seed(seed)
  
  if (!data_name %in% c("train", "test", "val")) {
    stop("'data_name' must be one of 'train', 'test', or 'val'.")
  }
  
  data_x <- as.data.frame(data)
  
  risk_score <- predict(fit_best_rsf, newdata = data_x, type = "risk")
  cat("Risk scores calculated for", data_name, "data.\n")
  
  risk_score <- as.vector(risk_score$predicted)
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  data <- cbind(data, risk_score, risk_group)
  
  colnames(data)[ncol(data) - 1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  fit_km_rsf <- survfit(Surv(data[[time_col]], data[[status_col]]) ~ risk_group, data = data)
  cat("Kaplan-Meier model fitted for", data_name, "data.\n")
  
  km_pval_rsf <- pchisq(survdiff(Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)$chisq, 1, lower.tail = FALSE)
  
  cox_rsf <- coxph(Surv(data[[time_col]], data[[status_col]]) ~ risk_group, data = data)
  km_hr_rsf <- exp(coef(cox_rsf))
  km_ci_rsf <- confint(cox_rsf)
  
  km_results_rsf <- data.frame(
    KM_HR = km_hr_rsf,
    KM_CI_lower = km_ci_rsf[1],
    KM_CI_upper = km_ci_rsf[2],
    KM_p_value = km_pval_rsf
  )
  cat("Kaplan-Meier results stored in a dataframe for", data_name, "data.\n")
  
  km_plot_rsf <- ggsurvplot(
    fit_km_rsf, data = data, conf.int = TRUE,
    conf.int.fill = fill_color, conf.int.alpha = 0.5,
    pval = TRUE, pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve for RSF Risk Groups (", data_name, ")", sep = ""),
    surv.median.line = "hv", risk.table = TRUE,
    xlab = "Follow-up Time (days)", legend = c(0.8, 0.75),
    legend.title = "Risk Group", legend.labs = unique(data$risk_group),
    break.x.by = 100, palette = palette_name, base_size = base_size
  )
  print(km_plot_rsf)
  
  cat("Kaplan-Meier plot generated for", data_name, "data.\n")
  
  surv_plot <- km_plot_rsf$plot + theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot_rsf$table + theme_prism(base_size = base_size) +
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
    
    ggsave(file.path(save_dir, paste0("rsf_time_distribution_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("rsf_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("rsf_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", save_dir, "\n")
    }
  
  results_rsf_km <- list(
    KM_test_results = km_results_rsf,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_rsf_km)
}




evaluate_km_rsf_model <- function(object,
                                  time_col = "time",
                                  status_col = "status",
                                  save_plot = TRUE,
                                  save_dir = here("PrognosiX", "rsf_model"),
                                  plot_width = 7,
                                  plot_height = 7,
                                  base_size = 14,
                                  data_name = "test",
                                  data=NULL) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")
    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_best_rsf <- slot(object, "survival.model")[["rsf_model"]][["model"]]
    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$train
    test_data <- data_sets$test
  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")
    train_data <- object$train
    test_data <- object$test
    fit_best_rsf <- object$fit_best_rsf
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_rsf_km <- evaluate_rsf_km(
    fit_best_rsf,
    data = test_data,
    data_name = data_name,
    time_col = time_col,
    status_col = status_col,
    save_plot = save_plot,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size
  )

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with rsf model Kaplan-Meier results...\n")
    object@survival.model[["rsf_model"]][["results_km"]]  <- results_rsf_km
    return(object)
  }

  return(results_rsf_km)
}


compute_hr_rsf <- function(rsf_model) {

  variable_importance <- randomForestSRC::vimp(rsf_model)

  feature_names <-variable_importance[["xvar.names"]]


  hr <- exp(variable_importance$importance)

  assumed_se <- 0.05
  ci_lower <- exp(variable_importance$importance - 1.96 * assumed_se)
  ci_upper <- exp(variable_importance$importance + 1.96 * assumed_se)

  hr_results <- data.frame(
    Variable = feature_names,
    Coefficient = variable_importance$importance,
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

rsf_compute_hr_and_ci <- function(object,
                                  time_col = "time",
                                  status_col = "status") {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    fit_best_cv_rsf <- slot(object, "survival.model")[["rsf_model"]][["model"]]

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {

    fit_best_cv_rsf <- object$fit_best_cv_rsf

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  cat("Evaluating hazard ratios and confidence intervals for RSF model...\n")

  hr_results <- compute_hr_rsf(rsf_model = fit_best_cv_rsf)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with RSF model results...\n")

    object@survival.model[["rsf_model"]][["hr_results"]] <- hr_results
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")

    return(object)
  }

  cat("Returning HR results as a data frame...\n")

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
    scale_color_manual(values =wes_palette(palette_name)) +
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



forest_plot_rsf_model <- function(object,
                                  plot_title = "Evaluation of Hazard Ratios and Confidence Intervals for rsf Model in Survival Analysis",
                                  time_col = "time",
                                  status_col = "status",
                                  var_col = NULL,
                                  palette_name = "AsteroidCity1",
                                  save_plot = TRUE,
                                  save_dir = here('PrognosiX', "rsf_model"),
                                  plot_width = 14,
                                  plot_height = 7,
                                  base_size = 14,
                                  use_subgroup_data = FALSE,
                                  hr_limit = c(0.1, 3)
) {

  if (inherits(object, 'PrognosiX')) {
    hr_results <- object@survival.model[["rsf_model"]][["hr_results"]]

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
    object@survival.model[["rsf_model"]][["forest_plot"]] <- p
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'survival.model' slot updated.\n")
    return(object)
  }

  return(p)
}



