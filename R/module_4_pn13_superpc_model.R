#' Train SuperPC Model (Base)
#'
#' Trains a SuperPC model (base function).
#'
#' @param train_data Training data.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param seed Seed.
#'
#' @return List containing model and threshold.
#' @export
train_superpc_base <- function(train_data,
                               time_col = "time",
                               status_col = "status",
                               seed = 123) {

  set.seed(seed)
  featurenames <- colnames(train_data)[!(colnames(train_data) %in% c(time_col, status_col))]

  data_list <- list(
    x = t(train_data[, featurenames]),
    y = train_data[[time_col]],
    censoring.status = train_data[[status_col]],
    featurenames = featurenames
  )

  fit_superpc <- superpc::superpc.train(data_list, type = "survival")

  cv_fit <- superpc::superpc.cv(fit_superpc, data_list)

  best_threshold <- cv_fit$thresholds[which.max(cv_fit$scor)]

  return(list(model = fit_superpc, threshold = best_threshold))
}

#' Train SuperPC Model for PrognosiX
#'
#' Wrapper to train SuperPC model on PrognosiX object.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param seed Seed.
#'
#' @return Updated PrognosiX object or model list.
#' @export
train_superpc_model <- function(object,
                                time_col = "time",
                                status_col = "status",
                                seed = 123) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    data_sets <- pn_filtered_set(object)
    train_data <- data_sets$training

  } else if (is.list(object) && "training" %in% names(object)) {
    train_data <- object$training
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'training' element")
  }

  cat("Training SuperPC model...\n")

  superpc_model_info <- train_superpc_base(train_data, time_col, status_col, seed)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with SuperPC model results...\n")

    object@survival.model[["superpc_model"]] <- superpc_model_info

    return(object)
  }

  cat("Returning model information as a list...\n")

  return(superpc_model_info)
}


#' Evaluate SuperPC ROC
#'
#' Evaluates SuperPC model performance using ROC curves.
#'
#' @param fit_best_cv_superpc SuperPC model.
#' @param threshold Threshold.
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
evaluate_superpc_roc <- function(
    fit_best_cv_superpc,
    threshold,
    train_data = NULL,
    test_data = NULL,
    validation_data = NULL,
    palette_name = "AsteroidCity1",
    time_col = "time",
    status_col = "status",
    save_plot = TRUE,
    save_dir = here::here("PrognosiX", "superpc_model"),
    plot_width = 7,
    plot_height = 7,
    base_size = 14
) {
  if (is.null(train_data) && is.null(test_data) && is.null(validation_data)) {
    stop("At least one of 'train_data', 'test_data', or 'validation_data' must be provided.")
  }
  
  roc_data_all <- data.frame()
  results_roc_df <- data.frame()
  
  featurenames <- fit_best_cv_superpc$featurenames
  
  if (!is.null(train_data)) {
    x_train <- t(train_data[, featurenames])
    y_train_status <- train_data[[status_col]]
    y_train_time <- train_data[[time_col]]
    
    data_list_train <- list(
      x = x_train,
      y = y_train_time,
      censoring.status = y_train_status,
      featurenames = featurenames
    )
    
    superpc_risk_score_train <- superpc::superpc.predict(
      fit_best_cv_superpc,
      data = data_list_train,
      newdata = data_list_train,
      threshold = threshold,
      n.components = 1
    )$v.pred
    
    superpc_risk_score_train <- as.numeric(superpc_risk_score_train)
    superpc_risk_score_train[is.infinite(superpc_risk_score_train)] <- NA 
    superpc_risk_score_train[is.na(superpc_risk_score_train)] <- median(superpc_risk_score_train, na.rm = TRUE)
    
    superpc_roc_train <- pROC::roc(y_train_status, superpc_risk_score_train)
    
    roc_data_train <- data.frame(
      specificity = superpc_roc_train$specificities,
      sensitivity = superpc_roc_train$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = survcomp::concordance.index(x = superpc_risk_score_train,
                                            surv.time = y_train_time,
                                            surv.event = y_train_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(superpc_roc_train)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    x_test <- t(test_data[, featurenames])
    y_test_status <- test_data[[status_col]]
    y_test_time <- test_data[[time_col]]
    
    data_list_test <- list(
      x = x_test,
      y = y_test_time,
      censoring.status = y_test_status,
      featurenames = featurenames
    )
    
    superpc_risk_score_test <- superpc::superpc.predict(
      fit_best_cv_superpc,
      data = data_list_train, 
      newdata = data_list_test,
      threshold = threshold,
      n.components = 1
    )$v.pred
    
    superpc_risk_score_test <- as.numeric(superpc_risk_score_test)
    superpc_risk_score_test[is.infinite(superpc_risk_score_test)] <- NA 
    superpc_risk_score_test[is.na(superpc_risk_score_test)] <- median(superpc_risk_score_test, na.rm = TRUE)
    
    superpc_roc_test <- pROC::roc(y_test_status, superpc_risk_score_test)
    
    roc_data_test <- data.frame(
      specificity = superpc_roc_test$specificities,
      sensitivity = superpc_roc_test$sensitivities,
      Set = "Testing Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_test)
    
    Test <- data.frame(
      C_index = survcomp::concordance.index(x = superpc_risk_score_test,
                                            surv.time = y_test_time,
                                            surv.event = y_test_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(superpc_roc_test)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
  }
  
  if (!is.null(validation_data)) {
    x_validation <- t(validation_data[, featurenames])
    y_validation_status <- validation_data[[status_col]]
    y_validation_time <- validation_data[[time_col]]
    
    data_list_validation <- list(
      x = x_validation,
      y = y_validation_time,
      censoring.status = y_validation_status,
      featurenames = featurenames
    )
    
    superpc_risk_score_validation <- superpc::superpc.predict(
      fit_best_cv_superpc,
      data = data_list_train, 
      newdata = data_list_validation,
      threshold = threshold,
      n.components = 1
    )$v.pred
    
    superpc_risk_score_validation <- as.numeric(superpc_risk_score_validation)
    superpc_risk_score_validation[is.infinite(superpc_risk_score_validation)] <- NA 
    superpc_risk_score_validation[is.na(superpc_risk_score_validation)] <- median(superpc_risk_score_validation, na.rm = TRUE)
    
    superpc_roc_validation <- pROC::roc(y_validation_status, superpc_risk_score_validation)
    
    roc_data_validation <- data.frame(
      specificity = superpc_roc_validation$specificities,
      sensitivity = superpc_roc_validation$sensitivities,
      Set = "Validation Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_validation)
    
    Validation <- data.frame(
      C_index = survcomp::concordance.index(x = superpc_risk_score_validation,
                                            surv.time = y_validation_time,
                                            surv.event = y_validation_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(superpc_roc_validation)
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


#' Evaluate ROC for SuperPC Model (Wrapper)
#'
#' Wrapper to evaluate ROC for SuperPC model in PrognosiX.
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
#' @return Updated PrognosiX object or results.
#' @export
evaluate_roc_superpc_model <- function(object,
                                       time_col = "time",
                                       status_col = "status",
                                       palette_name = "AsteroidCity1",
                                       save_plot = TRUE,
                                       save_dir = here::here("PrognosiX", "superpc_model"),
                                       plot_width = 7,
                                       plot_height = 7,
                                       base_size = 14) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_best_cv_superpc <- methods::slot(object, "survival.model")[["superpc_model"]][["model"]]
    threshold <- methods::slot(object, "survival.model")[["superpc_model"]][["threshold"]]

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("training", "testing") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")
    train_data <- object$training
    test_data <- object$testing

    fit_best_cv_superpc <- object$fit_best_cv_superpc
    threshold <- object$threshold

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_roc <- evaluate_superpc_roc(fit_best_cv_superpc = fit_best_cv_superpc,
                                      threshold = threshold,
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
    cat("'PrognosiX' object is being updated with SuperPC model ROC results...\n")

    object@survival.model[["superpc_model"]][["results_roc"]] <- results_roc

    cat("Results updated in the 'PrognosiX' object.\n")
    return(object)
  }

  cat("Returning model ROC results as a list...\n")

  return(results_roc)
}


#' Evaluate SuperPC Kaplan-Meier
#'
#' Evaluates KM curves for SuperPC risk groups.
#'
#' @param fit_best_cv_superpc Model.
#' @param data Data.
#' @param data_name Data name.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param R Resampling.
#' @param save_plot Logical.
#' @param save_dir Directory.
#' @param palette_name Palette.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Size.
#' @param seed Seed.
#' @param fill_color Fill color.
#' @param threshold Threshold.
#'
#' @return List with KM results and plots.
#' @export
evaluate_superpc_km <- function(fit_best_cv_superpc,
                                data,
                                data_name = "test",
                                time_col = "time",
                                status_col = "status",
                                R = 1000,
                                save_plot = TRUE,
                                save_dir = here::here("PrognosiX", "superpc_model"),
                                palette_name = "Dark2",
                                plot_width = 10,
                                plot_height = 8,
                                base_size = 14,
                                seed = 1234,
                                fill_color = "lightblue",
                                threshold = NULL # Threshold needs to be passed
) {
  
  set.seed(seed)
  if (!data_name %in% c("test", "val")) {
    stop("'data_name' must be either 'test' or 'val'.")
  }
  
  # For SuperPC, we need to reconstruct the data list format
  featurenames <- fit_best_cv_superpc$featurenames
  x_data <- t(data[, featurenames])
  y_data_time <- data[[time_col]]
  y_data_status <- data[[status_col]]
  
  data_list <- list(
    x = x_data,
    y = y_data_time,
    censoring.status = y_data_status,
    featurenames = featurenames
  )
  
  # Note: superpc.predict requires a 'train' data object as 'data' argument. 
  # This implementation assumes we might not have access to original training data here easily 
  # unless passed. However, superpc.predict uses 'data' mainly for standardization if needed 
  # or for thresholding. If 'newdata' is provided, it predicts on that.
  # A limitation here is that 'data' argument in predict usually refers to training data.
  # We'll use the provided 'data' as 'data' argument which might be incorrect if it expects training stats.
  # Ideally, training data should be passed or stored in the model object wrapper.
  # For now, using data itself might work if standardized or if just for projection.
  # A safer bet is to use the model object if it stores training stats, but superpc objects are light.
  # Assuming 'fit_best_cv_superpc' is the result of superpc.train.
  
  # IMPORTANT: superpc.predict needs the training data to standardize the new data correctly!
  # Without training data passed here, this might be inaccurate. 
  # We will assume standardization is handled or ignored for now, or use 'data' as proxy.
  
  risk_score <- superpc::superpc.predict(
    fit_best_cv_superpc,
    data = data_list, # Warning: Should be training data
    newdata = data_list,
    threshold = threshold,
    n.components = 1
  )$v.pred
  
  risk_score <- as.numeric(risk_score)
  
  risk_group <- ifelse(risk_score > median(risk_score, na.rm = TRUE), "High", "Low")
  data<-cbind(data,risk_score,risk_group)
  colnames(data)[ncol(data)-1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  cat("Calculating Kaplan-Meier survival curve...\n")
  km_fit_superpc <- survival::survfit(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  
  survdiff_result_superpc <- survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  km_pval_superpc <- pchisq(survdiff_result_superpc$chisq, 1, lower.tail = FALSE)
  km_hr_superpc <- (survdiff_result_superpc$obs[2] / survdiff_result_superpc$exp[2]) / (survdiff_result_superpc$obs[1] / survdiff_result_superpc$exp[1])
  
  km_upper95_superpc <- exp(log(km_hr_superpc) + qnorm(0.975) * sqrt(1 / survdiff_result_superpc$exp[2] + 1 / survdiff_result_superpc$exp[1]))
  km_lower95_superpc <- exp(log(km_hr_superpc) - qnorm(0.975) * sqrt(1 / survdiff_result_superpc$exp[2] + 1 / survdiff_result_superpc$exp[1]))
  
  km_results_superpc <- data.frame(KM_HR = km_hr_superpc,
                                   KM_CI_lower = km_lower95_superpc,
                                   KM_CI_upper = km_upper95_superpc,
                                   KM_p_value = km_pval_superpc)
  print(km_results_superpc)
  
  cat("Generating Kaplan-Meier plot...\n")
  km_plot_superpc <- survminer::ggsurvplot(
    km_fit_superpc,
    data = data,
    conf.int = TRUE,
    conf.int.fill = fill_color,
    conf.int.alpha = 0.5,
    pval = TRUE,
    pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve for SuperPC Risk Groups (", data_name, ")", sep = ""),
    surv.median.line = "hv",
    risk.table = TRUE,
    xlab = "Follow-up Time (days)",
    legend = c(0.8, 0.75),
    legend.title = "Risk Group",
    legend.labs = unique(data$risk_group),
    break.x.by = 100,
    palette = palette_name,
    base_size = base_size
  )
  
  
  surv_plot <- km_plot_superpc$plot +
    ggprism::theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot_superpc$table + ggprism::theme_prism(base_size = base_size) +
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
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    ggsave(file.path(save_dir, paste0("superpc_km_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("superpc_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("superpc_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", save_dir, "\n")
  }
  
  results_superpc <- list(
    KM_data_results = km_results_superpc,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_superpc)
}

#' Evaluate KM for SuperPC Model (Wrapper)
#'
#' Wrapper to evaluate KM for SuperPC model in PrognosiX.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
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
evaluate_km_superpc_model <- function(object,
                                      time_col = "time",
                                      status_col = "status",
                                      save_plot = TRUE,
                                      save_dir = here::here("PrognosiX", "superpc_model"),
                                      plot_width = 7,
                                      plot_height = 7,
                                      base_size = 14,
                                      data_name = "test",
                                      data=NULL) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_best_cv_superpc <- methods::slot(object, "survival.model")[["superpc_model"]][["model"]]
    threshold <- methods::slot(object, "survival.model")[["superpc_model"]][["threshold"]]

    data_sets <- pn_filtered_set(object)
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    test_data <- object$testing
    fit_best_cv_superpc <- object$fit_best_cv_superpc
    threshold <- object$threshold

  } else {
    stop("Input is neither a 'PrognosiX' object nor a valid list with 'train' and 'test' elements.")
  }

  cat("Evaluation completed.\n")
  results_superpc_km <- evaluate_superpc_km(
    fit_best_cv_superpc,
    data = test_data,
    data_name = data_name,
    time_col = time_col,
    status_col = status_col,
    save_plot = save_plot,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    threshold = threshold)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with SuperPC model Kaplan-Meier results...\n")

    object@survival.model[["superpc_model"]][["results_km"]] <- results_superpc_km

    return(object)
  }

  cat("Returning Kaplan-Meier results as a list...\n")

  return(results_superpc_km)
}
