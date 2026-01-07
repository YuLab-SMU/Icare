#' Train RSF Model (Base)
#'
#' Trains a Random Survival Forest (RSF) model (base function).
#'
#' @param train_data Training data.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param ntree Number of trees.
#' @param nodesize Node size.
#' @param seed Seed.
#'
#' @return List containing model.
#' @export
train_rsf_base <- function(train_data,
                           time_col = "time",
                           status_col = "status",
                           ntree = 1000,
                           nodesize = 5,
                           seed = 123) {

  set.seed(seed)
  formula <- as.formula(paste("Surv(", time_col, ", ", status_col, ") ~ ."))
  fit_rsf <- randomForestSRC::rfsrc(formula, data = train_data, ntree = ntree, nodesize = nodesize)

  return(list(model = fit_rsf))
}


#' Train RSF Model for PrognosiX
#'
#' Wrapper to train RSF model on PrognosiX object.
#'
#' @param object PrognosiX object.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param ntree Number of trees.
#' @param nodesize Node size.
#' @param seed Seed.
#'
#' @return Updated PrognosiX object or model list.
#' @export
train_rsf_model <- function(object,
                            time_col = "time",
                            status_col = "status",
                            ntree = 1000,
                            nodesize = 5,
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

  cat("Training RSF model...\n")

  rsf_model_info <- train_rsf_base(train_data, time_col, status_col, ntree, nodesize, seed)

  if (inherits(object, 'PrognosiX')) {
    cat("'PrognosiX' object is being updated with RSF model results...\n")

    object@survival.model[["rsf_model"]] <- rsf_model_info

    return(object)
  }

  cat("Returning model information as a list...\n")

  return(rsf_model_info)
}


#' Evaluate RSF ROC
#'
#' Evaluates RSF model performance using ROC curves.
#'
#' @param fit_rsf_model RSF model.
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
evaluate_rsf_roc <- function(
    fit_rsf_model,
    train_data = NULL,
    test_data = NULL,
    validation_data = NULL,
    palette_name = "AsteroidCity1",
    time_col = "time",
    status_col = "status",
    save_plot = TRUE,
    save_dir = here::here("PrognosiX", "rsf_model"),
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
    y_train_status <- train_data[[status_col]]
    y_train_time <- train_data[[time_col]]
    
    rsf_risk_score_train <- predict(fit_rsf_model, newdata = train_data)$predicted
    rsf_risk_score_train[is.infinite(rsf_risk_score_train)] <- NA 
    rsf_risk_score_train[is.na(rsf_risk_score_train)] <- median(rsf_risk_score_train, na.rm = TRUE)
    
    rsf_roc_train <- pROC::roc(y_train_status, as.numeric(rsf_risk_score_train))
    
    roc_data_train <- data.frame(
      specificity = rsf_roc_train$specificities,
      sensitivity = rsf_roc_train$sensitivities,
      Set = "Training Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_train)
    
    Train <- data.frame(
      C_index = survcomp::concordance.index(x = rsf_risk_score_train,
                                            surv.time = y_train_time,
                                            surv.event = y_train_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(rsf_roc_train)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Train", Train))
  }
  
  if (!is.null(test_data)) {
    y_test_status <- test_data[[status_col]]
    y_test_time <- test_data[[time_col]]
    
    rsf_risk_score_test <- predict(fit_rsf_model, newdata = test_data)$predicted
    rsf_risk_score_test[is.infinite(rsf_risk_score_test)] <- NA 
    rsf_risk_score_test[is.na(rsf_risk_score_test)] <- median(rsf_risk_score_test, na.rm = TRUE)
    
    rsf_roc_test <- pROC::roc(y_test_status, as.numeric(rsf_risk_score_test))
    
    roc_data_test <- data.frame(
      specificity = rsf_roc_test$specificities,
      sensitivity = rsf_roc_test$sensitivities,
      Set = "Testing Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_test)
    
    Test <- data.frame(
      C_index = survcomp::concordance.index(x = rsf_risk_score_test,
                                            surv.time = y_test_time,
                                            surv.event = y_test_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(rsf_roc_test)
    )
    results_roc_df <- rbind(results_roc_df, data.frame(Dataset = "Test", Test))
  }
  
  if (!is.null(validation_data)) {
    y_validation_status <- validation_data[[status_col]]
    y_validation_time <- validation_data[[time_col]]
    
    rsf_risk_score_validation <- predict(fit_rsf_model, newdata = validation_data)$predicted
    rsf_risk_score_validation[is.infinite(rsf_risk_score_validation)] <- NA 
    rsf_risk_score_validation[is.na(rsf_risk_score_validation)] <- median(rsf_risk_score_validation, na.rm = TRUE)
    
    rsf_roc_validation <- pROC::roc(y_validation_status, as.numeric(rsf_risk_score_validation))
    
    roc_data_validation <- data.frame(
      specificity = rsf_roc_validation$specificities,
      sensitivity = rsf_roc_validation$sensitivities,
      Set = "Validation Set"
    )
    roc_data_all <- rbind(roc_data_all, roc_data_validation)
    
    Validation <- data.frame(
      C_index = survcomp::concordance.index(x = rsf_risk_score_validation,
                                            surv.time = y_validation_time,
                                            surv.event = y_validation_status,
                                            method = "noether")$c.index,
      ROC_AUC = pROC::auc(rsf_roc_validation)
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


#' Evaluate ROC for RSF Model (Wrapper)
#'
#' Wrapper to evaluate ROC for RSF model in PrognosiX.
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
evaluate_roc_rsf_model <- function(object,
                                   time_col = "time",
                                   status_col = "status",
                                   palette_name = "AsteroidCity1",
                                   save_plot = TRUE,
                                   save_dir = here::here("PrognosiX", "rsf_model"),
                                   plot_width = 7,
                                   plot_height = 7,
                                   base_size = 14) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_rsf_model <- methods::slot(object, "survival.model")[["rsf_model"]][["model"]]

    data_sets <- pn_filtered_set(object)

    train_data <- data_sets$training
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("training", "testing") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")
    train_data <- object$training
    test_data <- object$testing

    fit_rsf_model <- object$fit_rsf_model

  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'train' and 'test' elements")
  }

  results_roc <- evaluate_rsf_roc(fit_rsf_model = fit_rsf_model,
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
    cat("'PrognosiX' object is being updated with RSF model ROC results...\n")

    object@survival.model[["rsf_model"]][["results_roc"]] <- results_roc

    cat("Results updated in the 'PrognosiX' object.\n")
    return(object)
  }

  cat("Returning model ROC results as a list...\n")

  return(results_roc)
}


#' Evaluate RSF Kaplan-Meier
#'
#' Evaluates KM curves for RSF risk groups.
#'
#' @param fit_best_rsf Model.
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
#'
#' @return List with KM results and plots.
#' @export
evaluate_rsf_km <- function(fit_best_rsf,
                            data,
                            data_name = "test",
                            time_col = "time",
                            status_col = "status",
                            R = 1000,
                            save_plot = TRUE,
                            save_dir = here::here("PrognosiX", "rsf_model"),
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
  
  risk_score <- predict(fit_best_rsf, newdata = data)$predicted
  
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  data<-cbind(data,risk_score,risk_group)
  colnames(data)[ncol(data)-1] <- "risk_score"
  colnames(data)[ncol(data)] <- "risk_group"
  
  cat("Calculating Kaplan-Meier survival curve...\n")
  km_fit_rsf <- survival::survfit(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  
  survdiff_result_rsf <- survival::survdiff(survival::Surv(data[[time_col]], as.numeric(data[[status_col]])) ~ risk_group, data = data)
  km_pval_rsf <- pchisq(survdiff_result_rsf$chisq, 1, lower.tail = FALSE)
  km_hr_rsf <- (survdiff_result_rsf$obs[2] / survdiff_result_rsf$exp[2]) / (survdiff_result_rsf$obs[1] / survdiff_result_rsf$exp[1])
  
  km_upper95_rsf <- exp(log(km_hr_rsf) + qnorm(0.975) * sqrt(1 / survdiff_result_rsf$exp[2] + 1 / survdiff_result_rsf$exp[1]))
  km_lower95_rsf <- exp(log(km_hr_rsf) - qnorm(0.975) * sqrt(1 / survdiff_result_rsf$exp[2] + 1 / survdiff_result_rsf$exp[1]))
  
  km_results_rsf <- data.frame(KM_HR = km_hr_rsf,
                               KM_CI_lower = km_lower95_rsf,
                               KM_CI_upper = km_upper95_rsf,
                               KM_p_value = km_pval_rsf)
  print(km_results_rsf)
  
  cat("Generating Kaplan-Meier plot...\n")
  km_plot_rsf <- survminer::ggsurvplot(
    km_fit_rsf,
    data = data,
    conf.int = TRUE,
    conf.int.fill = fill_color,
    conf.int.alpha = 0.5,
    pval = TRUE,
    pval.method = TRUE,
    title = paste("Kaplan-Meier Survival Curve for RSF Risk Groups (", data_name, ")", sep = ""),
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
  
  
  surv_plot <- km_plot_rsf$plot +
    ggprism::theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )
  
  risk_table <- km_plot_rsf$table + ggprism::theme_prism(base_size = base_size) +
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
    ggsave(file.path(save_dir, paste0("rsf_km_", data_name, ".pdf")), plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("rsf_curve_", data_name, ".pdf")), plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("rsf_risk_table_", data_name, ".pdf")), plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", save_dir, "\n")
  }
  
  results_rsf <- list(
    KM_data_results = km_results_rsf,
    combined_plot = combined_plot,
    data_risk = data
  )
  
  return(results_rsf)
}

#' Evaluate KM for RSF Model (Wrapper)
#'
#' Wrapper to evaluate KM for RSF model in PrognosiX.
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
evaluate_km_rsf_model <- function(object,
                                  time_col = "time",
                                  status_col = "status",
                                  save_plot = TRUE,
                                  save_dir = here::here("PrognosiX", "rsf_model"),
                                  plot_width = 7,
                                  plot_height = 7,
                                  base_size = 14,
                                  data_name = "test",
                                  data=NULL) {

  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
    fit_rsf_model <- methods::slot(object, "survival.model")[["rsf_model"]][["model"]]

    data_sets <- pn_filtered_set(object)
    test_data <- data_sets$testing

  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")

    test_data <- object$testing
    fit_rsf_model <- object$fit_rsf_model

  } else {
    stop("Input is neither a 'PrognosiX' object nor a valid list with 'train' and 'test' elements.")
  }

  cat("Evaluation completed.\n")
  results_rsf_km <- evaluate_rsf_km(
    fit_rsf_model,
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
    cat("'PrognosiX' object is being updated with RSF model Kaplan-Meier results...\n")

    object@survival.model[["rsf_model"]][["results_km"]] <- results_rsf_km

    return(object)
  }

  cat("Returning Kaplan-Meier results as a list...\n")

  return(results_rsf_km)
}


#' Compute HR and CI for RSF Coefficients
#'
#' Computes HR and CI for RSF coefficients (variable importance).
#'
#' @param object PrognosiX object or model.
#'
#' @return Updated PrognosiX object or results.
#' @export
rsf_compute_hr_and_ci <- function(object) {
  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting data...\n")

    fit_rsf <- methods::slot(object, "survival.model")[["rsf_model"]][["model"]]
  } else if (is.list(object) && "model" %in% names(object)) {
    fit_rsf <- object$model
  } else {
    stop("Input must be an object of class 'PrognosiX' or a list with 'model'.")
  }

  importance <- fit_rsf$importance
  feature_names <- names(importance)

  # RSF does not provide HR directly. Using Importance as proxy or placeholder.
  # A rigorous HR calculation from RSF is non-trivial.
  # We will output importance values.
  
  hr_results <- data.frame(
    Variable = feature_names,
    Importance = importance
  )
  
  # Dummy columns to match structure if needed, or just return importance
  hr_results$HR <- NA
  hr_results$CI_lower <- NA
  hr_results$CI_upper <- NA
  hr_results$HR_95CI <- NA


  print(hr_results)

  if (inherits(object, 'PrognosiX')) {
    object@survival.model[["rsf_model"]][["hr_results"]] <- hr_results
    cat("Updating 'PrognosiX' object...\n")
    return(object)
  }

  return(hr_results)
}


#' Create Forest Plot for RSF
#'
#' Internal function to create forest plot for RSF (Importance plot).
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
#' @param ci_range_limit CI limit.
#'
#' @return Plot.
#' @export
create_forest_plot_rsf <- function(hr_results,
                                     plot_title = "Variable Importance for RSF Model",
                                     save_plot = FALSE,
                                     save_dir = here::here('Prognosi'),
                                     plot_width = 11,
                                     plot_height = 5,
                                     palette_name = "AsteroidCity1",
                                     base_size = 14,
                                     hr_limit = c(0, 3),
                                     ci_range_limit = 1000) {

  # For RSF, we plot Importance instead of HR Forest Plot
  
  hr_results <- hr_results %>%
    dplyr::arrange(desc(Importance))
  
  p <- ggplot2::ggplot(hr_results, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = wesanderson::wes_palette(palette_name, 1)) +
    coord_flip() +
    labs(title = plot_title, x = "Variable", y = "Importance") +
    theme_classic(base_size = base_size)

  print(p)

  if (save_plot) {
    pdf_path <- file.path(save_dir, "forest_plot.pdf")
    ggsave(filename = pdf_path, plot = p, width = plot_width, height = plot_height)
    cat("Plot saved at: ", pdf_path, "\n")
  }
  return(p)
}


#' Forest Plot for RSF Model (Wrapper)
#'
#' Generates forest plot (importance) for RSF model in PrognosiX.
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
forest_plot_rsf_model <- function(object,
                                  plot_title = "Variable Importance for RSF Model",
                                  time_col = "time",
                                  status_col = "status",
                                  var_col = NULL,
                                  palette_name = "AsteroidCity1",
                                  save_plot = TRUE,
                                  save_dir = here::here('PrognosiX', "rsf_model"),
                                  plot_width = 14,
                                  plot_height = 7,
                                  base_size = 14,
                                  use_subgroup_data = FALSE,
                                  hr_limit = c(0.1, 3)
) {

  if (inherits(object, 'PrognosiX')) {
    hr_results <- methods::slot(object, "survival.model")[["rsf_model"]][["hr_results"]]

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

  p <- create_forest_plot_rsf(hr_results,
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
