#' Train All Models
#'
#' Trains multiple survival models specified in the list.
#'
#' @param object PrognosiX object.
#' @param model_list List of models to train ("ridge", "lasso", "pls", "coxboost", "coxph", "superpc", "rsf").
#'
#' @return Updated PrognosiX object.
#' @export
train_all_models <- function(object,
                             model_list = c("ridge", "lasso", "pls", "coxboost", "coxph", "superpc", "rsf")) {
  
  if ("ridge" %in% model_list) {
    tryCatch({
      library(survcomp)
      cat("Training Ridge model...\n")
      object <- train_ridge_model(object)
      object <- evaluate_roc_ridge_model(object)
      object <- evaluate_km_ridge_model(object)
      object <- ridge_compute_hr_and_ci(object)
      object<- forest_plot_ridge_model(object)
      cat("Ridge model training and evaluation completed.\n")
    }, error = function(e) {
      cat("Error in Ridge model training: ", e$message, "\n")
    })
  }
  
  if ("lasso" %in% model_list) {
    tryCatch({
      library(survcomp)
      cat("Training Lasso model...\n")
      object <- train_lasso_model(object)
      object <- evaluate_roc_lasso_model(object)
      object <- evaluate_km_lasso_model(object)
      object <- lasso_compute_hr_and_ci(object)
      object <- forest_plot_lasso_model(object)
      cat("Lasso model training and evaluation completed.\n")
    }, error = function(e) {
      cat("Error in Lasso model training: ", e$message, "\n")
    })
  }
  
  if ("pls" %in% model_list) {
    tryCatch({
      library(survcomp)
      cat("Training PLS model...\n")
      object <- train_pls_model(object)
      object <- evaluate_roc_pls_model(object)
      object <- evaluate_km_pls_model(object)
      object <- pls_compute_hr_and_ci(object)
      object <- forest_plot_pls_model(object)
      cat("PLS model training and evaluation completed.\n")
    }, error = function(e) {
      cat("Error in PLS model training: ", e$message, "\n")
    })
  }
  
  if ("coxboost" %in% model_list) {
    tryCatch({
      library(survcomp)
      cat("Training CoxBoost model...\n")
      object <- train_coxboost_model(object)
      object <- evaluate_roc_coxboost_model(object)
      object <- evaluate_km_coxboost_model(object)
      object <- coxboost_coefs_compute_hr_and_ci(object)
      object <- forest_plot_coxboost_model(object)
      cat("CoxBoost model training and evaluation completed.\n")
    }, error = function(e) {
      cat("Error in CoxBoost model training: ", e$message, "\n")
    })
  }
  
  if ("coxph" %in% model_list) {
    tryCatch({
      library(survcomp)
      cat("Training CoxPH model...\n")
      object <- train_coxph_model(object)
      object <- evaluate_roc_coxph_model(object)
      object <- evaluate_km_coxph_model(object)
      object <- coxph_compute_hr_and_ci(object)
      object <- forest_plot_coxph_model(object)
      cat("CoxPH model training and evaluation completed.\n")
    }, error = function(e) {
      cat("Error in CoxPH model training: ", e$message, "\n")
    })
  }
  
  if ("superpc" %in% model_list) {
    tryCatch({
      library(survcomp)
      cat("Training SuperPC model...\n")
      object <- train_superpc_model(object)
      object <- evaluate_roc_superpc_model(object)
      object <- evaluate_km_superpc_model(object)
      # SuperPC HR plotting not standard, skipping forest plot or reusing generic
      cat("SuperPC model training and evaluation completed.\n")
    }, error = function(e) {
      cat("Error in SuperPC model training: ", e$message, "\n")
    })
  }
  
  if ("rsf" %in% model_list) {
    tryCatch({
      library(survcomp)
      cat("Training RSF model...\n")
      object <- train_rsf_model(object)
      object <- evaluate_roc_rsf_model(object)
      object <- evaluate_km_rsf_model(object)
      object <- rsf_compute_hr_and_ci(object)
      object <- forest_plot_rsf_model(object)
      cat("RSF model training and evaluation completed.\n")
    }, error = function(e) {
      cat("Error in RSF model training: ", e$message, "\n")
    })
  }
  
  cat("Model training and evaluation completed.\n")
  return(object)
}

#' Combine Model ROC Results
#'
#' Combines ROC results from all trained models in a PrognosiX object.
#'
#' @param object PrognosiX object.
#'
#' @return Data frame of combined results.
#' @export
combine_model_roc_results <- function(object) {
  models <- setdiff(names(object@survival.model), "all_results")
  results_df <- data.frame()
  
  for (model_name in models) {
    roc <- object@survival.model[[model_name]][["results_roc"]][["results_roc"]]
    
    if (is.null(roc) || nrow(roc) == 0) {
      warning(paste("No ROC results found for model:", model_name))
      next
    }
    
    roc$Model <- model_name
    results_df <- rbind(results_df, roc)
  }
  
  return(results_df)
}

#' Extract Model Results
#'
#' Extracts and identifies the best model based on a metric.
#'
#' @param object PrognosiX object.
#' @param metric Metric to optimize ("C_index" or "ROC_AUC").
#'
#' @return Updated PrognosiX object.
#' @export
extract_model_results <- function(object, 
                                  metric = "C_index") {
  
  results_df <- combine_model_roc_results(object)
  
  if (nrow(results_df) == 0) {
    warning("No models have valid ROC results.")
    return(object)
  }
  
  if (!metric %in% c("C_index", "ROC_AUC")) {
    stop("Invalid metric. Choose either 'C_index' or 'ROC_AUC'.")
  }
  
  results_no_na <- na.omit(results_df)
  
  best_model <- NULL
  
  if (nrow(results_no_na) > 0) {
    best_model_train <- results_no_na[results_no_na$Dataset == "Train", ]
    best_model_test <- results_no_na[results_no_na$Dataset == "Test", ]
    
    if (nrow(best_model_train) > 0 && nrow(best_model_test) > 0) {
      best_model_train <- best_model_train[order(-best_model_train[[metric]]), ]
      best_model_train <- best_model_train[1, ]
      
      best_model_test <- best_model_test[best_model_test$Model == best_model_train$Model, ]
      
      if (nrow(best_model_test) > 0) {
        best_model <- rbind(best_model_train, best_model_test)
      } else {
        best_model <- best_model_train
      }
    }
  }
  
  if (is.null(best_model) || nrow(best_model) == 0) {
    warning("No valid best model found.")
    return(object)
  }
  
  best_model_name <- unique(best_model$Model)
  
  if (inherits(object, 'PrognosiX')) {  
    object@survival.model[["all_results"]] <- results_df
    model <- object@survival.model[[best_model_name]]
    object@best.model[["model"]] <- model
    object@best.model[["model_name"]] <- best_model_name
    object@best.model[["best_model_results"]] <- best_model
    
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'best.model' slot updated.\n")
    cat("- 'survival.model' slot updated.\n")
  }
  
  return(object)
}
