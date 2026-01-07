run_univariate_cox_analysis <- function(data,
                                        time_col = "time",
                                        status_col = "status",
                                        selected_vars = NULL,
                                        covariates = NULL,
                                        P_value =0.05) {
  cat("Starting univariate Cox regression analysis for selected variables...\n")
  
  if (is.null(selected_vars)) {
    selected_vars <- setdiff(names(data), c(status_col, time_col, covariates))
  }
  
  if (length(selected_vars) == 0) {
    cat("No variables specified for analysis. Skipping Cox analysis.\n")
    return(NULL)
  }
  
  results <- data.frame(
    Variable = character(0),
    HR = numeric(0),
    CI_lower = numeric(0),
    CI_upper = numeric(0),
    P_value = numeric(0),
    HR_95CI = character(0),
    se = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (var_col in selected_vars) {
    cat("Analyzing variable:", var_col, "\n")
    
    if (!(var_col %in% names(data))) {
      cat("Variable", var_col, "does not exist in the data frame. Skipping...\n")
      next
    }
    
    # Create formula with covariates if specified
    if (!is.null(covariates) && length(covariates) > 0) {
      formula_str <- paste("Surv(", time_col, ",", status_col, ") ~", var_col, "+", paste(covariates, collapse = " + "))
      cat("Adjusting for covariates:", paste(covariates, collapse = ", "), "\n")
    } else {
      formula_str <- paste("Surv(", time_col, ",", status_col, ") ~", var_col)
    }
    
    temp_data <- data %>% filter(!is.na(data[[var_col]]) & !is.na(data[[status_col]]))
    
    # Also check for NA in covariates if they exist
    if (!is.null(covariates) && length(covariates) > 0) {
      temp_data <- temp_data %>% filter(complete.cases(temp_data[, covariates]))
    }
    
    if (nrow(temp_data) == 0) {
      cat("Filtered data for", var_col, "has no valid rows. Skipping...\n")
      next
    }
    
    cox_model <- coxph(as.formula(formula_str), data = temp_data)
    cox_summary <- summary(cox_model)
    
    # Extract results for the variable of interest (not covariates)
    # Find the row corresponding to our variable of interest
    var_row <- which(rownames(cox_summary$coefficients) == var_col)
    
    if (length(var_row) == 0) {
      cat("Could not find results for variable", var_col, "in the model output. Skipping...\n")
      next
    }
    
    hr <- cox_summary$coefficients[var_row, "exp(coef)"]
    ci_lower <- cox_summary$conf.int[var_row, "lower .95"]
    ci_upper <- cox_summary$conf.int[var_row, "upper .95"]
    p_value <- cox_summary$coefficients[var_row, "Pr(>|z|)"]
    
    hr_95ci <- paste0(
      round(hr, 2),
      " (",
      round(ci_lower, 2),
      "-",
      round(ci_upper, 2),
      ")"
    )
    
    se <- (log(ci_upper) - log(hr)) / 1.96
    
    results <- rbind(results, data.frame(
      Variable = var_col,
      HR = hr,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      P_value = p_value,
      HR_95CI = hr_95ci,
      se = se,
      stringsAsFactors = FALSE
    ))
  }
  
  # Subset results where p-value < 0.05 (significant results)
  significant_results <- results[results$P_value <P_value, ]
  
  cat("Univariate Cox regression analysis for selected variables completed.\n")
  return(list(all_results = results, significant_results = significant_results))
}




Sur_cox_univariate_analysis <- function(object,
                                        formula = NULL,
                                        status_col = "status",
                                        time_col = "time",
                                        selected_vars = NULL,
                                        covariates = NULL,
                                        P_value =0.05,
                                        save_dir = here('SurObj', "univariate_analysis")) {
  
  if (inherits(object, 'SurObj')) {
    example_data <- slot(object, "survival.data")
    cat("Using survival data from SurObj object...\n")
  } else if (is.data.frame(object)) {
    example_data <- object
  } else {
    stop("Input must be an object of class 'SurObj' or a data frame.")
  }
  
  if (is.null(example_data) || nrow(example_data) == 0)
    stop("No valid data found in the input.")
  
  
  
  example_data[[status_col]] <- as.numeric(example_data[[status_col]])
  
  results <- run_univariate_cox_analysis(
    data = example_data,
    time_col = time_col,
    status_col = status_col,
    selected_vars = selected_vars,
    covariates =covariates,
    P_value=P_value
  )
  
  if (inherits(object, 'SurObj')) {
    cat("Updating 'SurObj' object...\n")
    
    object@univariate.analysis[["all_univariate_results"]] <- results
    
    cat("The 'SurObj' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }
  
  return(results)
}
