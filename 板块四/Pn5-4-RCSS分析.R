library(survival)
library(rms)
library(ggplot2)

run_regression_analysis <- function(data, 
                                    method = c("linear", "logistic", "cox"),
                                    y_label, 
                                    time_col = NULL,  
                                    x_label, 
                                    knots = 3, 
                                    adjust_vars = NULL) { 
  method <- match.arg(method)
  
  if (!y_label %in% names(data)) stop("Outcome variable not found in data.")
  if (!x_label %in% names(data)) stop("Continuous variable not found in data.")
  if (method == "cox" && !time_col %in% names(data)) stop("Time variable not found for Cox regression.")
  
  if (!is.null(adjust_vars)) {
    missing_vars <- setdiff(adjust_vars, names(data))
    if (length(missing_vars) > 0) stop("Missing variables in the data: ", paste(missing_vars, collapse = ", "))
  }
  
  vars_to_check <- c(y_label, x_label)
  if (!is.null(adjust_vars)) {
    vars_to_check <- c(vars_to_check, adjust_vars)
  }
  if (!is.null(time_col)) {
    vars_to_check <- c(vars_to_check, time_col)
  }
  data <- na.omit(data[, vars_to_check])
  
  cat("Data dimensions after filtering: ", dim(data), "\n")
  
  cat("Data structure before modeling:\n")
  print(str(data))
  
  formula_str <- paste(y_label, "~ rcs(", x_label, ",", knots, ")")
  if (!is.null(adjust_vars)) {
    formula_str <- paste(formula_str, "+", paste(adjust_vars, collapse = "+"))
  }
  cat("Formula: ", formula_str, "\n")
  formula <- as.formula(formula_str)
  
  if (method == "cox") {
    cox_formula_str <- paste("Surv(", time_col, ",", y_label, ") ~ rcs(", x_label, ",", knots, ")")
    if (!is.null(adjust_vars)) {
      cox_formula_str <- paste(cox_formula_str, "+", paste(adjust_vars, collapse = "+"))
    }
    cat("Cox Formula: ", cox_formula_str, "\n")
    cox_formula <- as.formula(cox_formula_str)
  }
  
  if (method == "linear") {
    fit <- rms::ols(formula, data = data)
    
  } else if (method == "logistic") {
    fit <- rms::lrm(formula, data = data)
    
  } else if (method == "cox") {
    fit <- rms::cph(cox_formula, data = data, x = TRUE, y = TRUE)
  }
  
  cat("Model Summary:\n")
  print(summary(fit))
  
  return(fit)
}

library(rcssci)
plot_regression_analysis <- function(data, 
                                     method = c("linear", "logistic", "cox"),
                                     y_label, 
                                     time_col = NULL,  
                                     x_label, 
                                     knots = 3,  
                                     adjust_vars = NULL,  
                                     prob = 0.1,  
                                     save_dir = here('PrognosiX', "univariate_analysis")) {  
  

    filepath <-file.path(save_dir) 
   
  
  if (!y_label %in% names(data)) stop("Outcome variable not found in data.")
  if (!x_label %in% names(data)) stop("Continuous variable not found in data.")
  if (method == "cox" && is.null(time_col)) stop("Time variable is required for Cox regression.")
  if (method == "cox" && !time_col %in% names(data)) stop("Time variable not found in data.")
  
  vars_to_check <- c(y_label, x_label, adjust_vars, time_col)
  data <- na.omit(data[, vars_to_check])
  data[[y_label]] <- as.factor(data[[y_label]])
  if (method == "cox") {
    result <- rcssci::rcssci_cox(data = data, 
                             y = y_label, 
                             x = x_label, 
                             covs = adjust_vars, 
                             time = time_col, 
                             prob = prob, 
                             filepath = filepath)
    
  } else if (method == "logistic") {
    result <- rcssci::rcssci_logistic(data = data, 
                              y = y_label, 
                              x = x_label, 
                              prob = prob, 
                              filepath = filepath)
    
  } else if (method == "linear") {
    result <- rcssci::rcssci_linear(data = data, 
                            y = y_label, 
                            x = x_label, 
                            prob = prob, 
                            filepath = filepath)
  }
  
  return(result)
}

Pron_univariate_regression <- function(object,  
                                       method = "cox",
                                       y_label, 
                                       time_col, 
                                       x_label, 
                                       knots = 3,  
                                       adjust_vars = NULL,  
                                       prob = 0.1,  
                                       save_dir = here('PrognosiX', "univariate_analysis")) {
  
  # Check if the object is a 'PrognosiX' object or a data frame
  if (inherits(object, 'PrognosiX')) {
    data <- slot(object, "survival.data")
    y_label <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
    
    if (is.null(data) || nrow(data) == 0) {
      stop("The survival.data in the PrognosiX object is empty.")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }
  
  # Call the plot_regression_analysis function
  plot_fit <- plot_regression_analysis(data = data, 
                                       method = method,
                                       y_label = y_label, 
                                       time_col = time_col,  
                                       x_label = x_label, 
                                       adjust_vars = adjust_vars,  
                                       prob = prob)
  
  # Get knots from plot_fit
  knots <- plot_fit[["kn"]]
  
  # Run the final regression
  final_fit <- run_regression_analysis(data = data, 
                                       method = method,
                                       y_label = y_label, 
                                       time_col = time_col,  
                                       x_label = x_label, 
                                       knots = knots,
                                       adjust_vars = adjust_vars)
  
  # Update the PrognosiX object with the results
  object@univariate.analysis[["rcs_analysis"]] <- list(plot = plot_fit, fit = final_fit)
  
  cat("The 'PrognosiX' object has been updated with the following slots:\n")
  cat("- 'univariate.analysis' slot updated.\n")
  
  return(object)
}

