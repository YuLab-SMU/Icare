#' Run Univariate Cox Regression Analysis
#'
#' This function performs univariate Cox regression analysis for a list of selected variables.
#'
#' @param data A data frame containing the survival data.
#' @param time_col The column name for time-to-event data.
#' @param status_col The column name for event status (0/1).
#' @param selected_vars A character vector of variable names to analyze. If NULL, all variables except time and status are used.
#'
#' @return A list containing a data frame of all results and a data frame of significant results (p < 0.05).
#' @export
run_univariate_cox_analysis <- function(data,
                                        time_col = "time",
                                        status_col = "status",
                                        selected_vars = NULL) {
  cat("Starting univariate Cox regression analysis for selected variables...\n")

  if (is.null(selected_vars)) {
    selected_vars <- setdiff(names(data), c(status_col, time_col))
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

    temp_data <- data %>% dplyr::filter(!is.na(data[[var_col]]) & !is.na(data[[status_col]]))
    if (nrow(temp_data) == 0) {
      cat("Filtered data for", var_col, "has no valid rows. Skipping...\n")
      next
    }

    cox_model <- survival::coxph(as.formula(paste("Surv(", time_col, ",", status_col, ") ~", var_col)), data = temp_data)
    cox_summary <- summary(cox_model)

    hr <- cox_summary$coefficients[1, "exp(coef)"]
    ci_lower <- cox_summary$conf.int[1, "lower .95"]
    ci_upper <- cox_summary$conf.int[1, "upper .95"]
    p_value <- cox_summary$coefficients[1, "Pr(>|z|)"]

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
  significant_results <- results[results$P_value < 0.05, ]

  cat("Univariate Cox regression analysis for selected variables completed.\n")
  return(list(all_results = results, significant_results = significant_results))
}

#' PrognosiX Univariate Analysis Wrapper
#'
#' Wrapper function to run univariate Cox analysis on a PrognosiX object or data frame.
#'
#' @param object A PrognosiX object or a data frame.
#' @param formula Optional formula (not currently used in logic).
#' @param status_col Status column name.
#' @param time_col Time column name.
#' @param selected_vars Variables to analyze.
#' @param response_var Optional response variable override.
#' @param save_plots Logical, whether to save plots (not used in this function but kept for consistency).
#' @param save_dir Directory to save results.
#' @param use_subgroup_data Logical, whether to use subgroup data from PrognosiX object.
#'
#' @return Updated PrognosiX object or results data frame.
#' @export
Prognos_cox_univariate_analysis <- function(object,
                                            formula = NULL,
                                            status_col = "status",
                                            time_col = "time",
                                            selected_vars = NULL,
                                            response_var = NULL,
                                            save_plots = TRUE,
                                            save_dir = here::here('PrognosiX', "univariate_analysis"),
                                            use_subgroup_data = FALSE) {

  if (inherits(object, 'PrognosiX')) {
    if (use_subgroup_data) {
      example_data <- methods::slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
    } else {
      example_data <- methods::slot(object, "survival.data")
      cat("Using original survival data...\n")
    }
    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")
  } else if (is.data.frame(object)) {
    example_data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  if (is.null(example_data) || nrow(example_data) == 0)
    stop("No valid data found in the input.")

  if (is.null(response_var)) {
    response_var <- status_col
  }
  example_data[[status_col]]<-as.numeric(example_data[[status_col]])
  results <- run_univariate_cox_analysis(
    data = example_data,
    time_col = time_col,
    status_col = response_var,
    selected_vars = selected_vars
  )

  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object...\n")

    object@univariate.analysis[["all_univariate_results"]] <- results


    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }

  return(results)
}
