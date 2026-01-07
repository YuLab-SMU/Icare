#' Evaluate and Save Model Performance Metrics
#'
#' This function evaluates model performance across training, testing, validation, 
#' and external validation datasets, and saves the results to CSV files.
#'
#' @param object A Best_Model object containing the model and datasets to evaluate
#' @param group_col Name of the response variable column (default: "group")
#' @param best_threshold Custom classification threshold (default: NULL uses model's optimal threshold)
#' @param save_data Logical indicating whether to save results to CSV (default: TRUE)
#' @param save_dir Directory path to save CSV files (default: here::here("ModelData", "performance_results"))
#'
#' @return The updated Best_Model object with performance results stored in @performance.result slot.
#'         When save_data=TRUE, also saves:
#'         - performance_summary.csv: Combined metrics across all datasets
#'         - individual CSV files for each dataset (training_metrics.csv, etc.)
#'
#' @examples
#' \dontrun{
#' # With custom parameters
#' object_best <- ModelPerformance(
#'   object = object_best,
#'   best_threshold = 0.4,
#'   save_dir = "model_performance"
#' )
#' }
#'
#' @export
ModelPerformance <- function(object,
                             group_col = "group",
                             best_threshold = NULL,
                             save_dir = here::here("ModelData", "best_model_result"),
                             save_data = TRUE,
                             csv_filename = "performance_summary.csv") {
  
  cat("\n=== Starting Model Performance Evaluation ===\n")
  
  if (!inherits(object, "Best_Model")) {
    stop("Input must be an object of class 'Best_Model'")
  }
  cat("Input validation passed - Best_Model object detected\n")
  
  best_model <- object@best.model[[1]]
  best_model_type <- object@best.model.type
  group_col <- object@group_col
  data_sets <- object@filtered.set
  
  if (is.null(best_model)) {
    stop("Best model not found in the object")
  }
  cat("Best model successfully extracted from object\n")
  
  if (is.null(best_threshold)) {
    best_threshold <- object@best_threshold[["best_threshold"]]
    cat("Using model's optimal threshold:", best_threshold, "\n")
  } else {
    cat("Using custom threshold:", best_threshold, "\n")
  }
  
  eval_results <- list()
  
  cat("\nBeginning evaluation across datasets...\n")
  
  for (set_name in c("training", "testing", "validation", "external_validation")) {
    data <- data_sets[[set_name]]
    
    if (!is.null(data) && nrow(data) > 0) {
      cat("\n[Evaluating] Dataset:", set_name, 
          "| Rows:", nrow(data), 
          "| Columns:", ncol(data), "\n")
      
      res <- evaluate_model_performance(
        data = data,
        model_result = best_model,
        group_col = group_col,
        custom_cutoff = best_threshold
      )
      
      eval_results[[set_name]] <- res
      cat("[Completed] Evaluation for", set_name, "dataset\n")
      
    } else {
      cat("\n[Skipped] Dataset:", set_name, "is NULL or empty\n")
    }
  }
  
  cat("\nCombining evaluation results...\n")
  summary_df <- do.call(rbind, lapply(names(eval_results), function(nm) {
    df <- eval_results[[nm]]
    df$Dataset <- nm
    df
  }))
  
  summary_df <- summary_df[, c("Dataset", setdiff(names(summary_df), "Dataset"))]
  
  cat("\n=== Performance Summary ===\n")
  print(summary_df)
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("Created directory:", save_dir, "\n")
    }
    
    summary_path <- file.path(save_dir, csv_filename)
    write.csv(summary_df, summary_path, row.names = FALSE)
    cat("Saved combined performance metrics to:", summary_path, "\n")
    
  }
  
  object@performance.result[["summary_performance"]] <- summary_df
  cat("\nUpdating 'Best_Model' object...\n")
  cat("The 'Best_Model' object has been updated with the following slots:\n")
  cat("- 'performance.result' slot updated.\n")
  
  cat("\n=== Evaluation Completed ===\n")
  cat("Total datasets evaluated:", length(eval_results), "\n")
  
  return(object)
}