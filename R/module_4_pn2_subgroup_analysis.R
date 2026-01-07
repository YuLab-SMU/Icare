#' Perform Subgroup Analysis for PrognosiX
#'
#' @param object PrognosiX object.
#' @param status_col Status col.
#' @param id_col ID col.
#' @param time_col Time col.
#' @param max_unique_values Max unique values.
#' @param methods Method.
#' @param seed Seed.
#' @param save_dir Save dir.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
PrognosiX_subgroup_analysis <- function(object,
                                      status_col="status",
                                      id_col = "id",
                                      time_col = "time",
                                      max_unique_values  = 5,
                                      methods = "roc",
                                      seed =123,
                                      save_dir = here::here("StatObject","Data"),
                                      save_data = TRUE,
                                      csv_filename = "clean_data.csv" ) {
  
  cat("Performing subgroup analysis...\n")
  
  if (inherits(object, 'PrognosiX')) {
    data <- slot(object, "survival.data")
    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("The survival.data in the PrognosiX object is empty.")
  }
  
  cat("Diagnosing variable types...\n")
  diagnostics <- diagnose_variable_type(data, 
                                        max_unique_values=max_unique_values)
  
  cat("Processing variables by method: ", methods, "\n")
  
  if ("maxstat" %in% methods) {
    new_data <- find_maxstat_cutoff(data, 
                                    variables=diagnostics$numeric_vars,
                                    time_col =time_col,
                                    status_col = status_col)  } 
  else if ("roc" %in% methods) {
    new_data <- find_roc_cutoff(data, 
                                variables=diagnostics$numeric_vars, 
                                status_col = status_col, 
                                
                                time_col=time_col)  } else {
                                  stop("Invalid method specified. Choose either 'roc' or 'maxstat'.")
                                }
  
  variable_types <- diagnose_variable_type(new_data, max_unique_values=max_unique_values)
  new_data <- convert_two_class_to_binary(new_data)
  new_data <- convert_variables(new_data, variable_types = variable_types)
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(new_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  
  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object...\n")
    
    object@sub.data <- new_data
    
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'sub.data' slot updated.\n")
    
    return(object)
  }
  
  return(new_data)
}
