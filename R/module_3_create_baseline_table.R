#' Sur Gaze Analysis
#'
#' @param object SurObj.
#' @param formula Formula.
#' @param status_col Status column.
#' @param digits Digits.
#' @param show.p Show P.
#' @param gaze_method Method.
#' @param response_var Response var.
#' @param save_word Save word.
#' @param use_subgroup_data Use subgroup data.
#' @param save_dir Save dir.
#' @export
Sur_gaze_analysis <- function(object,
                              formula = NULL,
                              status_col = "status",
                              digits = 1,
                              show.p = TRUE,
                              gaze_method = 3,
                              response_var = "group",
                              save_word = TRUE,   
                              use_subgroup_data = FALSE,
                              save_dir = here("SurObj", "gaze_baseline")) {
  
  cat("Starting gaze analysis...\n")
  
  if (inherits(object, "SurObj")) {
    if (use_subgroup_data) {
      gaze_data <- slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
    } else {
      gaze_data <- slot(object, "survival.data")
      cat("Using original survival data...\n")
    } 
  }
  if (is.null(gaze_data) || nrow(gaze_data) == 0) {
    stop("Input data is empty or NULL")
  }
  if (!is.null(response_var)) {
    if (!response_var %in% colnames(gaze_data)) {
      stop("Response variable '", response_var, "' not found in data columns: ", 
           paste(colnames(gaze_data), collapse=", "))
    }
    if (!is.factor(gaze_data[[response_var]])) {
      cat("Converting response variable '", response_var, "' to factor...\n", sep="")
      gaze_data[[response_var]] <- as.factor(gaze_data[[response_var]])
    }
  }
  
  if (is.null(formula)) {
    if (!is.null(response_var)) {
      formula <- as.formula(paste(response_var, "~ ."))
      cat("Using default formula: ", deparse(formula), "\n", sep="")
    } else {
      stop("Either 'formula' or 'response_var' must be specified.")
    }
  }
  
  tryCatch({
    cat("Performing gaze analysis (method ", gaze_method, ")...\n", sep = "")
    
    if (!exists("gaze_analysis")) {
      stop("gaze_analysis function not found. Please check if the required package is loaded.")
    }
    
    result <- gaze_analysis(
      data = gaze_data,
      formula = formula,
      group_cols = response_var,
      digits = digits,
      show.p = show.p,
      gaze_method = gaze_method,
      save_word = FALSE
    )
    
    if (is.null(result)) {
      stop("Gaze analysis returned NULL result")
    }
    
    if (save_word) {
      cat("Generating Word report...\n")
      
      if (!requireNamespace("officer", quietly = TRUE)) {
        stop("Package 'officer' required for Word output. Please install it.")
      }
      if (!requireNamespace("flextable", quietly = TRUE)) {
        stop("Package 'flextable' required for Word output. Please install it.")
      }
      
      doc <- officer::read_docx()
      
      doc <- tryCatch({
        doc %>%
          officer::body_add_par("Gaze Analysis Results", style = "heading 1") %>%
          flextable::body_add_flextable(result)
      }, error = function(e) {
        warning("Failed to add table to Word document: ", e$message)
        doc %>%
          officer::body_add_par("Gaze Analysis Results", style = "heading 1") %>%
          officer::body_add_par("Table could not be added due to formatting issues")
      })
      
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
        cat("Created output directory: '", save_dir, "'\n", sep="")
      }
      
      word_filename <- file.path(save_dir, "gaze_analysis_results.docx")
      print(doc, target = word_filename)
      cat("Results saved to: '", word_filename, "'\n", sep="")
    }
    
    if (inherits(object, "SurObj")) {
      cat("Updating SurObj object with gaze results...\n")
      object@baseline.table <- result
      cat("The 'SurObj' object has been updated with the following slots:\n")
      cat("- 'baseline.table' slot updated.\n")
      return(object)
    }
    
    return(result)
    
  }, error = function(e) {
    stop("Gaze analysis failed: ", e$message)
  })
}
