library(autoReg)
library(dplyr)
library(gridExtra)
library(grid)
library(flextable)
library(here)
library(webshot)
library(ggprism)
library(officer)
gaze_analysis <- function(data,
                          formula = NULL,
                          group_cols = NULL,
                          digits = 1,
                          show.p = TRUE,
                          gaze_method = 3,
                          save_word = TRUE,
                          save_dir = here("PrognosiX", "gaze_baseline")) {
  
  if (!is.data.frame(data)) stop("The input 'data' must be a data frame.")
  
  if (is.null(formula)) {
    if (!is.null(group_cols)) {
      if (length(group_cols) == 1 && group_cols %in% colnames(data)) {
        formula <- as.formula(paste(group_cols, "~ ."))
        cat("Using formula with one group:", deparse(formula), "\n")
      } else if (length(group_cols) > 1 && all(group_cols %in% colnames(data))) {
        formula <- as.formula(paste(paste(group_cols, collapse = " + "), "~ ."))
        cat("Using formula with multiple groups:", deparse(formula), "\n")
      } else {
        stop("Group columns not found in data.")
      }
    } else {
      formula <- as.formula("~ .")
      cat("Using default formula: ~ .\n")
    }
  } else if (!inherits(formula, "formula")) {
    stop("The input 'formula' must be a valid formula.")
  }
  
  if (!is.numeric(digits) || digits < 0 || digits != as.integer(digits))
    stop("The input 'digits' must be a non-negative integer.")
  
  if (!is.logical(show.p))
    stop("The input 'show.p' must be a logical value (TRUE or FALSE).")
  
  if (!is.numeric(gaze_method) || gaze_method < 1 || gaze_method > 5 || gaze_method != as.integer(gaze_method))
    stop("The input 'gaze_method' must be an integer between 1 and 5.")
  
  tryCatch({
    cat("Running gaze analysis with method:", gaze_method, "\n")
    result <- gaze(formula, data, digits = digits, show.p = show.p, method = gaze_method)
    
    cat("Result type:", class(result), "\n")
    if (is.data.frame(result) || is.matrix(result)) {
      result <- myft(result)
      cat("Gaze analysis completed successfully.\n")
      
      if (save_word) {
        doc <- read_docx()
        doc <- doc %>%
          body_add_flextable(result) %>%
          body_add_par("Gaze Analysis Results", style = "heading 1")
        
        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }
        
        word_filename <- file.path(save_dir, "gaze_analysis.docx")
        print(doc, target = word_filename)
        cat("Word file saved to:", word_filename, "\n")
      }
      
      return(result)
    } else {
      stop("The result is not a data frame or matrix.")
    }
  }, error = function(e) {
    stop("An error occurred while performing the gaze analysis: ", e$message)
  })
}

Prognos_gaze_analysis <- function(object,
                                  formula = NULL,
                                  status_col = "status",
                                  digits = 1,
                                  show.p = TRUE,
                                  gaze_method = 3,
                                  response_var = NULL,
                                  save_word = TRUE,   
                                  save_dir = here("PrognosiX", "gaze_baseline"),
                                  use_subgroup_data = FALSE) {
  
  cat("Starting Prognos_gaze_analysis function...\n")
  
  if (inherits(object, "PrognosiX")) {
    if (use_subgroup_data) {
      gaze_data <- slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
    } else {
      gaze_data <- slot(object, "survival.data")
      cat("Using original survival data...\n")
    }
    status_col <- slot(object, "status_col")
    time_col <- slot(object, "time_col")
  } else if (is.data.frame(object)) {
    gaze_data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }
  
  if (is.null(gaze_data) || nrow(gaze_data) == 0) {
    stop("No valid data found in the input.")
  }
  
  if (!is.null(response_var)) {
    if (!response_var %in% colnames(gaze_data)) {
      stop("Response variable not found in data.")
    }
    if (!is.factor(gaze_data[[response_var]])) {
      cat("Converting response variable to factor...\n")
      gaze_data[[response_var]] <- as.factor(gaze_data[[response_var]])
    }
  }
  
  if (is.null(formula)) {
    formula <- as.formula(paste(response_var, "~ ."))
    cat("Using formula:", deparse(formula), "\n")
  }
  
  tryCatch({
    cat("Running gaze analysis with method:", gaze_method, "\n")
    result <- gaze_analysis(data = gaze_data,
                            formula = formula,
                            group_cols = response_var,
                            digits = digits,
                            show.p = show.p,
                            gaze_method = gaze_method,
                            save_word = FALSE,  # Disable saving in gaze_analysis
                            save_dir = save_dir)
    
    if (save_word) {
      cat("Saving results as Word document...\n")
      doc <- read_docx()
      doc <- doc %>%
        body_add_flextable(result) %>%
        body_add_par("Gaze Analysis Results", style = "heading 1")
      
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }
      
      word_filename <- file.path(save_dir, "gaze_analysis.docx")
      print(doc, target = word_filename)
      cat("Word file saved to:", word_filename, "\n")
    }
    
    if (inherits(object, "PrognosiX")) {
      cat("Updating 'PrognosiX' object...\n")
      object@baseline.table <- result
      cat("The 'PrognosiX' object has been updated.\n")
      return(object)
    }
    
    return(result)
    
  }, error = function(e) {
    stop("An error occurred during gaze analysis: ", e$message)
  })
}