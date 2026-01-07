#' Prepare Data For PrognosiX
#'
#' @param object PrognosiX object.
#' @param group_col Group col.
#' @param max_unique_values Max unique values.
#' @export
PrepareDataForPrognosiX <- function(object ,
                                    group_col = NULL,
                                    max_unique_values = 5) {
  cat("Starting transformation of factors to dummy variables... \n")
  
  if (inherits(object, 'PrognosiX')) {
    cat("Input is a 'PrognosiX' object. Extracting the data slot... \n")
    survival_data<-slot(object, "survival.data")
    status_col <- slot(object, "status_col")
    sub_data<-slot(object, "sub.data")
  }
  
  if (nrow(survival_data) == 0 && nrow(sub_data) == 0) {
    cat("Both of the extracted datasets are empty. Skipping transformation.\n")
    return(object)
  }
  
  survival_data <- convert_to_numeric(survival_data)
  sub_data <- convert_to_numeric(sub_data)
  
  survival_data <- convert_two_class_to_binary(survival_data)
  sub_data <- convert_two_class_to_binary(sub_data)
  
  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object with transformed data... \n")
    slot(object, "survival.data") <- survival_data
    slot(object, "sub.data") <-sub_data
    cat("The 'PrognosiX' object has been updated with the following:\n")
    cat("- 'survival.data' slot updated.\n")
    cat("- 'sub.data' slot updated.\n")
    
    return(object)
  }
  
  cat("Returning transformed data frame.  \n")
  return(data)
}
