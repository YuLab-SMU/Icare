#' Best_Model Class
#'
#' An S4 class to store information about the best performing model from training,
#' including model details, performance metrics, and interpretation results.
#'
#' @slot best.model.type Character specifying the type of best model (e.g., "randomForest")
#' @slot filtered.set List containing filtered training and testing datasets
#' @slot group_col Character specifying the column name of the response variable
#' @slot best.model The actual best model object (can be any model type)
#' @slot best_threshold List containing optimal classification threshold information
#' @slot performance.result List containing model performance metrics
#' @slot shap.result List containing SHAP (SHapley Additive exPlanations) analysis results
#' @slot feature.importance List containing feature importance results
#' @slot cm.result List containing confusion matrix results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a Best_Model object
#' best_model <- new("Best_Model",
#'                  best.model.type = "randomForest",
#'                  group_col = "response")
#' }
Best_Model <- setClass(
  Class = 'Best_Model',
  slots = c(
    best.model.type = 'character',  
    filtered.set = 'list',
    group_col = 'character', 
    process.info = 'list',
    best.model = 'ANY', 
    best_threshold = 'list',
    performance.result = 'list',
    shap.result = 'list',
    feature.importance = 'list',
    cm.result = 'list'
  ),
  prototype = list(
    best.model.type = "",
    filtered.set = list(),
    process.info = list(),
    group_col = "group",
    best.model = NULL,  
    performance.result = list(),
    shap.result = list(),
    feature.importance = list(),
    cm.result = list()
  )
)

#' Create Best_Model Object
#'
#' Constructs a Best_Model object either from scratch or by extracting information
#' from a Train_Model object. Handles both direct input of model components and
#' extraction from existing objects.
#'
#' @param filtered.set Optional list containing filtered training and testing datasets.
#'                    Must include 'training' and 'testing' elements.
#' @param object Optional Train_Model or Best_Model object to extract information from.
#' @param group_col Character specifying the response variable column name.
#'                 Defaults to "group" if not provided.
#' @param best_model Optional pre-trained model object to include.
#' @param best.model.type Optional character specifying the model type.
#' @param ... Additional arguments passed to methods.
#'
#' @return A Best_Model object containing model information and results.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From Train_Model object:
#' object_best <- CreateBestModel(object = object_model)
#' }
CreateBestModel <- function(
    filtered.set = NULL,
    object = NULL,
    group_col = NULL,
    best_model = NULL,
    best.model.type = NULL,
    process.info = list(),
    ...
) {
  if (!is.null(filtered.set) && !is.null(object)) {
    stop("Only one of 'filtered.set' and 'object' should be provided.")
  }
  
  if (is.null(filtered.set) && is.null(object)) {
    stop("Either 'filtered.set' or 'object' must be provided.")
  }
  
  model_type <- if (!is.null(best.model.type)) best.model.type else ""
  
  if (!is.null(object)) {
    if (inherits(object, "Train_Model")) {
      filtered.set <- object@filtered.set
      group_col <- object@group_col
      
      if (length(object@best.model.result) > 0) {
        model_type <- object@best.model.result[["model_type"]]
        best_model <- object@best.model.result[["model"]][1]
        
      } else {
        warning("No best model result found in the Train_Model object")
        best_model <- list()
      }
      mapping_info <- if (!is.null(object@process.info[["mapping"]])) {
        object@process.info[["mapping"]]
      }
      
      missing_removal <- if (!is.null(object@process.info[["missing_removal"]])) {
        object@process.info[["missing_removal"]]
      }
      
      missing_info <- if (!is.null(object@process.info[["missing_info"]])) {
        object@process.info[["missing_info"]]
      }
      
      outlier_info <- if (!is.null(object@process.info[["outlier_result"]][["handle_outliers_train"]])) {
        object@process.info[["outlier_result"]]
      }
      
      normalization_info <- if (!is.null(object@split.scale.data[["normalization_info"]])) {
        object@split.scale.data[["normalization_info"]]
      }      
      process.info <- list(
        mapping_info = mapping_info,
        missing_removal = missing_removal,
        missing_info = missing_info,
        outlier_info = outlier_info,
        normalization_info=normalization_info
      )
    } else if (inherits(object, "Best_Model")) {
      
      cat("Input is already a Best_Model object")
      return(object)
    } else {
      stop("Input object must be either 'Train_Model' or 'Best_Model' class")
    }
  } else {
    if (!is.list(filtered.set) || 
        !all(c("training", "testing") %in% names(filtered.set))) {
      stop("filtered.set must be a list with 'training' and 'testing' elements")
    }
    
    if (is.null(group_col)) {
      group_col <- "group"
      warning("No group_col provided, using default 'group'")
    }
    
    if (!is.null(best_model) && !is.list(best_model)) {
      best_model <- list(best_model)
      if (!is.null(best.model.type)) {
        names(best_model) <- best.model.type
      }
    }
  }
  
  Best_Model_instance <- new(
    'Best_Model',
    best.model.type = model_type,
    filtered.set = filtered.set,
    group_col = group_col,
    process.info =process.info,
    best.model = if (!is.null(best_model)) best_model else list()
  )
  
  cat("Best_Model object created successfully.\n")
  return(Best_Model_instance)
}