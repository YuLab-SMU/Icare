#' SurObj Class
#'
#' @import methods
#' @slot clean.data Data frame.
#' @slot info.data Data frame.
#' @slot clustered.data Data frame.
#' @slot survival.data Data frame.
#' @slot sub.data Data frame.
#' @slot baseline.table Any.
#' @slot univariate.analysis List.
#' @slot survival.analysis List.
#' @slot rcss.results Any.
#' @export
SurObj <- setClass(
  Class = 'SurObj',
  slots = c(
    clean.data = 'data.frame',
    info.data = 'data.frame',
    clustered.data = 'data.frame',
    survival.data = 'data.frame',
    sub.data = 'data.frame',
    baseline.table = 'ANY',
    univariate.analysis = 'list',
    survival.analysis = 'list',
    rcss.results = 'ANY'
  ),
  prototype = list(
    clean.data = data.frame(),
    info.data = data.frame(),
    clustered.data = data.frame(),
    survival.data = data.frame(),
    sub.data = data.frame(),
    baseline.table = NULL,
    univariate.analysis = list(),
    survival.analysis = list(),
    rcss.results = NULL
  )
)

#' Create SurObj
#'
#' @param clean.data Clean data.
#' @param info.data Info data.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param clustered.data Clustered data.
#' @param survival.data Survival data.
#' @param sub.data Sub data.
#' @param baseline.table Baseline table.
#' @param univariate.analysis Univariate analysis.
#' @param survival.analysis Survival analysis.
#' @param rcss.results RCSS results.
#' @param data_type Data type.
#' @param object Input object.
#' @export
CreateSurObj <- function(
    clean.data = NULL,
    info.data = data.frame(),
    time_col = "time",  
    status_col = "status",  
    clustered.data = data.frame(),
    survival.data = data.frame(),
    sub.data = data.frame(),
    baseline.table = NULL,
    univariate.analysis = list(),
    survival.analysis = list(),
    rcss.results = NULL,
    data_type = "train",
    object = NULL
) {
  
  if (is.null(clean.data) && is.null(object) && is.null(clustered.data) && is.null(survival.data)) {
    stop("At least one data source must be provided from: 'clean.data', 'object', 'clustered.data', or 'survival.data'")
  }
  
  if (!is.null(object)) {
    cat("Extracting data from provided object...\n")
    
    if (!inherits(object, c("Stat", "SurObj", "Subtyping", "PrognosiX"))) {
      stop("Input object must be instance of class 'Stat', 'SurObj', 'Subtyping' or 'PrognosiX'")
    }
    
    if (inherits(object, "Stat")) {
      clean.data <- ExtractCleanData(object)
      info.data <- ExtractInfoData(object)
      
      if (is.null(clean.data) || nrow(clean.data) == 0) {
        stop("Failed to extract valid clean data from Stat object")
      }
      
      if (is.null(info.data) || nrow(info.data) == 0) {
        info.data <- data.frame(row.names = rownames(clean.data))
        cat("Created empty info.data with", nrow(info.data), "rows\n")
      }
      
    } else if (inherits(object, "SurObj")) {
      clean.data <- object@clean.data
      info.data <- object@info.data
      clustered.data <- object@clustered.data
      survival.data <- object@survival.data
      sub.data <- object@sub.data
      baseline.table <- object@baseline.table
      univariate.analysis <- object@univariate.analysis
      survival.analysis <- object@survival.analysis
      rcss.results <- object@rcss.results
      
    } else if (inherits(object, "Subtyping")) {
      clustered_data <- slot(object, "clustered.data")
      info.data <- object@info.data
      if (!is.null(clustered_data) && nrow(clustered_data) > 0) {
        clean.data <- clustered_data
        cat("Using clustered.data (", nrow(clustered_data), " rows) as clean.data\n")
      } else {
        clean_data <- slot(object, "clean.data")
        if (!is.null(clean_data) && nrow(clean_data) > 0) {
          clean.data <- clean_data
          cat("Using clean.data (", nrow(clean_data), " rows)\n")
        } else {
          stop("No valid data found in Subtyping object")
        }
      }
      
      info.data <- slot(object, "info.data")
      if (is.null(info.data) || nrow(info.data) == 0) {
        warning("info.data is empty in Subtyping object")
        info.data <- data.frame(row.names = rownames(clean.data))
      }
    } else if(inherits(object, "PrognosiX")) {
      clean.data  <- object@split.data[[data_type]] 
      info.data <- object@info.data[[data_type]] 
      time_col <- object@time_col
      status_col <- object@status_col
      survival.data <- object@survival.data[[data_type]] 
      sub.data <- object@sub.data[[data_type]]
    }
  }
  
  prepare_data <- function(data, data_name) {
    if (is.null(data) || nrow(data) == 0) {
      return(data.frame())
    }
    
    if (anyDuplicated(rownames(data))) {
      warning("Duplicate row names in ", data_name, " - making unique")
      rownames(data) <- make.unique(rownames(data))
    }
    
    if (anyDuplicated(colnames(data))) {
      warning("Duplicate column names in ", data_name, " - making unique")
      colnames(data) <- make.unique(colnames(data))
    }
    
    return(data)
  }
  
  # 准备各数据
  clean.data <- prepare_data(clean.data, "clean.data")
  info.data <- prepare_data(info.data, "info.data")
  sub.data <- prepare_data(sub.data, "sub.data")
  
  # 确保clean.data存在
  if (is.null(clean.data) || nrow(clean.data) == 0) {
    stop("No valid clean.data available")
  }
  
  if (nrow(survival.data) == 0) {
    if (!all(c(time_col, status_col) %in% colnames(info.data))) {
      if (all(c(time_col, status_col) %in% colnames(clean.data))) {
        info.data <- cbind(info.data, clean.data[, c(time_col, status_col)])
        clean.data <- clean.data[, !colnames(clean.data) %in% c(time_col, status_col)]
      } else {
        stop("Time/status columns not found in data")
      }
    }
    
    # 标准化列名
    colnames(info.data)[colnames(info.data) == time_col] <- "time"
    colnames(info.data)[colnames(info.data) == status_col] <- "status"
    
    survival.data <- cbind(clean.data, info.data[, c("time", "status")])
    survival.data <- survival.data[complete.cases(survival.data[, c("time", "status")]), ]
  }
  
  SurObj<-new(
    Class = 'SurObj',
    clean.data = clean.data,
    info.data = info.data,
    clustered.data = clustered.data,
    survival.data = survival.data,
    sub.data = sub.data,
    baseline.table = baseline.table,
    univariate.analysis = univariate.analysis,
    survival.analysis = survival.analysis,
    rcss.results = rcss.results
  )
  
  cat("SurObj object created successfully.\n")
  return(SurObj)
}
