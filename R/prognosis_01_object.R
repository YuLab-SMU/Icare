#' PrognosiX Class
#'
#' An S4 class to store survival data, clinical metadata, and various results from the
#' prognostic analysis pipeline, including univariate analysis, feature selection,
#' survival modeling, and subgroup risk assessment.
#' 
#' @import methods
#' @slot clean.data Data frame.
#' @slot info.data Data frame.
#' @slot time_col Any.
#' @slot status_col Any.
#' @slot survival.data Data frame.
#' @slot baseline.table Any.
#' @slot variable.types List.
#' @slot survival.var List.
#' @slot sub.data Data frame.
#' @slot univariate.analysis List.
#' @slot split.data List.
#' @slot feature.result List.
#' @slot filtered.set List.
#' @slot survival.model Any.
#' @slot best.model List.
#' @slot subgroup.risk List.
#' @export
#' @aliases PrognosiX-class
#' @aliases PrognosiX
#' @exportClass PrognosiX
#' @examples
#' \dontrun{
#'   # Create a PrognosiX object
#'   obj <- new("PrognosiX",
#'              clean.data = data.frame(gene1 = c(1, 2, 3)),
#'              info.data = data.frame(time = c(10, 20, 30), status = c(1, 0, 1)),
#'              time_col = "time",
#'              status_col = "status")
#' }
PrognosiX <- setClass(
  Class = 'PrognosiX',
  slots = c(
    clean.data = 'data.frame',
    info.data = 'data.frame',
    time_col = 'ANY',
    status_col = 'ANY',
    survival.data ='data.frame',
    baseline.table = 'ANY',
    variable.types = 'list',
    survival.var = 'list',
    sub.data = 'data.frame',
    univariate.analysis= 'list',
    split.data = 'list',
    feature.result = 'list',
    filtered.set = 'list',
    survival.model = 'ANY',
    best.model = 'list',
    subgroup.risk = 'list'
  ),
  prototype = list(
    clean.data = data.frame(),
    info.data = data.frame(),
    time_col = NULL,
    status_col = NULL,
    survival.data = data.frame(),
    baseline.table = NULL,
    variable.types = list(),
    survival.var = list(),
    sub.data = data.frame(),
    univariate.analysis= list(),
    split.data = list(),
    filtered.set = list(),
    feature.result = list(),
    survival.model = NULL,
    best.model = list(),
    subgroup.risk = list()
  ),
  validity = function(object) {
    # Check data frame slots
    if (!is.data.frame(object@clean.data)) {
      return("Slot 'clean.data' must be a data.frame")
    }
    if (!is.data.frame(object@info.data)) {
      return("Slot 'info.data' must be a data.frame")
    }
    if (!is.data.frame(object@survival.data)) {
      return("Slot 'survival.data' must be a data.frame")
    }
    if (!is.data.frame(object@sub.data)) {
      return("Slot 'sub.data' must be a data.frame")
    }
    
    # Check character columns
    if (!is.null(object@time_col) && !is.character(object@time_col)) {
      return("Slot 'time_col' must be a character string or NULL")
    }
    if (!is.null(object@status_col) && !is.character(object@status_col)) {
      return("Slot 'status_col' must be a character string or NULL")
    }
    
    # Check list slots
    if (!is.list(object@variable.types)) {
      return("Slot 'variable.types' must be a list")
    }
    if (!is.list(object@survival.var)) {
      return("Slot 'survival.var' must be a list")
    }
    if (!is.list(object@univariate.analysis)) {
      return("Slot 'univariate.analysis' must be a list")
    }
    if (!is.list(object@split.data)) {
      return("Slot 'split.data' must be a list")
    }
    if (!is.list(object@feature.result)) {
      return("Slot 'feature.result' must be a list")
    }
    if (!is.list(object@filtered.set)) {
      return("Slot 'filtered.set' must be a list")
    }
    if (!is.list(object@best.model)) {
      return("Slot 'best.model' must be a list")
    }
    if (!is.list(object@subgroup.risk)) {
      return("Slot 'subgroup.risk' must be a list")
    }
    
    # Check survival data has required columns
    if (nrow(object@survival.data) > 0 && !is.null(object@time_col)) {
      if (!object@time_col %in% colnames(object@survival.data)) {
        return(sprintf("Time column '%s' not found in survival.data", object@time_col))
      }
    }
    if (nrow(object@survival.data) > 0 && !is.null(object@status_col)) {
      if (!object@status_col %in% colnames(object@survival.data)) {
        return(sprintf("Status column '%s' not found in survival.data", object@status_col))
      }
    }
    
    # Check row name consistency
    if (nrow(object@clean.data) > 0 && nrow(object@info.data) > 0) {
      if (!setequal(rownames(object@clean.data), rownames(object@info.data))) {
        return("Row names of 'clean.data' and 'info.data' should match")
      }
    }
    
    return(TRUE)
  }
)


`%||%` <- function(a,b) if(!is.null(a)&&length(a)>0) a else b

# Helper to save plots - renamed as requested
save_plot_sur <- function(p, name, w=8, h=6) {
  if(inherits(p, "ggplot")) ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, w, h)
}


#' Create PrognosiX Object
#'
#' Constructs a PrognosiX S4 object. This function ensures that:
#' \itemize{
#'   \item The time and status columns are stored in \code{info.data}.
#'   \item \code{clean.data} contains only numeric feature columns.
#'   \item All other columns (clinical, IDs, etc.) are preserved in \code{info.data}.
#' }
#' It handles input from raw data frames or by converting existing \code{Stat},
#' \code{Subtyping}, or \code{PrognosiX} objects.
#'
#' @param clean.data Data frame of numeric features (samples x features).
#' @param info.data Data frame of metadata (including time, status, clinical vars).
#' @param time_col Character; name of the survival time column (default "time").
#' @param status_col Character; name of the event status column (default "status").
#' @param baseline.table Optional baseline table (reserved).
#' @param variable.types List of variable types (reserved).
#' @param survival.var List of survival variables (reserved).
#' @param survival.data Pre-computed survival data (if NULL, built automatically).
#' @param sub.data Data frame for subgroup analysis (reserved).
#' @param univariate.analysis Results of univariate analysis (reserved).
#' @param split.data Train/test split info (reserved).
#' @param feature.result Feature selection results (reserved).
#' @param filtered.set Filtered data sets (reserved).
#' @param survival.model Fitted survival model (reserved).
#' @param best.model Best model info (reserved).
#' @param subgroup.risk Subgroup risk results (reserved).
#' @param object Optional input object of class \code{Stat}, \code{Subtyping}, or \code{PrognosiX}.
#' @return A \code{PrognosiX} S4 object.
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' veteran$celltype <- as.character(veteran$celltype)
#' prog <- CreatePrognosiXObject(
#'   clean.data = veteran[, c("karno", "diagtime", "age")],
#'   info.data = veteran[, c("time", "status", "celltype", "trt")],
#'   time_col = "time", status_col = "status"
#' )
#' class(prog)
#' }
CreatePrognosiXObject <- function(
    clean.data = NULL,
    info.data = data.frame(),
    time_col = "time",
    status_col = "status",
    baseline.table = NULL,
    variable.types = list(),
    survival.var = list(),
    survival.data = data.frame(),
    sub.data = data.frame(),
    univariate.analysis = list(),
    split.data = list(),
    feature.result = list(),
    filtered.set = list(),
    survival.model = NULL,
    best.model = list(),
    subgroup.risk = list(),
    object = NULL
) {
  
  # ---- 1. Input validation ----
  if (is.null(clean.data) && is.null(object)) {
    stop("At least one of 'clean.data' or 'object' must be provided.")
  }
  
  # ---- 2. Extract data from object if given ----
  if (!is.null(object)) {
    cat("Extracting data from provided object...\n")
    if (!inherits(object, c("Stat", "PrognosiX", "Subtyping"))) {
      stop("'object' must be of class 'Stat', 'PrognosiX', or 'Subtyping'.")
    }
    if (inherits(object, "Stat")) {
      clean.data <- ExtractCleanData(object)
      info.data <- ExtractInfoData(object)
      if (is.null(clean.data) || nrow(clean.data) == 0) {
        stop("Failed to extract valid clean data from 'Stat' object.")
      }
      if (is.null(info.data) || nrow(info.data) == 0) {
        info.data <- data.frame(row.names = rownames(clean.data))
        cat("info.data was created from clean.data with", nrow(info.data), "rows.\n")
      }
    } else if (inherits(object, "PrognosiX")) {
      clean.data <- object@clean.data
      info.data <- object@info.data
      sub.data <- object@sub.data
      time_col <- object@time_col
      status_col <- object@status_col
    } else if (inherits(object, "Subtyping")) {
      clean.data <- if (!is.null(slot(object, "clustered.data")) && nrow(slot(object, "clustered.data")) > 0) {
        slot(object, "clustered.data")
      } else {
        slot(object, "clean.data")
      }
      info.data <- slot(object, "info.data")
    }
  }
  
  if (is.null(clean.data) || nrow(clean.data) == 0) {
    stop("'clean.data' must be provided and not empty.")
  }
  
  # ---- 3. Ensure info.data exists and has row names ----
  if (nrow(info.data) == 0) {
    info.data <- data.frame(row.names = rownames(clean.data))
  }
  
  # ---- 4. Synchronise row names ----
  common <- intersect(rownames(clean.data), rownames(info.data))
  if (length(common) == 0) {
    warning("No matching row names between clean.data and info.data. Recreating info.data from clean.data.")
    info.data <- data.frame(row.names = rownames(clean.data))
  } else {
    if (length(common) < nrow(clean.data)) {
      warning("Subsetting clean.data to match info.data rows.")
      clean.data <- clean.data[common, , drop = FALSE]
    }
    info.data <- info.data[rownames(clean.data), , drop = FALSE]
  }
  
  # ---- 5. Extract time and status columns ----
  # Define helper to find columns either in clean or info
  .find_col <- function(col_name) {
    if (col_name %in% colnames(info.data)) {
      return(list(source = "info", col = col_name))
    } else if (col_name %in% colnames(clean.data)) {
      return(list(source = "clean", col = col_name))
    } else {
      return(NULL)
    }
  }
  
  time_info <- .find_col(time_col)
  status_info <- .find_col(status_col)
  
  if (is.null(time_info) || is.null(status_info)) {
    stop(sprintf("Columns '%s' and/or '%s' not found in clean.data or info.data.", time_col, status_col))
  }
  
  # Move time/status to info.data if they are currently in clean.data
  if (time_info$source == "clean") {
    info.data[[time_col]] <- clean.data[[time_col]]
    clean.data[[time_col]] <- NULL
  }
  if (status_info$source == "clean") {
    info.data[[status_col]] <- clean.data[[status_col]]
    clean.data[[status_col]] <- NULL
  }
  
  # ---- 6. Move all non-numeric columns from clean.data to info.data ----
  # Identify non-numeric columns (character, factor, logical, etc.)
  non_numeric <- names(which(!sapply(clean.data, is.numeric)))
  if (length(non_numeric) > 0) {
    for (col in non_numeric) {
      info.data[[col]] <- clean.data[[col]]
      clean.data[[col]] <- NULL
    }
    cat("Moved non-numeric columns to info.data:", paste(non_numeric, collapse = ", "), "\n")
  }
  
  # ---- 7. Ensure clean.data contains only numeric columns ----
  if (ncol(clean.data) == 0) {
    warning("clean.data has no columns after moving non-numeric data. This may be acceptable if you only have clinical metadata.")
  }
  
  # ---- 8. Prepare column names (sanitise) ----
  .prepare_colnames <- function(df) {
    if (nrow(df) == 0 || ncol(df) == 0) return(df)
    # Use make.names to ensure valid R variable names
    colnames(df) <- make.names(colnames(df), unique = TRUE)
    df
  }
  clean.data <- .prepare_colnames(clean.data)
  info.data  <- .prepare_colnames(info.data)
  
  # ---- 9. Build survival.data (features + time/status) ----
  # Ensure time and status are numeric
  info.data[[time_col]] <- as.numeric(as.character(info.data[[time_col]]))
  info.data[[status_col]] <- as.numeric(as.character(info.data[[status_col]]))
  
  # Remove rows with missing time/status
  valid <- complete.cases(info.data[, c(time_col, status_col)])
  if (any(!valid)) {
    warning(sprintf("Removing %d rows with missing time or status.", sum(!valid)))
    clean.data <- clean.data[valid, , drop = FALSE]
    info.data <- info.data[valid, , drop = FALSE]
  }
  
  # Combine features and survival columns
  survival.data <- cbind(clean.data, info.data[, c(time_col, status_col), drop = FALSE])
  
  # ---- 10. Create the S4 object ----
  obj <- new(
    Class = "PrognosiX",
    clean.data = clean.data,
    info.data = info.data,
    sub.data = sub.data,
    univariate.analysis = univariate.analysis,
    time_col = time_col,
    status_col = status_col,
    survival.data = survival.data,
    baseline.table = baseline.table,
    variable.types = variable.types,
    survival.var = survival.var,
    split.data = split.data,
    feature.result = feature.result,
    filtered.set = filtered.set,
    survival.model = survival.model,
    best.model = best.model,
    subgroup.risk = subgroup.risk
  )
  
  cat("PrognosiX object created successfully.\n")
  cat(sprintf("  Features: %d | Samples: %d\n", ncol(clean.data), nrow(clean.data)))
  cat(sprintf("  Metadata columns: %d (including time/status)\n", ncol(info.data)))
  return(obj)
}
