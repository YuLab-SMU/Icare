#' Extract Info Data
#'
#' @param object An object containing info.data slot.
#' @return A data.frame containing the info.data slot, or NULL if extraction fails.
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'sub_obj_test' is a valid Subtyping object with an info.data slot
#'   info <- ExtractInfoData(sub_obj_test)
#' }
ExtractInfoData <- function(object) {
  
  data <- tryCatch(
    slot(object, "info.data"),
    error = function(e) {
      warning("Error extracting 'info.data' from the object: ", e$message)
      return(NULL)
    }
  )
  return(data)
}

#' Subtyping Class
#'
#' @import methods
#' @slot clean.data Data frame.
#' @slot info.data Data frame.
#' @slot scale.data Data frame.
#' @slot Optimal.cluster Any.
#' @slot cluster.results Any.
#' @slot visualization.results Any.
#' @slot clustered.data Data frame.
#' @slot evaluation_results List.
#' @export
#' @examples
#' \dontrun{
#'   obj <- new("Subtyping",
#'              clean.data = data.frame(x = 1:3, y = 4:6),
#'              info.data = data.frame(id = letters[1:3]),
#'              scale.data = data.frame(x = c(0, 0.5, 1), y = c(0, 0.5, 1)))
#' }
Subtyping <- setClass(
  Class = 'Subtyping',
  slots = c(
    clean.data = 'data.frame',
    info.data = 'data.frame',
    scale.data = 'data.frame',
    Optimal.cluster = 'ANY',
    cluster.results = 'ANY',
    visualization.results = 'ANY',
    clustered.data = 'data.frame',
    evaluation_results = 'list'
  ),
  prototype = list(
    clean.data = data.frame(),
    info.data = data.frame(),
    scale.data = data.frame(),
    Optimal.cluster = NULL,
    cluster.results = NULL,
    visualization.results = NULL,
    clustered.data = data.frame(),
    evaluation_results = list()
  ),
  validity = function(object) {
    # Check data frame slots
    if (!is.data.frame(object@clean.data)) {
      return("Slot 'clean.data' must be a data.frame")
    }
    if (!is.data.frame(object@info.data)) {
      return("Slot 'info.data' must be a data.frame")
    }
    if (!is.data.frame(object@scale.data)) {
      return("Slot 'scale.data' must be a data.frame")
    }
    if (!is.data.frame(object@clustered.data)) {
      return("Slot 'clustered.data' must be a data.frame")
    }
    
    # Check list slot
    if (!is.list(object@evaluation_results)) {
      return("Slot 'evaluation_results' must be a list")
    }
    
    # Check row name consistency
    if (nrow(object@clean.data) > 0 && nrow(object@info.data) > 0) {
      if (!setequal(rownames(object@clean.data), rownames(object@info.data))) {
        return("Row names of 'clean.data' and 'info.data' must match")
      }
    }
    
    # Check cluster column if clustered.data exists
    if (nrow(object@clustered.data) > 0) {
      if (!"cluster" %in% colnames(object@clustered.data)) {
        return("Slot 'clustered.data' must contain a 'cluster' column")
      }
    }
    
    return(TRUE)
  }
)

#' Create Subtyping Object
#'
#' @param clean.data Clean data.
#' @param info.data Info data.
#' @param scale.data Scale data.
#' @param Optimal.cluster Optimal cluster.
#' @param cluster.results Cluster results.
#' @param visualization.results Visualization results.
#' @param clustered.data Clustered data.
#' @param evaluation_results Evaluation results.
#' @param object Input object (Stat, Subtyping, PrognosiX, Model_data).
#' @param convert_factors Logical; Default FALSE – only enable if you are certain the variable is ordinal or binary.
#' @param na.action Character string specifying how to handle NA values. 
#'   Options are "omit" (remove rows with NA, default for subtyping), 
#'   "allow" (keep NA values), or "error" (stop if NA found).
#' @export
#' @examples
#' \dontrun{
#'   # Create from data frames
#'   clean <- data.frame(gene1 = c(1.2, 3.4, 5.6), gene2 = c(2.3, 4.5, 6.7),
#'                       row.names = c("S1", "S2", "S3"))
#'   info <- data.frame(age = c(25, 30, 35), row.names = c("S1", "S2", "S3"))
#'   obj <- CreateSubtypingObject(clean.data = clean, info.data = info)
#' }
CreateSubtypingObject <- function(
    clean.data = NULL,
    info.data = data.frame(),
    scale.data = data.frame(),
    Optimal.cluster = NULL,
    cluster.results = NULL,
    visualization.results = NULL,
    clustered.data = data.frame(),
    evaluation_results = list(),
    object = NULL,
    convert_factors = FALSE,
    na.action = c("omit", "allow", "error")
) {
  
  if (!requireNamespace("janitor", quietly = TRUE)) {
    message("Installing janitor package for column name cleaning...")
    install.packages("janitor")
  }
  
  na.action <- match.arg(na.action)
  if (is.null(clean.data) && is.null(object)) {
    stop("At least one of 'clean.data' or 'object' must be provided.")
  }
  
  if (!is.null(object)) {
    if (!inherits(object, "Stat") && !inherits(object, "Subtyping") && 
        !inherits(object, "PrognosiX") && !inherits(object, "Model_data") && !inherits(object, "Train_Model")) {
      stop("The 'object' parameter must be an instance of class 'Stat', 'Subtyping', 'Model_data'/'Train_Model' or 'PrognosiX'.")
    }
    
    if (inherits(object, "Stat")) {
      clean.data <- ExtractCleanData(object)   # Note: ExtractCleanData must be defined elsewhere
      info.data <- ExtractInfoData(object)
      scale.data <- object@scale.data
      if (is.null(clean.data)) {
        stop("Failed to extract clean data from the provided 'Stat' object.")
      }
    }
    
    else if (inherits(object, "Subtyping")) {
      clean.data <- object@clean.data
      info.data <- object@info.data
      scale.data <- object@scale.data
      Optimal.cluster <- object@Optimal.cluster
      cluster.results <- object@cluster.results
      visualization.results <- object@visualization.results
      clustered.data <- object@clustered.data
      evaluation_results <- object@evaluation_results
    }
    
    else if (inherits(object, "PrognosiX")) {
      clean.data <- object@clean.data
      info.data <- object@info.data
      scale.data <- object@scale.data
    }
    else if (inherits(object, "Model_data") || inherits(object, "Train_Model")) {
      clean.data <- object@clean.df
      # scale.data might be different in Train_Model
      if(.hasSlot(object, "scale.data")) scale.data <- object@scale.data
    }
  }
  

  if (is.data.frame(clean.data) && ncol(clean.data) > 0) {
    colnames(clean.data) <- janitor::make_clean_names(colnames(clean.data))
  }
  if (is.data.frame(info.data) && ncol(info.data) > 0) {
    colnames(info.data) <- janitor::make_clean_names(colnames(info.data))
  }
  if (is.data.frame(scale.data) && ncol(scale.data) > 0) {
    colnames(scale.data) <- janitor::make_clean_names(colnames(scale.data))
  }
  if (is.data.frame(clustered.data) && ncol(clustered.data) > 0) {
    colnames(clustered.data) <- janitor::make_clean_names(colnames(clustered.data))
  }
  
  ensure_numeric_data <- function(data, data_name, convert_factors = TRUE) {
    if (is.null(data) || nrow(data) == 0) {
      return(data.frame())
    }
    
    row_names <- rownames(data)
    orig_colnames <- colnames(data)
    
    # Helper function to safely convert character to numeric
    .safe_as_numeric <- function(x, col_name) {
      converted <- suppressWarnings(as.numeric(x))
      if (all(is.na(converted)) && any(!is.na(x))) {
        if (any(grepl("%$", x[!is.na(x)]))) {
          x <- gsub("%$", "", x)
          converted <- suppressWarnings(as.numeric(x))
        }
        else if (any(grepl(",", x[!is.na(x)]) && !any(grepl("\\.", x[!is.na(x)])))) {
          x <- gsub(",", ".", x)
          converted <- suppressWarnings(as.numeric(x))
        }
      }
      return(converted)
    }
    
    converted_list <- lapply(seq_along(data), function(i) {
      x <- data[[i]]
      col_name <- orig_colnames[i]
      
      if (is.factor(x)) {
        if (convert_factors) {
          levels <- levels(x)
          if (length(levels) == 2) {
            return(as.numeric(x) - 1)
          } else {
            numeric_levels <- .safe_as_numeric(levels, col_name)
            if (!all(is.na(numeric_levels))) {
              warning(paste("Factor variable '", col_name, "' in", data_name, 
                            "has", length(levels), "levels. Converting level names to numeric values."))
              return(numeric_levels[as.numeric(x)])
            } else {
              warning(paste("Factor variable '", col_name, "' in", data_name, 
                            "has", length(levels), "levels that cannot be converted to numeric. Using integer codes."))
              return(as.numeric(x))
            }
          }
        } else {
          warning(paste("Factor variable '", col_name, "' in", data_name, "was not converted. Use convert_factors=TRUE to convert."))
          return(x)
        }
      } else if (is.character(x)) {
        if (convert_factors) {
          converted <- .safe_as_numeric(x, col_name)
          if (!all(is.na(converted))) {
            return(converted)
          }
          x_factor <- as.factor(x)
          levels <- levels(x_factor)
          if (length(levels) == 2) {
            warning(paste("Character variable '", col_name, "' in", data_name, 
                          "converted to binary numeric (0/1)."))
            return(as.numeric(x_factor) - 1)
          } else {
            warning(paste("Character variable '", col_name, "' in", data_name, 
                          "has", length(levels), "unique values. Converting to numeric codes."))
            return(as.numeric(x_factor))
          }
        } else {
          warning(paste("Character variable '", col_name, "' in", data_name, "was not converted. Use convert_factors=TRUE to convert."))
          return(x)
        }
      } else if (is.numeric(x)) {
        return(x)
      } else if (is.logical(x)) {
        return(as.numeric(x))
      } else {
        warning(paste("Column '", col_name, "' in", data_name, "is of type", class(x)[1], 
                      "and cannot be converted to numeric. It will be kept as-is."))
        return(x)
      }
    })
    
    names(converted_list) <- orig_colnames
    data <- as.data.frame(converted_list, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(data) <- row_names
    
    na_columns <- colnames(data)[apply(data, 2, function(x) all(is.na(x)))]
    if (length(na_columns) > 0) {
      warning(paste("The following columns in", data_name, "were entirely NA after conversion and will be removed:", paste(na_columns, collapse = ", ")))
      data <- data[, !colnames(data) %in% na_columns, drop = FALSE]
    }
    
    non_numeric_cols <- colnames(data)[!sapply(data, is.numeric)]
    if (length(non_numeric_cols) > 0) {
      warning(paste("The following columns in", data_name, "are not numeric and may cause errors:", 
                    paste(non_numeric_cols, collapse = ", ")))
    }
    
    return(data)
  }
  
  prepare_data <- function(data, data_name) {
    if (is.null(data) || nrow(data) == 0) {
      return(data.frame())
    }
    
    # Handle NA values based on na.action
    if (any(is.na(data))) {
      switch(na.action,
             "allow" = {
               message("Note: NA values found in ", data_name, " but kept (na.action='allow')")
             },
             "omit" = {
               warning("NA values found in ", data_name, ". Removing ", sum(is.na(data)), " NA values.")
               data <- na.omit(data)
               cat(paste("Removed rows with NA values from", data_name, ". New dimensions:", nrow(data), "rows and", ncol(data), "columns.\n"))
             },
             "error" = {
               stop("NA values found in ", data_name, ". Set na.action='omit' to remove or na.action='allow' to keep.")
             }
      )
    }
    
    if (anyDuplicated(rownames(data))) {
      warning(paste("Duplicate row names found in", data_name, "; they have been made unique."))
      rownames(data) <- make.unique(rownames(data))
    }
    
    if (anyDuplicated(colnames(data))) {
      warning(paste("Duplicate column names found in", data_name, "; they have been made unique."))
      colnames(data) <- make.unique(colnames(data))
    }
    
    if (is.null(colnames(data))){
      stop(paste(data_name, "is missing column names."))
    }
    
    cat(paste("Data prepared for", data_name, "with", nrow(data), "rows and", ncol(data), "columns.\n"))
    
    return(data)
  }
  
  clean.data <- ensure_numeric_data(clean.data, "clean.data", convert_factors)
  scale.data <- ensure_numeric_data(scale.data, "scale.data", convert_factors)
  clustered.data <- ensure_numeric_data(clustered.data, "clustered.data", convert_factors)
  
  clean.data <- prepare_data(clean.data, "clean.data")
  scale.data <- prepare_data(scale.data, "scale.data")
  clustered.data <- prepare_data(clustered.data, "clustered.data")
  
  # Match info.data rows with clean.data rows with safety check
  if (nrow(info.data) > 0) {
    common_rows <- intersect(rownames(clean.data), rownames(info.data))
    missing_rows <- setdiff(rownames(clean.data), rownames(info.data))
    
    if (length(missing_rows) > 0) {
      warning("The following rows in clean.data are not found in info.data (", 
              length(missing_rows), " rows): ", 
              paste(head(missing_rows, 5), collapse = ", "),
              ifelse(length(missing_rows) > 5, " ...", ""))
    }
    
    if (length(common_rows) == 0) {
      warning("No matching rows found between clean.data and info.data. Creating empty info.data.")
      matched_info <- data.frame(row.names = rownames(clean.data))
    } else {
      if (length(common_rows) < nrow(clean.data)) {
        warning("Only ", length(common_rows), " of ", nrow(clean.data), 
                " clean.data rows found in info.data. Subsetting clean.data to match.")
        clean.data <- clean.data[common_rows, , drop = FALSE]
      }
      matched_info <- info.data[rownames(clean.data), , drop = FALSE]
    }
  } else {
    matched_info <- data.frame(row.names = rownames(clean.data))
  }
  
  if (is.null(clean.data)) {
    stop("The 'clean.data' must be provided.")
  }
  
  # Use matched_info (the matched version) instead of original info.data
  if (nrow(matched_info) == 0) {
    matched_info <- data.frame(row.names = row.names(clean.data))
    cat("info.data was created from clean.data with", nrow(matched_info), "rows.\n")
  } else if (nrow(matched_info) > 0 && nrow(clean.data) > 0) {
    if (!identical(rownames(clean.data), rownames(matched_info))) {
      warning("Row names of 'clean.data' and 'info.data' are not identical. They will be unified.")
      rownames(matched_info) <- rownames(clean.data)
      cat("Row names of info.data have been unified with clean.data.\n")
    }
  }
  
  Subtyping <- new(
    Class = 'Subtyping',
    clean.data = clean.data,
    info.data = matched_info,   # store the matched version
    scale.data = scale.data,
    Optimal.cluster = Optimal.cluster,
    cluster.results = cluster.results,
    visualization.results = visualization.results,
    clustered.data = clustered.data,
    evaluation_results = evaluation_results
  )
  
  cat("Subtyping object created successfully.\n")
  
  return(Subtyping)
}


#' Split Subtyping object into training and test sets
#'
#' @param object A Subtyping object.
#' @param p Proportion for training set (default 0.7).
#' @param stratify_by Column name in info.data for stratified split (optional).
#' @param seed Random seed.
#' @return List with $train and $test Subtyping objects.
#' @export
#' @examples
#' set.seed(1)
#' demo_df <- data.frame(
#'   id    = paste0("S", 1:60),
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2)
#' )
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#'
#' split <- SplitSubtypingObject(sub_obj, p = 0.7, stratify_by = "group")
#' nrow(split$train@clean.data)
#' nrow(split$test@clean.data)
SplitSubtypingObject <- function(object, p = 0.7, stratify_by = NULL, seed = 123) {
  
  if (!inherits(object, "Subtyping"))
    stop("'object' must be a Subtyping object.")
  
  set.seed(seed)
  n <- nrow(object@clean.data)
  
  if (!is.null(stratify_by) && nrow(object@info.data) > 0 && stratify_by %in% colnames(object@info.data)) {
    strat_var <- object@info.data[[stratify_by]]
  } else {
    strat_var <- NULL
  }
  
  if (is.null(strat_var)) {
    train_idx <- caret::createDataPartition(seq_len(n), p = p, list = FALSE)[, 1]
  } else {
    train_idx <- caret::createDataPartition(strat_var, p = p, list = FALSE)[, 1]
  }
  test_idx <- setdiff(seq_len(n), train_idx)
  
  train_obj <- CreateSubtypingObject(
    clean.data = object@clean.data[train_idx, , drop = FALSE],
    info.data  = if (nrow(object@info.data) > 0) object@info.data[train_idx, , drop = FALSE] else data.frame()
  )
  test_obj <- CreateSubtypingObject(
    clean.data = object@clean.data[test_idx, , drop = FALSE],
    info.data  = if (nrow(object@info.data) > 0) object@info.data[test_idx, , drop = FALSE] else data.frame()
  )
  
  cat("Split completed: training (", length(train_idx), "), test (", length(test_idx), ")\n")
  return(list(train = train_obj, test = test_obj))
}

