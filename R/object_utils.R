## ============================================================
##  object_utils.R
##  Utility functions for Stat / Train_Model / Subtyping / PrognosiX
##  -- ConvertObject   (any-to-any conversion)
##  -- SubsetObject    (samples + features)
##  -- FilterByMeta    (filter by metadata value / range)
##  -- FilterByFeature (filter by expression threshold)
##  -- SplitByMeta     (split into list by metadata column)
##  -- DownsampleObject
##  -- SelectFeatures / RemoveFeatures / RenameFeatures
##  -- AddMetadata
##  -- MergeObjects    (same-class merge)
##  -- InspectObject   (summary printer)
## ============================================================

## -- internal helpers ---------------------------------------------------------

.check_class <- function(object, allowed) {
  if (!inherits(object, allowed))
    stop("Object must be one of: ", paste(allowed, collapse = ", "),
         ". Got: ", class(object)[1])
}

.subset_df <- function(df, rows = NULL, cols = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!is.null(rows)) df <- df[rownames(df) %in% rows, , drop = FALSE]
  if (!is.null(cols)) df <- df[, colnames(df) %in% cols,  drop = FALSE]
  df
}

.safe_rbind <- function(df1, df2) {
  if (is.null(df1) || nrow(df1) == 0) return(df2)
  if (is.null(df2) || nrow(df2) == 0) return(df1)
  shared <- intersect(colnames(df1), colnames(df2))
  if (length(shared) == 0) {
    warning("No shared columns; returning first data frame only.")
    return(df1)
  }
  if (length(shared) < max(ncol(df1), ncol(df2)))
    message(length(shared), " shared columns kept (of ",
            ncol(df1), " / ", ncol(df2), ").")
  out <- rbind(df1[, shared, drop = FALSE], df2[, shared, drop = FALSE])
  rownames(out) <- make.unique(rownames(out))
  out
}

.object_clean_data <- function(object) {
  if (inherits(object, "Stat"))        return(object@clean.data)
  if (inherits(object, "Subtyping"))   {
    cd <- object@clustered.data
    if (!is.null(cd) && nrow(cd) > 0) return(cd)
    return(object@clean.data)
  }
  if (inherits(object, "PrognosiX"))   return(object@clean.data)
  if (inherits(object, "Train_Model")) return(object@clean.df)
  stop("Unsupported object class: ", class(object)[1])
}

.object_info_data <- function(object) {
  if (inherits(object, "Train_Model")) return(data.frame())
  tryCatch(slot(object, "info.data"), error = function(e) data.frame())
}

## -- 1. ConvertObject ---------------------------------------------------------
#' Convert any package object to another class
#'
#' A single entry-point for all cross-class conversions.
#'
#' Supported conversions
#' \tabular{ll}{
#'   \strong{From}  \tab \strong{To} \cr
#'   Stat           \tab Train_Model, Subtyping, PrognosiX \cr
#'   Subtyping      \tab Stat, PrognosiX, Train_Model \cr
#'   PrognosiX      \tab Stat, Subtyping \cr
#'   Train_Model    \tab Stat, Subtyping \cr
#' }
#'
#' @param object     Source object (any of the four classes).
#' @param to         Target class name as a string: \code{"Stat"},
#'   \code{"Train_Model"}, \code{"Subtyping"}, or \code{"PrognosiX"}.
#' @param time_col   Survival time column (required when \code{to = "PrognosiX"}
#'   and the column is not yet named \code{"time"}).
#' @param status_col Event status column (required when \code{to = "PrognosiX"}
#'   and the column is not yet named \code{"status"}).
#' @param group_col  Label column (used when converting to \code{Train_Model}
#'   from a class that has no \code{group_col} slot).
#' @returns An object of the requested class.
#' @export
#' @examples
#' \dontrun{
#' Subtem=ConvertObject(stat_obj_test,to='Subtyping')
#' Protem=ConvertObject(stat_obj_test,to='PrognosiX')
#' Traintem=ConvertObject(stat_obj_test,to='Train_Model')
#' set.seed(1)
#' demo_df <- data.frame(
#'   id    = paste0("S", 1:60),
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2),
#'   feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#' )
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' class(sub_obj)
#' }
ConvertObject <- function(object,
                          to,
                          time_col   = "time",
                          status_col = "status",
                          group_col  = NULL) {
  
  allowed_from <- c("Stat", "Subtyping", "PrognosiX", "Train_Model")
  allowed_to   <- c("Stat", "Subtyping", "PrognosiX", "Train_Model")
  .check_class(object, allowed_from)
  
  if (!to %in% allowed_to)
    stop("'to' must be one of: ", paste(allowed_to, collapse = ", "))
  
  from <- class(object)[1]
  
  if (from == to) {
    message("Source and target class are identical; returning object unchanged.")
    return(object)
  }
  
  # -- helper: build a minimal Stat from any object ---------------------
  # Track conversion chain to prevent infinite recursion
  .to_stat <- function(obj, grp = group_col, visited = NULL) {
    # Initialize visited set if not provided
    if (is.null(visited)) {
      visited <- character()
    }
    
    # Check for circular reference
    obj_id <- paste0(class(obj)[1], "_", digest::digest(obj))
    if (obj_id %in% visited) {
      stop("Circular reference detected in object conversion. Conversion chain: ", 
           paste(visited, collapse = " -> "), " -> ", obj_id)
    }
    visited <- c(visited, obj_id)
    
    cd   <- .object_clean_data(obj)
    info <- .object_info_data(obj)
    if (is.null(grp) && .hasSlot(obj, "group_col")) grp <- obj@group_col
    CreateStatObject(clean.data = cd,
                     info.data  = if (nrow(info) > 0) info else data.frame(),
                     group_col  = if (!is.null(grp)) grp else "group")
  }
  
  # -- dispatch table ----------------------------------------------------
  result <- switch(
    
    paste(from, to, sep = "_to_"),
    
    ## -- Stat -> * -----------------------------------------------------
    Stat_to_Train_Model = {
      grp <- if (!is.null(group_col)) group_col else object@group_col
      CreateModelObject(object = object, group_col = grp)
    },
    
    Stat_to_Subtyping = {
      # Determine the group column name from the Stat object
      grp_col <- if (!is.null(group_col)) group_col else object@group_col
      
      # -- Extract group labels BEFORE CreateSubtypingObject -----------
      # CreateSubtypingObject calls ensure_numeric_data() internally, which
      # converts character/factor group labels (e.g. "A", "B") into integer
      # codes (0, 1, 2 ...).  We must capture the original labels here, from
      # the unmodified Stat clean.data, so that info.data receives the human-
      # readable group names rather than numeric codes.
      grp_values <- if (!is.null(grp_col) && grp_col %in% colnames(object@clean.data))
        object@clean.data[[grp_col]]
      else
        NULL
      
      # Build the Subtyping object from the Stat object
      orig_colnames <- colnames(object@clean.data)
      sub_obj <- CreateSubtypingObject(object = object)
      
      # Restore column names preserved before CreateSubtypingObject may reformat them
      colnames(sub_obj@clean.data) <- orig_colnames
      if (ncol(sub_obj@scale.data) > 0) {
        colnames(sub_obj@scale.data) <- orig_colnames
      }
      
      # -- Move group_col from clean.data -> info.data ------------------
      # If the group column is present in clean.data, migrate it to info.data
      # and strip it from the numeric feature matrices (clean.data, scale.data).
      if (!is.null(grp_values) && grp_col %in% colnames(sub_obj@clean.data)) {
        
        # Initialise info.data if it is empty, preserving sample row names
        if (is.null(sub_obj@info.data) || nrow(sub_obj@info.data) == 0) {
          sub_obj@info.data <- data.frame(row.names = rownames(sub_obj@clean.data))
        }
        
        # Write original (non-numeric-coded) group labels into info.data
        sub_obj@info.data[[grp_col]] <- grp_values
        
        # Remove group column from clean.data
        sub_obj@clean.data <- sub_obj@clean.data[
          , colnames(sub_obj@clean.data) != grp_col, drop = FALSE]
        
        # Remove group column from scale.data if present
        if (ncol(sub_obj@scale.data) > 0 &&
            grp_col %in% colnames(sub_obj@scale.data)) {
          sub_obj@scale.data <- sub_obj@scale.data[
            , colnames(sub_obj@scale.data) != grp_col, drop = FALSE]
        }
        
        message("Stat \u2192 Subtyping: column '", grp_col,
                "' moved from clean.data to info.data.")
      }
      
      sub_obj
    },
    
    Stat_to_PrognosiX = {
      CreatePrognosiXObject(object     = object,
                            time_col   = time_col,
                            status_col = status_col)
    },
    
    ## -- Subtyping -> * ------------------------------------------------
    Subtyping_to_Stat = {
      .to_stat(object, group_col)
    },
    
    Subtyping_to_PrognosiX = {
      CreatePrognosiXObject(object     = object,
                            time_col   = time_col,
                            status_col = status_col)
    },
    
    Subtyping_to_Train_Model = {
      stat_tmp <- .to_stat(object, group_col)
      grp <- if (!is.null(group_col)) group_col else "group"
      CreateModelObject(object = stat_tmp, group_col = grp)
    },
    
    ## -- PrognosiX -> * ------------------------------------------------
    PrognosiX_to_Stat = {
      .to_stat(object, group_col)
    },
    
    PrognosiX_to_Subtyping = {
      CreateSubtypingObject(object = object)
    },
    
    PrognosiX_to_Train_Model = {
      stat_tmp <- .to_stat(object, group_col)
      grp <- if (!is.null(group_col)) group_col else "group"
      CreateModelObject(object = stat_tmp, group_col = grp)
    },
    
    ## -- Train_Model -> * ----------------------------------------------
    Train_Model_to_Stat = {
      grp <- if (!is.null(group_col)) group_col else object@group_col
      CreateStatObject(clean.data = object@clean.df, group_col = grp)
    },
    
    Train_Model_to_Subtyping = {
      stat_tmp <- CreateStatObject(clean.data = object@clean.df,
                                   group_col  = object@group_col)
      CreateSubtypingObject(object = stat_tmp)
    },
    
    Train_Model_to_PrognosiX = {
      stat_tmp <- CreateStatObject(clean.data = object@clean.df,
                                   group_col  = object@group_col)
      CreatePrognosiXObject(object     = stat_tmp,
                            time_col   = time_col,
                            status_col = status_col)
    },
    
    ## -- unsupported ---------------------------------------------------
    stop("Conversion from '", from, "' to '", to, "' is not supported.")
  )
  
  cat(from, "->", to, "conversion complete.\n")
  return(result)
}

## -- 2. SubsetObject ----------------------------------------------------------

#' Subset a package object by sample and / or feature names
#'
#' A unified wrapper that dispatches to the correct slot logic for each class.
#'
#' @param object   Any of \code{Stat}, \code{Train_Model}, \code{Subtyping},
#'   \code{PrognosiX}.
#' @param samples  Character vector of sample (row) names to keep. \code{NULL}
#'   keeps all samples.
#' @param features Character vector of feature (column) names to keep.
#'   \code{NULL} keeps all features.
#' @returns A subsetted object of the same class.
#' @export
#' @examples
#' \dontrun{
#' sub_obj <- SubsetObject(stat_obj_test, samples = rownames(stat_obj_test@clean.data)[1:5])
#' }
SubsetObject <- function(object, samples = NULL, features = NULL) {
  
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  if (inherits(object, "Stat")) {
    object@raw.data   <- .subset_df(object@raw.data,   samples, features)
    object@clean.data <- .subset_df(object@clean.data, samples, features)
    object@info.data  <- .subset_df(object@info.data,  samples, NULL)
    object@scale.data <- .subset_df(object@scale.data, samples, features)
    if (length(object@meta.featurename) > 0 && !is.null(features))
      object@meta.featurename <- intersect(object@meta.featurename, features)
    cat("Stat subsetted:", nrow(object@clean.data), "x",
        ncol(object@clean.data), "\n")
  }
  
  else if (inherits(object, "Subtyping")) {
    object@clean.data     <- .subset_df(object@clean.data,     samples, features)
    object@scale.data     <- .subset_df(object@scale.data,     samples, features)
    object@info.data      <- .subset_df(object@info.data,      samples, NULL)
    object@clustered.data <- .subset_df(object@clustered.data, samples, features)
    cat("Subtyping subsetted:", nrow(object@clean.data), "x",
        ncol(object@clean.data), "\n")
  }
  
  else if (inherits(object, "PrognosiX")) {
    object@clean.data    <- .subset_df(object@clean.data,    samples, features)
    object@info.data     <- .subset_df(object@info.data,     samples, NULL)
    object@survival.data <- .subset_df(object@survival.data, samples, NULL)
    object@sub.data      <- .subset_df(object@sub.data,      samples, NULL)
    cat("PrognosiX subsetted:", nrow(object@clean.data), "x",
        ncol(object@clean.data), "\n")
  }
  
  else if (inherits(object, "Train_Model")) {
    object@data.df  <- .subset_df(object@data.df,  samples, features)
    object@clean.df <- .subset_df(object@clean.df, samples, features)
    # invalidate downstream slots
    object@split.data        <- list()
    object@split.scale.data  <- list()
    object@train.models      <- list()
    object@all.results       <- list()
    object@best.model.result <- list()
    cat("Train_Model subsetted:", nrow(object@clean.df), "x",
        ncol(object@clean.df), "\n")
    message("Note: split / model slots reset after subsetting.")
  }
  
  return(object)
}


## -- 3. FilterByMeta ----------------------------------------------------------

#' Filter object samples by a metadata column
#'
#' Keeps only samples whose value in a specified \code{info.data} column
#' matches a set of allowed values, or falls within a numeric range.
#'
#' @param object  A \code{Stat}, \code{Subtyping}, or \code{PrognosiX} object.
#' @param col     Column name in \code{info.data} (or \code{clean.data} for
#'   \code{Train_Model}).
#' @param values  Character / factor vector of allowed values. Ignored when
#'   \code{range} is supplied.
#' @param range   Numeric vector of length 2 \code{c(min, max)} for numeric
#'   columns (inclusive). \code{values} is ignored when \code{range} is set.
#' @param invert  Logical. When \code{TRUE} the matching samples are
#'   \emph{excluded} rather than kept. Default \code{FALSE}.
#' @returns An object of the same class containing only the filtered samples.
#' @export
#' @examples
#' \dontrun{
#' sub_obj <- FilterByMeta(stat_obj_test, col = "SWAB", values = "0")
#' }
FilterByMeta <- function(object, col, values = NULL, range = NULL,
                         invert = FALSE) {
  
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  # locate the metadata data frame
  meta_df <- tryCatch(.object_info_data(object), error = function(e) data.frame())
  if (nrow(meta_df) == 0 || !col %in% colnames(meta_df)) {
    # fall back to clean data for Train_Model
    if (inherits(object, "Train_Model")) {
      meta_df <- object@clean.df
    } else {
      stop("Column '", col, "' not found in info.data.")
    }
  }
  
  if (!is.null(range)) {
    if (length(range) != 2 || !is.numeric(range))
      stop("'range' must be a numeric vector of length 2: c(min, max).")
    vec <- suppressWarnings(as.numeric(meta_df[[col]]))
    keep_flag <- !is.na(vec) & vec >= range[1] & vec <= range[2]
  } else {
    if (is.null(values))
      stop("Provide either 'values' or 'range'.")
    keep_flag <- meta_df[[col]] %in% values
  }
  
  if (invert) keep_flag <- !keep_flag
  kept_samples <- rownames(meta_df)[keep_flag]
  
  if (length(kept_samples) == 0)
    warning("No samples matched the filter; returning original object.")
  
  cat("FilterByMeta: keeping", sum(keep_flag), "of",
      nrow(meta_df), "samples.\n")
  SubsetObject(object, samples = kept_samples)
}


## -- 4. FilterByFeature -------------------------------------------------------

#' Filter samples by a feature value threshold
#'
#' Keeps samples where the value of a specific feature in \code{clean.data}
#' is above or below a given threshold.
#'
#' @param object    Any of the four object classes.
#' @param feature   Feature (column) name in \code{clean.data}.
#' @param threshold Numeric cut-off value.
#' @param direction \code{"above"} (default) keeps samples >= threshold;
#'   \code{"below"} keeps samples <= threshold; \code{"equal"} keeps exact
#'   matches.
#' @param invert    Logical. Inverts the selection when \code{TRUE}.
#' @returns An object of the same class.
#' @export
#' @examples
#' \dontrun{
#' stat_obj <- FilterByFeature(stat_obj_test, feature = "AGE", threshold = 10, direction = "above")
#' }
FilterByFeature <- function(object, feature, threshold,
                            direction = c("above", "below", "equal"),
                            invert = FALSE) {
  
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  direction <- match.arg(direction)
  
  cd <- .object_clean_data(object)
  if (!feature %in% colnames(cd))
    stop("Feature '", feature, "' not found in clean data.")
  
  vec <- suppressWarnings(as.numeric(cd[[feature]]))
  keep_flag <- switch(direction,
                      above = !is.na(vec) & vec >= threshold,
                      below = !is.na(vec) & vec <= threshold,
                      equal = !is.na(vec) & vec == threshold
  )
  if (invert) keep_flag <- !keep_flag
  kept <- rownames(cd)[keep_flag]
  
  cat("FilterByFeature: keeping", sum(keep_flag), "of",
      nrow(cd), "samples.\n")
  SubsetObject(object, samples = kept)
}


## -- 5. SplitByMeta -----------------------------------------------------------

#' Split an object into a list by a metadata column
#'
#' Splits the object into a named list of objects, one per unique value in the
#' specified \code{info.data} column.  Equivalent to Seurat's
#' \code{SplitObject}.
#'
#' @param object  A \code{Stat}, \code{Subtyping}, or \code{PrognosiX} object.
#' @param col     Column name in \code{info.data}.
#' @returns A named list of objects, one element per unique level of \code{col}.
#' @export
#' @examples
#' \dontrun{
#' split_list <- SplitByMeta(stat_obj_test, col = "SWAB")
#' }
SplitByMeta <- function(object, col) {
  
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  meta_df <- tryCatch(.object_info_data(object), error = function(e) data.frame())
  if (nrow(meta_df) == 0 || !col %in% colnames(meta_df))
    stop("Column '", col, "' not found in info.data.")
  
  levels_vec <- unique(meta_df[[col]])
  result <- lapply(levels_vec, function(lv) {
    ids <- rownames(meta_df)[meta_df[[col]] == lv]
    SubsetObject(object, samples = ids)
  })
  names(result) <- as.character(levels_vec)
  
  cat("SplitByMeta: split into", length(result),
      "objects by column '", col, "'.\n")
  return(result)
}


## -- 6. DownsampleObject ------------------------------------------------------

#' Downsample an object to a fixed number of samples
#'
#' Randomly retains \code{n} samples, optionally balanced per group.
#'
#' @param object     Any of the four object classes.
#' @param n          Number of samples to keep (total, or per group when
#'   \code{per_group = TRUE}).
#' @param per_group  Logical. When \code{TRUE}, \code{n} samples are drawn
#'   \emph{per unique level} of \code{group_col}. Default \code{FALSE}.
#' @param group_col  Column to use for per-group sampling. When \code{NULL}
#'   the slot \code{group_col} from the object is used.
#' @param seed       Random seed for reproducibility. Default \code{42}.
#' @returns A downsampled object of the same class.
#' @export
#' @examples
#' \dontrun{
#' down <- DownsampleObject(stat_obj_test, n = 10)
#' }
DownsampleObject <- function(object, n, per_group = FALSE,
                             group_col = NULL, seed = 42) {
  
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  set.seed(seed)
  
  all_samples <- rownames(.object_clean_data(object))
  
  if (!per_group) {
    n <- min(n, length(all_samples))
    kept <- sample(all_samples, n)
  } else {
    # determine grouping vector
    grp_col <- group_col
    if (is.null(grp_col) && .hasSlot(object, "group_col"))
      grp_col <- object@group_col
    if (is.null(grp_col))
      stop("Provide 'group_col' for per-group downsampling.")
    
    meta_df <- tryCatch(.object_info_data(object), error = function(e) data.frame())
    cd      <- .object_clean_data(object)
    src_df  <- if (nrow(meta_df) > 0 && grp_col %in% colnames(meta_df)) meta_df else cd
    if (!grp_col %in% colnames(src_df))
      stop("Column '", grp_col, "' not found.")
    
    groups <- src_df[[grp_col]]
    levels_vec <- unique(groups)
    kept <- unlist(lapply(levels_vec, function(lv) {
      ids <- all_samples[groups == lv]
      sample(ids, min(n, length(ids)))
    }))
  }
  
  cat("DownsampleObject: kept", length(kept), "of",
      length(all_samples), "samples.\n")
  SubsetObject(object, samples = kept)
}


## -- 7. SelectFeatures / RemoveFeatures / RenameFeatures ----------------------

#' Keep a set of features in an object
#'
#' Keeps only the specified features (columns) in all relevant data slots.
#' Supports direct name vectors or regex patterns.
#'
#' @param object   Any of the four object classes.
#' @param features Character vector of feature names to keep. Ignored when
#'   \code{pattern} is given.
#' @param pattern  A regular expression; all matching column names are kept.
#' @returns Object with only the selected features.
#' @export
#' @examples
#' \dontrun{
#' selected <- SelectFeatures(stat_obj_test, features = c("AGE"))
#' }
SelectFeatures <- function(object, features = NULL, pattern = NULL) {
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  all_cols <- colnames(.object_clean_data(object))
  if (!is.null(pattern)) {
    features <- grep(pattern, all_cols, value = TRUE)
    if (length(features) == 0)
      warning("Pattern '", pattern, "' matched no features.")
  }
  if (is.null(features)) stop("Provide 'features' or 'pattern'.")
  
  cat("SelectFeatures: keeping", length(features), "of",
      length(all_cols), "features.\n")
  SubsetObject(object, features = features)
}


#' Remove specific features from an object
#'
#' Drops the specified features (columns) from all relevant data slots.
#' Supports direct name vectors or regex patterns.
#'
#' @param object   Any of the four object classes.
#' @param features Character vector of feature names to remove. Ignored when
#'   \code{pattern} is given.
#' @param pattern  A regular expression; all matching column names are removed.
#' @returns Object without the removed features.
#' @export
#' @examples
#' \dontrun{
#' selected <- RemoveFeatures(stat_obj_test, features = c("AGE"))
#' }
RemoveFeatures <- function(object, features = NULL, pattern = NULL) {
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  all_cols <- colnames(.object_clean_data(object))
  if (!is.null(pattern)) {
    features <- grep(pattern, all_cols, value = TRUE)
    if (length(features) == 0) {
      warning("Pattern '", pattern, "' matched no features; returning unchanged.")
      return(object)
    }
  }
  if (is.null(features)) stop("Provide 'features' or 'pattern'.")
  
  keep <- setdiff(all_cols, features)
  cat("RemoveFeatures: removing", length(features), "features;",
      length(keep), "remain.\n")
  SubsetObject(object, features = keep)
}


#' Rename features in an object
#'
#' Renames columns in all data slots simultaneously to keep the object
#' consistent.
#'
#' @param object Any of the four object classes.
#' @param old    Character vector of current feature names.
#' @param new    Character vector of replacement names (same length as
#'   \code{old}).
#' @returns Object with renamed features.
#' @export
#' @examples
#' \dontrun{
#' selected <- RenameFeatures(stat_obj_test, old = "AGE", new = "age")
#' }
RenameFeatures <- function(object, old, new) {
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  if (length(old) != length(new))
    stop("'old' and 'new' must be the same length.")
  
  rename_cols <- function(df, old, new) {
    if (is.null(df) || nrow(df) == 0) return(df)
    idx <- match(old, colnames(df))
    found <- !is.na(idx)
    colnames(df)[idx[found]] <- new[found]
    df
  }
  
  if (inherits(object, "Stat")) {
    object@raw.data   <- rename_cols(object@raw.data,   old, new)
    object@clean.data <- rename_cols(object@clean.data, old, new)
    object@scale.data <- rename_cols(object@scale.data, old, new)
    object@meta.featurename <- colnames(object@clean.data)
  } else if (inherits(object, "Subtyping")) {
    object@clean.data     <- rename_cols(object@clean.data,     old, new)
    object@scale.data     <- rename_cols(object@scale.data,     old, new)
    object@clustered.data <- rename_cols(object@clustered.data, old, new)
  } else if (inherits(object, "PrognosiX")) {
    object@clean.data    <- rename_cols(object@clean.data,    old, new)
    object@survival.data <- rename_cols(object@survival.data, old, new)
  } else if (inherits(object, "Train_Model")) {
    object@data.df  <- rename_cols(object@data.df,  old, new)
    object@clean.df <- rename_cols(object@clean.df, old, new)
  }
  
  cat("RenameFeatures: renamed", sum(old %in% colnames(.object_clean_data(object)) |
                                       new %in% colnames(.object_clean_data(object))),
      "features.\n")
  return(object)
}


## -- 8. AddMetadata -----------------------------------------------------------

#' Add or overwrite a column in info.data
#'
#' Appends (or replaces) a metadata column in the \code{info.data} slot.
#' Checks that \code{values} length matches the number of samples.
#'
#' @param object Any of the four object classes.
#' @param col    Name of the new (or existing) metadata column.
#' @param values Vector of values; must match \code{nrow(info.data)} or the
#'   number of samples in the object.
#' @returns Object with the updated metadata.
#' @export
#' @examples
#' \dontrun{
#'stat_obj <- AddMetadata(stat_obj_test, col = "new_meta", values = 
#'sample(1:2, nrow(stat_obj_test@clean.data), replace=TRUE))
#'}
AddMetadata <- function(object, col, values) {
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  cd <- .object_clean_data(object)
  n  <- nrow(cd)
  if (length(values) != n)
    stop("'values' length (", length(values), ") must equal number of samples (",
         n, ").")
  
  if (inherits(object, "Train_Model")) {
    object@clean.df[[col]] <- values
    cat("AddMetadata: column '", col, "' added to clean.df.\n")
  } else {
    if (nrow(object@info.data) == 0)
      object@info.data <- data.frame(row.names = rownames(cd))
    object@info.data[[col]] <- values
    cat("AddMetadata: column '", col, "' added to info.data.\n")
  }
  return(object)
}


## -- 9. MergeObjects ----------------------------------------------------------

#' Merge two or more objects of the same class
#'
#' Row-binds clean data, info data and other shared slots. Downstream
#' results (models, clustering) are reset for objects where these would be
#' invalidated.
#'
#' @param x One object (any of the four classes) or a list of objects of the
#'   same class.
#' @param y A single object or a list of objects of the same class as \code{x}.
#'   Ignored when \code{x} is already a list.
#' @returns A merged object of the same class.
#' @export
#' @examples
#' \dontrun{
#'stat_obj <- MergeObjects(stat_obj_test, stat_obj_test)
#'}
MergeObjects <- function(x, y = NULL) {
  
  # normalise to a flat list
  if (is.list(x) && !isS4(x)) {
    obj_list <- x
  } else {
    obj_list <- c(list(x), if (is.list(y) && !isS4(y)) y else list(y))
    obj_list <- obj_list[!sapply(obj_list, is.null)]
  }
  if (length(obj_list) < 2)
    stop("Provide at least two objects to merge.")
  
  from <- class(obj_list[[1]])[1]
  if (!all(sapply(obj_list, function(o) inherits(o, from))))
    stop("All objects must be of the same class.")
  .check_class(obj_list[[1]],
               c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  merged_clean <- Reduce(.safe_rbind, lapply(obj_list, .object_clean_data))
  merged_info  <- Reduce(.safe_rbind, lapply(obj_list, .object_info_data))
  if (nrow(merged_info) > 0)
    rownames(merged_info) <- make.unique(rownames(merged_info))
  
  result <- switch(from,
                   
                   Stat = {
                     merged_raw <- Reduce(.safe_rbind,
                                          lapply(obj_list, function(o) o@raw.data))
                     CreateStatObject(
                       raw.data   = merged_raw,
                       clean.data = merged_clean,
                       info.data  = merged_info,
                       group_col  = obj_list[[1]]@group_col
                     )
                   },
                   
                   Subtyping = {
                     obj <- CreateSubtypingObject(
                       clean.data = merged_clean,
                       info.data  = merged_info
                     )
                     message("Subtyping merge: cluster.results and visualization.results reset.")
                     obj
                   },
                   
                   PrognosiX = {
                     obj <- CreatePrognosiXObject(
                       clean.data = merged_clean,
                       info.data  = merged_info,
                       time_col   = obj_list[[1]]@time_col,
                       status_col = obj_list[[1]]@status_col
                     )
                     message("PrognosiX merge: model slots reset; re-run survival modelling.")
                     obj
                   },
                   
                   Train_Model = {
                     obj <- CreateModelObject(data      = merged_clean,
                                              group_col = obj_list[[1]]@group_col)
                     message("Train_Model merge: split / model slots reset.")
                     obj
                   }
  )
  
  cat("MergeObjects:", length(obj_list), from, "objects merged ->",
      nrow(.object_clean_data(result)), "samples total.\n")
  return(result)
}


## -- 10. InspectObject --------------------------------------------------------

#' Print a structured summary of any package object
#'
#' Displays slot-level dimensions, key settings and the fill state of
#' list/result slots so you can quickly see how far processing has progressed.
#'
#' @param object Any of the four object classes.
#' @returns \code{invisible(object)} (called for side-effects).
#' @export
#' @examples
#' \dontrun{
#'stat_obj <- InspectObject(stat_obj_test)
#'}
InspectObject <- function(object) {
  .check_class(object, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  
  .dim_str <- function(df) {
    if (is.null(df) || nrow(df) == 0) return("empty")
    paste0(nrow(df), " x ", ncol(df))
  }
  .list_str <- function(lst) {
    if (length(lst) == 0) return("(empty)")
    paste0("(", length(lst), " element(s): ", paste(names(lst), collapse = ", "), ")")
  }
  
  cat("==========================================\n")
  cat("  Object class :", class(object)[1], "\n")
  cat("==========================================\n")
  
  if (inherits(object, "Stat")) {
    cat("  raw.data           :", .dim_str(object@raw.data), "\n")
    cat("  clean.data         :", .dim_str(object@clean.data), "\n")
    cat("  info.data          :", .dim_str(object@info.data), "\n")
    cat("  scale.data         :", .dim_str(object@scale.data), "\n")
    cat("  group_col          :", object@group_col, "\n")
    cat("  n features         :", length(object@meta.featurename), "\n")
    cat("  variable.types     :", .list_str(object@variable.types), "\n")
    cat("  compute.descriptive:", .list_str(object@compute.descriptive), "\n")
    cat("  corr.result        :", .list_str(object@corr.result), "\n")
    cat("  process.info       :", .list_str(object@process.info), "\n")
  }
  
  else if (inherits(object, "Train_Model")) {
    cat("  data.df            :", .dim_str(object@data.df), "\n")
    cat("  clean.df           :", .dim_str(object@clean.df), "\n")
    cat("  group_col          :", as.character(object@group_col), "\n")
    cat("  split.data         :", .list_str(object@split.data), "\n")
    cat("  feature.selection  :", .list_str(object@feature.selection), "\n")
    cat("  train.models       :", .list_str(object@train.models), "\n")
    cat("  all.results        :", .list_str(object@all.results), "\n")
    cat("  best.model.result  :", .list_str(object@best.model.result), "\n")
  }
  
  else if (inherits(object, "Subtyping")) {
    cat("  clean.data         :", .dim_str(object@clean.data), "\n")
    cat("  info.data          :", .dim_str(object@info.data), "\n")
    cat("  scale.data         :", .dim_str(object@scale.data), "\n")
    cat("  clustered.data     :", .dim_str(object@clustered.data), "\n")
    cat("  Optimal.cluster    :", as.character(object@Optimal.cluster), "\n")
    cat("  cluster.results    :",
        if (is.null(object@cluster.results)) "(NULL)" else "(set)", "\n")
    cat("  visualization.res  :",
        if (is.null(object@visualization.results)) "(NULL)" else "(set)", "\n")
    cat("  evaluation_results :", .list_str(object@evaluation_results), "\n")
  }
  
  else if (inherits(object, "PrognosiX")) {
    cat("  clean.data         :", .dim_str(object@clean.data), "\n")
    cat("  info.data          :", .dim_str(object@info.data), "\n")
    cat("  survival.data      :", .dim_str(object@survival.data), "\n")
    cat("  sub.data           :", .dim_str(object@sub.data), "\n")
    cat("  time_col           :", as.character(object@time_col), "\n")
    cat("  status_col         :", as.character(object@status_col), "\n")
    cat("  baseline.table     :",
        if (is.null(object@baseline.table)) "(NULL)" else "(set)", "\n")
    cat("  univariate.analysis:", .list_str(object@univariate.analysis), "\n")
    cat("  split.data         :", .list_str(object@split.data), "\n")
    cat("  feature.result     :", .list_str(object@feature.result), "\n")
    cat("  survival.model     :",
        if (is.null(object@survival.model)) "(NULL)" else "(set)", "\n")
    cat("  best.model         :", .list_str(object@best.model), "\n")
    cat("  subgroup.risk      :", .list_str(object@subgroup.risk), "\n")
  }
  
  cat("==========================================\n")
  invisible(object)
}
