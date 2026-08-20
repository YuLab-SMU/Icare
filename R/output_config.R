#' Global Output Path Configuration for icare Package
#' 
#' This file provides functions to configure and manage output paths
#' for all modules (m0-m4). Users can set a custom output directory
#' that will be used instead of the default code folder paths.
#'
#' @name output_config
#' @keywords internal

# Global environment to store configuration
.icare_config <- new.env(parent = emptyenv())

# Default: NULL until the user explicitly calls set_output_root() or
# use_local_output(). get_output_dir() falls back to tempdir() when this
# is NULL, so loading the package never touches the user's working directory.
.icare_config$output_root <- NULL

#' Set Global Output Root Directory
#'
#' @description
#' Sets the root directory for all module outputs. All results will be saved
#' under this directory in module-specific subfolders.
#'
#' @param path Character string. The root output directory path.
#'   If NULL, resets to default (./icare_output).
#'
#' @return Invisible NULL
#' @export
#'
#' @examples
#' \dontrun{
#' # Set output to a custom local folder
#' set_output_root("~/Documents/icare_results")
#' 
#' # Or use absolute path
#' set_output_root("/Users/username/Desktop/analysis_results")
#' 
#' # Reset to default (./icare_output)
#' set_output_root(NULL)
#' }
set_output_root <- function(path = NULL) {
  if (is.null(path)) {
    .icare_config$output_root <- path.expand("./icare_output")
    message("[*] Output root reset to default: ./icare_output")
  } else {
    # Expand user path (e.g., ~ to home directory)
    path <- path.expand(path)
    
    # Create directory if it doesn't exist
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      message("[*] Created output directory: ", path)
    }
    
    .icare_config$output_root <- path
    message("[*] Output root set to: ", path)
    message("[*] Module outputs will be saved to:")
    message("    - m1 (Statistics):     ", file.path(path, "m1"))
    message("    - m2 (ML Models):      ", file.path(path, "m2"))
    message("    - m3 (Subtyping):      ", file.path(path, "m3"))
    message("    - m4 (Prognosis):      ", file.path(path, "m4"))
    message("    - Figures:             ", file.path(path, "Figures"))
  }
  invisible(NULL)
}

#' Get Global Output Root Directory
#'
#' @return Character string with the current output root (default: ./icare_output)
#' @export
#'
#' @examples
#' \dontrun{
#' get_output_root()
#' }
get_output_root <- function() {
  return(.icare_config$output_root)
}

#' Get Module-Specific Output Directory
#'
#' @description
#' Returns the appropriate output directory for a specific module.
#' If a global output root has been explicitly set via \code{set_output_root()}
#' or \code{use_local_output()}, returns a path under that root. Otherwise,
#' falls back to a subdirectory under \code{tempdir()} for the current R
#' session, so the package never writes to the user's working directory
#' unless the user has explicitly opted in.
#'
#' @param module Character string. Module name: "m1", "m2", "m3", "m4", 
#'   "Figures", "StatObject", "ModelData", or "Subtyping".
#' @param subdir Character string. Optional subdirectory path within the module folder.
#'
#' @return Character string with the full output path
#' @export
#'
#' @examples
#' \dontrun{
#' # Get m1 output directory
#' get_output_dir("m1")
#' 
#' # Get m1 output with subdirectory
#' get_output_dir("m1", "categorical_descriptive")
#' 
#' # Get Figures directory
#' get_output_dir("Figures", "PrognosiX")
#' }
get_output_dir <- function(module = c("m1", "m2", "m3", "m4", "Figures", 
                                      "PrognosiX", "StatObject", "ModelData", "Subtyping"),
                           subdir = NULL) {
  module <- match.arg(module)
  
  root <- .icare_config$output_root
  
  # If the user has never explicitly called set_output_root()/use_local_output(),
  # root will be NULL -- fall back to the R session's temp directory instead of
  # the current working directory, so the package never writes to the user's
  # working directory without explicit consent.
  if (is.null(root)) {
    root <- file.path(tempdir(), "icare_output")
  }
  
  base_path <- file.path(root, module)
  
  # Add subdirectory if specified
  if (!is.null(subdir)) {
    base_path <- file.path(base_path, subdir)
  }
  
  # Create directory on-demand if it doesn't exist
  if (!dir.exists(base_path)) {
    dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  return(base_path)
}

#' Initialize Output Directory Structure
#'
#' @description
#' Creates the complete directory structure for all modules.
#' Useful when setting up a new output location.
#'
#' @param root_path Character string. The root output directory.
#'
#' @return Invisible NULL
#' @export
#'
#' @examples
#' \dontrun{
#' init_output_structure("~/Documents/icare_results")
#' }
init_output_structure <- function(root_path) {
  root_path <- path.expand(root_path)
  
  # Only create the root directory, subdirectories will be created on-demand
  if (!dir.exists(root_path)) {
    dir.create(root_path, recursive = TRUE, showWarnings = FALSE)
    message("[*] Created output root directory: ", root_path)
  }
  
  message("[*] Output directory initialized at: ", root_path)
  message("[*] Subdirectories will be created automatically when needed.")
  invisible(NULL)
}

#' Quick Setup for Local Output
#'
#' @description
#' Convenience function to quickly set output to a local directory
#' outside the code folder. Creates the directory structure automatically.
#'
#' @param path Character string. The local output directory path.
#'   Defaults to "~/icare_output" if not specified.
#'
#' @return Invisible NULL
#' @export
#'
#' @examples
#' \dontrun{
#' # Use default location in home directory
#' use_local_output()
#' 
#' # Use custom location
#' use_local_output("~/Documents/my_analysis")
#' use_local_output("/Volumes/external_drive/results")
#' }
use_local_output <- function(path = "./icare_output") {
  path <- path.expand(path)
  
  # Initialize structure
  init_output_structure(path)
  
  # Set as global root
  set_output_root(path)
  
  invisible(NULL)
}
