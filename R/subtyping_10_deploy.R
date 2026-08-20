#' @title Subtyping Deployment Dispatcher
#' @description 
#' Standardized interface for deploying subtyping models. It handles 
#' feature alignment, normalization synchronization, and space projection.
#'
#' @param sub_train_obj A trained S4 Subtyping object.
#' @param raw_newdata A data frame of raw clinical values (unnormalized).
#' @param method Stratification logic: "nmf", "lpa", or "kmeans".
#'
#' @return A Subtyping object containing predicted cluster labels in info.data.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load your trained object
#' # sub_train <- readRDS("sub_train.rds")
#' # Prepare raw data
#' # test_data <- data.frame(ID = "P001", wbc = 5.2, crp = 12.1, platelets = 250)
#' # results <- Sub_deploy_dispatcher(sub_train, test_data, method = "nmf")
#' }
Sub_deploy_dispatcher <- function(sub_train_obj, raw_newdata, method = "nmf") {
  
  # Ensure input is a data frame
  raw_newdata <- as.data.frame(raw_newdata)
  
  # Feature alignment: extract required variables from the training object
  train_features <- colnames(sub_train_obj@clean.data)
  
  if (!all(train_features %in% colnames(raw_newdata))) {
    missing_vars <- train_features[!train_features %in% colnames(raw_newdata)]
    stop("Input data is missing required features: ", paste(missing_vars, collapse = ", "))
  }
  
  # Normalization Memory Sync
  norm_params <- Sub_extract_norm_params(sub_train_obj)
  
  # Create a temporary Subtyping object for the new data
  # Keep rownames as Sample IDs
  new_obj <- CreateSubtypingObject(
    clean.data = raw_newdata[, train_features, drop = FALSE],
    info.data  = data.frame(SampleID = rownames(raw_newdata), row.names = rownames(raw_newdata))
  )
  
  # Apply normalization parameters from training set
  new_obj <- Sub_apply_norm_params(new_obj, norm_params = norm_params, verbose = FALSE)
  
  # Execute subtype projection using the existing workhorse function
  res_obj <- Sub_predict_subtypes(
    object       = new_obj, 
    train_object = sub_train_obj, 
    method       = method, 
    verbose      = FALSE
  )
  
  return(res_obj)
}


#' @title Create a Subtyping Deployment Manager
#'
#' @description
#' This function serves as a constructor for the \code{Sub_Manager} class. It encapsulates 
#' a trained \code{Subtyping} object and provides a standardized prediction interface 
#' that ensures new data is normalized and processed using the training set's parameters.
#'
#' @param sub_train_obj A trained S4 object of class \code{Subtyping}. This object 
#' should contain the trained models (NMF, LPA, etc.) and normalization parameters.
#'
#' @return A list of class \code{Sub_Manager} containing:
#' \itemize{
#'   \item \code{trained_obj}: The original trained S4 object.
#'   \item \code{sub_predict}: A scoped function for predicting subtypes on new data frames.
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' library(NMF)
#' set.seed(1)
#' demo_df <- data.frame(
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2),
#'   feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#' )
#' rownames(demo_df) <- paste0("S", 1:60)
#'
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' sub_obj <- Sub_normalize_process(sub_obj, normalize_method = "min_max")
#'
#' sub_obj <- Sub_nmf_estimate(sub_obj, rank_range = 2:3, nrun = 2, seed = 1)
#' sub_obj <- Sub_nmf_best_rank(sub_obj, nrun = 2)
#' sub_obj <- Sub_nmf_assign_subtypes(sub_obj)
#' sub_obj <- Sub_nmf_train_model(
#'   sub_obj, best_k = sub_obj@Optimal.cluster, nrun = 2, model_name = "demo_model",save_dir ="./")
#'
#' # Wrap the trained object in a deployment manager for reuse on new samples.
#' manager <- New_Sub_Manager(sub_obj)
#'
#' # Console-only prediction on a few "new" samples (here, reusing training
#' # rows purely to demonstrate the call -- in practice supply genuinely new data).
#' new_samples <- head(sub_obj@clean.data, 5)
#' pred <- manager$sub_predict(new_samples, method = "nmf")
#' pred@info.data$cluster_nmf
#' }
New_Sub_Manager <- function(sub_train_obj) {
  if (!inherits(sub_train_obj, "Subtyping")) {
    stop("The provided object must be of class 'Subtyping'.")
  }
  
  manager <- list(
    trained_obj = sub_train_obj,
    sub_predict = function(newdata, method = "nmf") {
      # This calls the dispatcher which ensures feature alignment and normalization
      Sub_deploy_dispatcher(
        sub_train_obj = sub_train_obj,
        raw_newdata   = newdata,
        method        = method
      )
    }
  )
  class(manager) <- "Sub_Manager"
  return(manager)
}

#' @title Launch the Clinlabomics Subtyping Terminal
#'
#' @description
#' Starts an interactive, high-contrast Shiny dashboard for clinical laboratory 
#' data stratification. The application supports both manual single-sample 
#' entry and batch CSV uploads, providing real-time subtype prediction and 
#' data export capabilities.
#'
#' @param sub_manager A manager object of class \code{Sub_Manager} created 
#' by \code{\link{New_Sub_Manager}}.
#' @param title A character string. The title to be displayed on the dashboard 
#' header. Default is "Subtyping Terminal".
#' @param var_dict A \code{data.frame} (optional). Contains three columns: 
#' 'Feature', 'Description', and 'Units' for documenting the required input variables.
#' @param project_info A \code{list} (optional). Contains elements like 
#' \code{abstract} and \code{citation} to describe the project background in the UI.
#'
#' @return This function starts a Shiny app instance and does not return a value.
#'
#' @import shiny
#' @import shinydashboard
#' @importFrom DT datatable formatStyle
#' 
#' @details 
#' The app uses a dark theme optimized for clinical environments. Batch 
#' uploads require a CSV file where the first column contains unique Sample IDs.
#'
#' @author Huaichao Luo
#'
#' @seealso \code{\link{New_Sub_Manager}}
#'
#' @examples
#' \dontrun{
#' launch_sub_deploy_app(manager, title = "TEP-LUAD Prediction")
#' }
#'
#' @export
launch_sub_deploy_app <- function(sub_manager, title = "Subtyping Terminal", 
                                  var_dict = NULL, project_info = NULL) {
  # --- Data Extraction ---
  # Safely get training features and their medians for default inputs
  train_data <- sub_manager$trained_obj@clean.data
  required_vars <- colnames(train_data)
  default_vals <- sapply(train_data, median, na.rm = TRUE)
  
  # --- UI Definition ---
  ui <- dashboardPage(
    skin = "black",
    dashboardHeader(title = span(icon("dna"), title)),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Analysis Portal", tabName = "portal", icon = icon("desktop")),
        menuItem("Documentation", tabName = "docs", icon = icon("info-circle"))
      )
    ),
    dashboardBody(
      tags$head(tags$style(HTML("
        .content-wrapper { background-color: #121212 !important; }
        .box { background-color: #1e1e1e !important; color: white !important; border-top: 3px solid #00dfc0 !important; }
        label { color: #00dfc0 !important; font-weight: bold; }
        .form-control { background-color: #333 !important; color: white !important; border: 1px solid #444 !important; }
        .btn-run { background-color: #00dfc0 !important; color: #121212 !important; font-size: 24px !important; font-weight: bold !important; padding: 15px !important; border-radius: 12px !important; margin: 20px 0; border: none; width: 60%; }
        .btn-run:hover { background-color: #00ffa2 !important; box-shadow: 0 0 20px #00dfc0; }
        .result-text { font-size: 50px !important; font-weight: 900 !important; color: #00dfc0 !important; }
        .id-label { font-size: 18px; color: #888; font-family: monospace; }
        
        /* DT Pagination visibility fix */
        .dataTables_wrapper .dataTables_paginate .paginate_button { 
          color: #00dfc0 !important; background: #222 !important; border: 1px solid #444 !important; margin: 2px;
        }
        .dataTables_wrapper .dataTables_paginate .paginate_button.current { 
          background: #00dfc0 !important; color: #121212 !important; 
        }
        .dataTables_wrapper .dataTables_info { color: #666 !important; }
        h4 { color: #00dfc0; font-weight: bold; border-left: 4px solid #00dfc0; padding-left: 12px; }
      "))),
      
      tabItems(
        tabItem(tabName = "portal",
                fluidRow(
                  box(width = 12, title = "Project Overview", status = "primary", solidHeader = TRUE,
                      column(8, h4("Abstract"), p(project_info$abstract %||% "Clinical laboratory-based stratification terminal.")),
                      column(4, h4("Citation"), p(project_info$citation %||% "Luo et al. 2026"))
                  )
                ),
                fluidRow(
                  box(title = "1. Configuration & Input", width = 12,
                      column(3, 
                             selectInput("meth", "Algorithm:", choices = c("NMF"="nmf", "LPA"="lpa", "K-means"="kmeans")),
                             hr(),
                             helpText("Batch: First column must be SampleID.")
                      ),
                      column(9, tabsetPanel(id = "input_mode",
                                            tabPanel("Batch (CSV)", br(), 
                                                     fileInput("up", "Upload Sample Data", accept = ".csv"),
                                                     downloadButton("dl_tpl", "Download Template", class="btn-xs")),
                                            tabPanel("Single Sample", br(), 
                                                     fluidRow(
                                                       column(4, textInput("sid", "ID:", "SAMPLE_001")),
                                                       lapply(required_vars, function(v) {
                                                         column(4, numericInput(paste0("in_", v), v, value = round(default_vals[[v]], 2)))
                                                       })
                                                     ))
                      ))
                  )
                ),
                fluidRow(column(12, align="center", actionButton("go", "EXECUTE PREDICTION", icon=icon("play"), class="btn-run"))),
                fluidRow(
                  box(title = "2. Stratification Results", width = 12,
                      column(4, align="center", uiOutput("main_res_ui")),
                      column(8, 
                             div(style="display: flex; justify-content: space-between; align-items: center;",
                                 h4("ID-Linked Data Table"),
                                 downloadButton("dl_res", "Export CSV", class="btn-success btn-xs")
                             ),
                             DT::dataTableOutput("res_table"))
                  )
                )
        ),
        tabItem(tabName = "docs", box(title = "Glossary", width = 12, DT::dataTableOutput("doc_table")))
      )
    )
  )
  
  # --- Server Logic ---
  server <- function(input, output, session) {
    
    # Process inputs into a data frame with IDs as row names
    input_data <- eventReactive(input$go, {
      if (input$input_mode == "Batch (CSV)") {
        req(input$up)
        # Force the first column to be SampleID/RowNames
        df <- read.csv(input$up$datapath, row.names = 1, check.names = FALSE)
        return(df)
      } else {
        vals <- sapply(required_vars, function(v) input[[paste0("in_", v)]])
        df <- as.data.frame(t(vals))
        rownames(df) <- input$sid
        return(df)
      }
    })
    
    # Core Prediction - Validated against the dispatcher
    results <- reactive({
      req(input_data())
      sub_manager$sub_predict(input_data(), method = input$meth)
    })
    
    # Result Summary Display
    output$main_res_ui <- renderUI({
      # Guard against switching methods before result is updated
      validate(need(results(), "Calculating..."))
      res_df <- results()@info.data
      target_col <- paste0("cluster_", input$meth)
      
      if (!target_col %in% colnames(res_df)) return(div("Processing..."))
      
      if (nrow(res_df) == 1) {
        tagList(
          div(class="id-label", paste("ID:", rownames(res_df)[1])),
          div(class="result-text", res_df[[target_col]][1])
        )
      } else {
        tagList(
          div(class="id-label", paste("Processed:", nrow(res_df), "Samples")),
          div(class="result-text", style="font-size:32px;", "Batch Completed")
        )
      }
    })
    
    # Data Table Rendering
    output$res_table <- DT::renderDataTable({
      validate(need(results(), "No data."))
      res_df <- results()@info.data
      target_col <- paste0("cluster_", input$meth)
      
      if (!target_col %in% colnames(res_df)) return(NULL)
      
      # Assemble the output table with SampleID as the first column
      display_df <- data.frame(
        SampleID = rownames(res_df),
        Subtype  = res_df[[target_col]],
        stringsAsFactors = FALSE
      )
      
      datatable(display_df, rownames = FALSE, options = list(
        dom = 'tp',
        pageLength = 5,
        # Using backticks for reserved keyword 'next'
        language = list(paginate = list(previous = "Prev", `next` = "Next"))
      )) %>%
        formatStyle('Subtype', color = '#00dfc0', fontWeight = 'bold')
    })
    
    # --- Export Results (Fixed) ---
    output$dl_res <- downloadHandler(
      filename = function() { paste0("Prediction_Results_", Sys.Date(), ".csv") },
      content = function(file) {
        # Retrieve the most recent info.data
        req(results())
        final_df <- results()@info.data
        # Ensure SampleID is included in the CSV file
        write.csv(final_df, file, row.names = TRUE)
      }
    )
    
    # Template Download
    output$dl_tpl <- downloadHandler(
      filename = "Input_Template.csv",
      content = function(file) { write.csv(head(train_data, 5), file, row.names = TRUE) }
    )
    
    # Documentation Table
    output$doc_table <- DT::renderDataTable({
      df <- if(is.null(var_dict)) data.frame(Feature=required_vars, Desc="Parameter") else var_dict
      datatable(df, rownames = FALSE, options = list(dom = 't')) %>% 
        formatStyle(columns = colnames(df), color = 'white')
    })
  }
  
  shinyApp(ui, server)
}
