library(shiny);
library(SingleCellExperiment);
library(Seurat);

# User Interface for Seurat Workflow ---
Seurat_UI <- function(id) {
    ns <- NS(id);
    fluidPage(
        sidebarPanel(fileInput(inputId = ns("sce_rds_file"), label = "Input SingleCellExperiment RDS file", multiple = FALSE, accept = c(".rds"), buttonLabel = "Browse")),
        mainPanel(textOutput(outputId = ns("sce_summary_output")))
    );
}
# ----

# Server for Seurat Workflow ---
Seurat_Server <- function(input, output, session, x) {
    ns <- session$ns

    output$sce_summary_output <- renderText({
        if (!is.null(input$sce_rds_file)) {
            capture.output(SCEtoSeurat(input$sce_rds_file$datapath))
        }
        })
}
# ----

# Helper/Wrapper Functions ---

#' RDStoSCE
#' Reads RDS file and loads into SingleCellExperiment object
#' @param filePath (path of the RDS file)
#'
#' @return SingleCellExperiment object
#' @export
#'
#' @examples
RDStoSCE <- function(filePath) {
    sce <- readRDS(filePath)
    return(sce)
}

#' SCEtoSeurat
#' Converts a SingleCellExperiment object to Seurat object
#' @param filePath (path of the RDS file)
#' 
#' @return Seurat object
#' @export
#'
#' @examples
SCEtoSeurat <- function(filePath) {
    seuratObject <- CreateSeuratObject(counts = counts(RDStoSCE(filePath)))
    return(seuratObject)
}
# ----