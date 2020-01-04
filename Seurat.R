library(shiny);
library(shinyWidgets);
library(SingleCellExperiment);
library(Seurat);

# User Interface for Seurat Workflow ---
Seurat_UI <- function(id) {
    ns <- NS(id);
    fluidPage(
        #Inputs
        sidebarPanel(
                        #data input ui
                        h4("Data Input"),
                        fileInput(inputId = ns("sce_rds_file"), label = "Input SingleCellExperiment RDS file", multiple = FALSE, accept = c(".rds"), buttonLabel = "Browse"),

                        #visualize data ui
                        h4("Visualize Data"),
                        selectInput(inputId = ns("plot_type"), label = "Select plot: ", c("","Violin Plot", "Feature Scatter Plot")),
                        conditionalPanel(condition = sprintf("input['%s'] == 'Violin Plot'", ns("plot_type")), panel(heading = "Violin Plot", selectInput(inputId = ns("violin_plot_select_feature"), label = "Select features: ", choices = c(), multiple = TRUE))),
                        conditionalPanel(condition = sprintf("input['%s'] == 'Feature Scatter Plot'", ns("plot_type")), panel(heading = "Feature Scatter Plot", selectizeInput(inputId = ns("feature_scatter_plot_select_feature"), label = "Select features: ", choices = c(), multiple = TRUE, options = list(maxItems = 2)))),
                        actionButton(inputId = ns("plot_button"), "Plot"),
                        br(),
                        br(),

                        #normalize data ui
                        h4("Normalize Data"),
                        selectInput(inputId = ns("normalization_method"), label = "Select normalization method: ", choices = c("LogNormalize", "CLR", "RC")),
                        textInput(inputId = ns("scale_factor"), label = "Set scaling factor: ", value = "10000"),
                        actionButton(inputId = ns("normalize_button"), "Normalize"),
                        br(),
                        br(),

                        #find highly variable genes ui
                        h4("Highly Variable Genes"),
                        selectInput(inputId = ns("hvg_method"), label = "Select HVG method: ", choices = c("vst", "mean.var.plot", "dispersion")),
                        textInput(inputId = ns("hvg_no_features"), label = "Select number of features: ", value = "2000"),
                        actionButton(inputId = ns("find_hvg_button"), "Find HVG"),
                        actionButton(inputId = ns("plot_hvg_button"), "Plot HVG"),
                        br(),
                        br(),

                        #download data ui
                        h4("Data Download"),
                        downloadButton(outputId = ns("download_seurat_object"), label = "Download Seurat Object"),
                        downloadButton(outputId = ns("download_sce_object"), label = "Download SCE Object"),

                        ),
        #Outputs
        mainPanel(
                    #textouput to show summary of seurat object in use
                    h4("Summary: "),
                    verbatimTextOutput(outputId = ns("sce_summary_output")),
                    #plotoutput to display different plots in use
                    h4("Plot: "),
                    plotOutput(outputId = ns("plot"))
        )
    );
}
# ----

# Server for Seurat Workflow ---
Seurat_Server <- function(input, output, session, x) {
    ns <- session$ns

    #seuratObect as reactiveVal which can be reassigned/updated as workflow progresses
    seuratObject <- reactiveVal()

    #plotObject as reactiveVal which is used throughout the workflow for the plot
    plotObject <- reactiveVal()

    observeEvent(input$sce_rds_file, {
        seuratObject(SCEtoSeurat(input$sce_rds_file$datapath))
        showNotification("Upload Complete")
    })

    observe({
        if (!is.null(input$sce_rds_file)) {
            updateSelectInput(session = session, inputId = "violin_plot_select_feature", choices = getFeatureNamesSeurat(seuratObject = seuratObject()))
            updateSelectInput(session = session, inputId = "feature_scatter_plot_select_feature", choices = getFeatureNamesSeurat(seuratObject = seuratObject()))
        }
    })

    observeEvent(input$plot_button, {
        if (!is.null(input$sce_rds_file) && input$plot_type == "Violin Plot") {
            plotObject(ViolinPlot(seuratObject = seuratObject(), input$violin_plot_select_feature))
        }
        else if (!is.null(input$sce_rds_file) && input$plot_type == "Feature Scatter Plot") {
            plotObject(FeatureScatterPlot(seuratObject = seuratObject(), input$feature_scatter_plot_select_feature))
        }
    })

    observeEvent(input$plot_hvg_button, {
        plotObject(PlotHVG(seuratObject = seuratObject()))
    })

    observeEvent(input$normalize_button, {
    if (!is.null(input$sce_rds_file)) {
        withProgress(message = "Normalizing", max = 1, value = 1, {
            seuratObject(Normalization(seuratObject = seuratObject(), input$normalization_method, as.numeric(input$scale_factor)))
            })
        showNotification("Normalization Complete")
        }
    })

    observeEvent(input$find_hvg_button, {
    if (!is.null(input$sce_rds_file)) {
        withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
            seuratObject(FindHVG(seuratObject = seuratObject(), input$hvg_method, as.numeric(input$hvg_no_features)))
            print(top10 <- head(VariableFeatures(seuratObject()), 10))
            })
        showNotification("Find HVG Complete")
        }
    })


    output$sce_summary_output <- renderText({
        if (!is.null(input$sce_rds_file)) {
            capture.output(seuratObject())
        }
    })

    output$download_seurat_object <- downloadHandler(
        filename = function() {
            "seuratObject.rds"
        },
        content = function(file) {
            saveRDS(seuratObject(), file)
        }
        )

    output$download_sce_object <- downloadHandler(
        filename = function() {
            "sceObject.rds"
        },
        content = function(file) {
            saveRDS(input$sce_rds_file, file)
        }
        )

    output$plot <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject()
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
#' Converts a SingleCellExperiment object to Seurat object using RDS filepath
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

#' ViolinPlot
#' Generates a violin plot against each selected attribute/feature of the seurat object
#' @param seuratObject
#' @param features (list of features selected from seurat object)
#' @return plot object
#' @export
#'
#' @examples
ViolinPlot <- function(seuratObject, features) {
    VlnPlot(seuratObject, features = features, ncol = length(features))
}

#' getFeatureNamesSeurat
#' Extract feature names from Seurat object available in the meta.data slot
#' @param seuratObject
#'
#' @return character list containing feature names
#' @export
#'
#' @examples
getFeatureNamesSeurat <- function(seuratObject) {
    return(names(seuratObject@meta.data))
}


#' FeatureScatterPlot
#' Generates a scatter plot against 'two' selected attributes/features of the seurat object
#' @param seuratObject
#' @param features (list of features selected from seurat object, max=2)
#'
#' @return plot object
#' @export
#'
#' @examples
FeatureScatterPlot <- function(seuratObject, features) {
    FeatureScatter(seuratObject, feature1 = features[1], feature2 = features[2])
}

#' Normalization
#' Normalize the data using one of the selected methods
#' @param seuratObject
#' @param normalizationMethod
#' @param scaleFactor
#'
#' @return seurat object
#' @export
#'
#' @examples
Normalization <- function(seuratObject, normalizationMethod, scaleFactor) {
    return(NormalizeData(seuratObject, normalization.method = normalizationMethod, scale.factor = scaleFactor))
}


#' FindHVG
#' Find highly variable genes using of the selected methods
#' @param seuratObject
#' @param hvgMethod
#' @param hvgNumber
#'
#' @return seurat object
#' @export
#'
#' @examples
FindHVG <- function(seuratObject, hvgMethod, hvgNumber) {
    return(FindVariableFeatures(seuratObject, selection.method = hvgMethod, nfeatures = hvgNumber))
}


#' PlotHVG
#' Plot the highly variable genes
#' @param seuratObject
#'
#' @return plot object
#' @export
#'
#' @examples
PlotHVG <- function(seuratObject) {
    return(VariableFeaturePlot(seuratObject))
}

# ----