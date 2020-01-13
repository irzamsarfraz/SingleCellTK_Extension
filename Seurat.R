library(shiny);
library(shinyWidgets);
library(SingleCellExperiment);
library(Seurat);

# User Interface for Seurat Workflow ---
Seurat_UI <- function(id) {
    ns <- NS(id);
    fluidPage(
        #Inputs
        sidebarPanel(style = "overflow-y:scroll; max-height: 500px; position:relative;",
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

                        #scale data ui
                        h4("Scale Data"),
                        actionButton(inputId = ns("scale_button"), "Scale"),
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

                        #principal component analysis
                        h4("Principal Component Analysis"),
                        actionButton(inputId = ns("run_pca_button"), "Run PCA"),
                        actionButton(inputId = ns("plot_pca_button"), "Plot PCA"),
                        br(),
                        br(),

                        #find clusters
                        h4("Clustering"),
                        actionButton(inputId = ns("find_clusters_button"), "Find Clusters"),
                        br(),
                        br(),

                        #tSNE
                        h4("tSNE"),
                        actionButton(inputId = ns("run_tsne_button"), "Run tSNE"),
                        actionButton(inputId = ns("plot_tsne_button"), "Plot tSNE"),
                        br(),
                        br(),

                        #download data ui
                        h4("Data Download"),
                        downloadButton(outputId = ns("download_seurat_object"), label = "Download Seurat Object"),
                        downloadButton(outputId = ns("download_sce_object"), label = "Download SCE Object"),

                        ),
        #Outputs
        mainPanel(style = "overflow-y:scroll; max-height: 500px; position:relative;",

                    #textouput to show summary of current seurat object in use
                    h4("Summary Seurat: "),
                    verbatimTextOutput(outputId = ns("seurat_summary_output")),
                    verbatimTextOutput(outputId = ns("seurat_metadata_output")),
                    #textoutput to show summary of current sce object in use
                    h4("Summary SCE: "),
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

    #seuratObect; a Seurat object as reactiveValue which can be reassigned/updated as workflow progresses
    seuratObject <- reactiveVal()

    #sceObject; a SingleCellExperiment object as reactiveValue which is used throughout the workflow 
    sceObject <- reactiveVal()

    #plotObject for plotOutput
    plotObject <- reactiveVal()

    #rowNames/geneNames; store rownames to avoid issues during sce to seurat conversion and vice versa
    geneNamesSCE <- reactiveVal()
    geneNamesSeurat <- reactiveVal()


    observeEvent(input$sce_rds_file, {
    #create seurat object where needed #remove below line and integrate with UpdateSeurat function
        seuratObject(SCEtoSeurat(input$sce_rds_file$datapath))
        sceObject(RDStoSCE(input$sce_rds_file$datapath))
        #
        geneNamesSCE(RowNamesSCE(sceObject()))
        geneNamesSeurat(RowNamesSeurat(seuratObject()))
        #
        showNotification("Upload Complete")
    })

    observe({
        if (!is.null(input$sce_rds_file)) {
            updateSelectInput(session = session, inputId = "violin_plot_select_feature", choices = getFeatureNamesSeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat())))
            updateSelectInput(session = session, inputId = "feature_scatter_plot_select_feature", choices = getFeatureNamesSeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat())))
        }
    })

    observeEvent(input$plot_button, {
        if (!is.null(input$sce_rds_file) && input$plot_type == "Violin Plot") {
            plotObject(ViolinPlot(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), input$violin_plot_select_feature))
        }
        else if (!is.null(input$sce_rds_file) && input$plot_type == "Feature Scatter Plot") {
            plotObject(FeatureScatterPlot(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), input$feature_scatter_plot_select_feature))
        }
    })

    observeEvent(input$plot_hvg_button, {
        plotObject(PlotHVG(seuratObject()))
    })

    observeEvent(input$plot_pca_button, {
        
        plotObject(PlotReductionSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), "pca"))
    })

    observeEvent(input$plot_tsne_button, {
        plotObject(PlotReductionSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), "tsne"))
    })

    observeEvent(input$normalize_button, {
    if (!is.null(input$sce_rds_file)) {
        withProgress(message = "Normalizing", max = 1, value = 1, {
            seuratObject(Normalization(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), input$normalization_method, as.numeric(input$scale_factor)))
            })
        sceObject(UpdateSCE(sceObject(), geneNamesSeurat(), seuratObject(), "seuratNormalizedData", "data"))
        sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
        showNotification("Normalization Complete")
        }
    })

    observeEvent(input$scale_button, {
    if (!is.null(input$sce_rds_file)) {
        withProgress(message = "Scaling", max = 1, value = 1, {
            seuratObject(ScaleDataSeurat(UpdateSeurat(sceObject(), geneNamesSeurat())))
            })
        sceObject(UpdateSCE(sceObject(), geneNamesSeurat(), seuratObject(), "seuratScaledData", "scale.data"))
        sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
        showNotification("Scale Complete")
        }
    })

    observeEvent(input$run_pca_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Running PCA", max = 1, value = 1, {
                seuratObject(PCASeurat(UpdateSeurat(sceObject(), geneNamesSeurat())))
            })
        sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
        showNotification("PCA Complete")
        }
    })

    observeEvent(input$find_hvg_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Finding highly variable genes", max = 1, value = 1, {       
                seuratObject(FindHVG(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), input$hvg_method, as.numeric(input$hvg_no_features)))
            })
        sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
        showNotification("Find HVG Complete")
        }
    })

    observeEvent(input$find_clusters_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Finding clusters", max = 1, value = 1, {
                seuratObject(FindClustersSeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat())))
            })
        sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
        showNotification("Find Clusters Complete")
        }
    })

    observeEvent(input$run_tsne_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Running tSNE", max = 1, value = 1, {
                seuratObject(RunTSNESeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat())))
            })
        sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
        showNotification("tSNE Complete")
        }
    })

 


    output$seurat_metadata_output <- renderText({
        if (!is.null(input$sce_rds_file)) {
            #capture.output(seuratObject())
            capture.output(names(seuratObject()@meta.data))
        }
    })


    output$seurat_summary_output <- renderText({
    if (!is.null(input$sce_rds_file)) {
        capture.output(seuratObject())
        #capture.output(names(seuratObject()@meta.data))
    }
 })

    output$sce_summary_output <- renderText({
    if (!is.null(input$sce_rds_file)) {
        capture.output(sceObject())
    }
})

    output$download_seurat_object <- downloadHandler(
        filename = function() {
            "seuratObject.rds"
        },
        content = function(file) {
            saveRDS(UpdateSeurat(sceObject(),geneNamesSeurat()), file)
        }
        )

    output$download_sce_object <- downloadHandler(
        filename = function() {
            "sceObject.rds"
        },
        content = function(file) {
            saveRDS(sceObject(), file)
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

ScaleDataSeurat <- function(seuratObject) {
    return(ScaleData(seuratObject))
}

PCASeurat <- function(seuratObject) {
    return(RunPCA(seuratObject))
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

PlotReductionSeurat <- function(seuratObject, reduction) {

    plot <- DimPlot(seuratObject, reduction = reduction)
    if ("ident" %in% names(plot$data) && "seurat_clusters" %in% names(seuratObject@meta.data)) {
        plot$data$ident <- seuratObject@meta.data$seurat_clusters
    }
    return(plot)
}


UpdateSCE <- function(sce, geneNames, seuratObject, assaySlotSCE, assaySlotSeurat) {
    assay(sce, assaySlotSCE) <- NULL
    assay(sce, assaySlotSCE) <- slot(seuratObject@assays$RNA, assaySlotSeurat)
    rownames(sce) <- geneNames
    return(sce)
}

UpdateSeurat <- function(sce, geneNames) {
    seuratObject <- CreateSeuratObject(counts = counts(sce))
    if ("seuratNormalizedData" %in% names(assays(sce))) {
        seuratObject@assays$RNA@data <- assay(sce, "seuratNormalizedData")
        rownames(seuratObject@assays$RNA@data) <- geneNames
    }
    if ("seuratScaledData" %in% names(assays(sce))) {
        seuratObject@assays$RNA@scale.data <- assay(sce, "seuratScaledData")
        rownames(seuratObject@assays$RNA@data) <- geneNames
    }
    if (!is.null(sce@metadata[["seurat"]]) && length(sce@metadata[["seurat"]]@assays$RNA@var.features) > 0) {
        seuratObject@assays$RNA@var.features <- sce@metadata[["seurat"]]@assays$RNA@var.features
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$pca)) {
        seuratObject@reductions$pca <- sce@metadata[["seurat"]]@reductions$pca
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$tsne)) {
        seuratObject@reductions$tsne <- sce@metadata[["seurat"]]@reductions$tsne
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@meta.data)) {
        seuratObject@meta.data <- sce@metadata[["seurat"]]@meta.data
    }

    if (!is.null(seuratObject@meta.data$seurat_clusters)) {
      #print(head(seuratObject@meta.data$seurat_clusters)) #finalize this
    }

    return(seuratObject)
}

AddSeuratToMetaDataSCE <- function(sce, seuratObject) {
    seuratObject@assays$RNA@counts <- new("dgCMatrix")
    seuratObject@assays$RNA@data <- new("dgCMatrix")
    seuratObject@assays$RNA@scale.data <- matrix()
    #update one metadata or all seurat objects everytime?
    sce@metadata[["seurat"]] <- seuratObject
  
    return(sce)
}

RowNamesSeurat <- function(seuratObject) {
    return(rownames(seuratObject))
}

RowNamesSCE <- function(sce) {
    return(rownames(sce))
}


FindClustersSeurat <- function(seuratObject) {
    seuratObject <- FindNeighbors(seuratObject)
    return(FindClusters(seuratObject))
}

RunTSNESeurat <- function(seuratObject) {
    seurat <- RunTSNE(seuratObject)
    return(seurat)
}



# ----