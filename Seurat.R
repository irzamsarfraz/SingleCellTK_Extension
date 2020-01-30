library(shiny);
library(shinyWidgets);
library(shinyBS);
library(SingleCellExperiment);
library(Seurat);

# User Interface for Seurat Workflow ---
Seurat_UI <- function(id) {
    ns <- NS(id);
    fluidPage(
        bsCollapse(id = "SeuratUI", open = "Data Input",
            bsCollapsePanel("Data Input",
                fluidRow(
                    column(4,
                        panel(heading = "Input",
                            fileInput(inputId = ns("sce_rds_file"), label = "Input SingleCellExperiment RDS file", multiple = FALSE, accept = c(".rds"), buttonLabel = "Browse")
                            )
                        ),
                    column(8,
                        panel(heading = "Summary",
                            verbatimTextOutput(outputId = ns("seurat_summary_output"), placeholder = TRUE),
                            verbatimTextOutput(outputId = ns("seurat_metadata_output"), placeholder = TRUE)
                            )
                          )
                        ),
                            style = "primary"
                            ),

            bsCollapsePanel("Normalize Data",
                fluidRow(
                    column(4,
                        panel(
                            selectInput(inputId = ns("normalization_method"), label = "Select normalization method: ", choices = c("LogNormalize", "CLR", "RC")),
                            textInput(inputId = ns("scale_factor"), label = "Set scaling factor: ", value = "10000"),
                            actionButton(inputId = ns("normalize_button"), "Normalize")
                             )
                          )
                        ),
                            style = "primary"
                            ),

            bsCollapsePanel("Scale Data",
                fluidRow(
                    column(4,
                        panel(
                            selectInput(inputId = ns("model.use"), label = "Select model for scaling: ", choices = c("linear", "poisson", "negbinom")),
                            materialSwitch(inputId = ns("do.scale"), label = "Scale data?", value = TRUE),
                            materialSwitch(inputId = ns("do.center"), label = "Center data?", value = TRUE),
                            textInput(inputId = ns("scale.max"), label = "Max value for scaled data: ", value = "10"),
                            actionButton(inputId = ns("scale_button"), "Scale")
                            )
                          )
                        ),
                            style = "primary"
                            ),

            bsCollapsePanel("Highly Variable Genes",
                fluidRow(
                    column(4,
                        fluidRow(
                            column(12,
                                panel(heading = "Compute HVG",
                                    selectInput(inputId = ns("hvg_method"), label = "Select HVG method: ", choices = c("vst", "mean.var.plot", "dispersion")),
                                    textInput(inputId = ns("hvg_no_features"), label = "Select number of features to find: ", value = "2000"),
                                    actionButton(inputId = ns("find_hvg_button"), "Find HVG")
                                     )
                                  )
                                ),
                        br(),
                        fluidRow(
                            column(12,
                                panel(heading = "Display HVG",
                                    textInput(inputId = ns("hvg_no_features_view"), label = "Select number of features to display: ", value = "100"),
                                    verbatimTextOutput(outputId = ns("hvg_output"), placeholder = TRUE)
                                     )
                                  )
                                )
                          ),
                     column(8,
                        fluidRow(
                            column(12,
                                panel(heading = "Plot",
                                    plotOutput(outputId = ns("plot_hvg"))
                                     )
                                  )
                                )
                           )      
                    ),
                    style = "primary"),

            bsCollapsePanel("Dimensionality Reduction",
                tabsetPanel(type = "tabs",
                    tabPanel("PCA",
                        br(),
                        fluidRow(
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "PCA",
                                            textInput(inputId = ns("pca_no_components"), label = "Select number of components to compute: ", value = "50"),
                                            actionButton(inputId = ns("run_pca_button"), "Run PCA")
                                             )
                                          )
                                        )
                                  ),
                            column(8,
                                fluidRow(
                                    column(12,
                                        tabsetPanel(type = "tabs",
                                            tabPanel(title = "PCA Plot",
                                                panel(heading = "PCA Plot",
                                                    plotOutput(outputId = ns("plot_pca"))
                                                     )
                                                    ),
                                            tabPanel(title = "Elbow Plot",
                                                panel(heading = "Elbow Plot",
                                                    plotOutput(outputId = ns("plot_elbow"))
                                                     )
                                                    )
                                                   )
                                          )

                                        )
                                  )
                               )

                            ),
                    tabPanel("ICA",
                        br(),
                        fluidRow(
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "ICA",
                                            textInput(inputId = ns("ica_no_components"), label = "Select number of components to compute: ", value = "50"),
                                            actionButton(inputId = ns("run_ica_button"), "Run ICA")
                                             )
                                          )
                                        ),
                                  ),
                            column(8,
                                fluidRow(
                                    column(12,
                                        panel(heading = "Plot",
                                            plotOutput(outputId = ns("plot_ica"))))

                                        )
                                  )
                                )
                            )
                    ),
                    style = "primary"),



            bsCollapsePanel("tSNE/UMAP",
                tabsetPanel(type = "tabs",
                    tabPanel("tSNE",
                        br(),
                        fluidRow(
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "tSNE",
                                            selectInput(inputId = ns("reduction_tsne_method"), label = "Select reduction method: ", choices = c("pca", "ica")),
                                            textInput(inputId = ns("reduction_tsne_count"), label = "Select number of reduction components: ", value = "20"),
                                            actionButton(inputId = ns("run_tsne_button"), "Run tSNE")
                                             )
                                          )
                                        )
                                  ),
                            column(8,
                                fluidRow(
                                    panel(heading = "Plot",
                                        column(12,
                                            plotOutput(outputId = ns("plot_tsne"))
                                              )
                                         )
                                        )
                                  )
                                )
                            ),
                    tabPanel("UMAP",
                        br(),
                        fluidRow(
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "UMAP",
                                            selectInput(inputId = ns("reduction_umap_method"), label = "Select reduction method: ", choices = c("pca", "ica")),
                                            textInput(inputId = ns("reduction_umap_count"), label = "Select number of reduction components: ", value = "20"),
                                            actionButton(inputId = ns("run_umap_button"), "Run UMAP")
                                            )
                                          )
                                        )
                                 ),
                            column(8,
                                fluidRow(
                                    panel(heading = "Plot",
                                        column(12,
                                            plotOutput(outputId = ns("plot_umap"))
                                              )
                                         )
                                        )
                                  )
                                )
                            )
                    ),
                    style = "primary"),

            bsCollapsePanel("Clustering",
                fluidRow(
                    column(4,
                        fluidRow(
                            column(12,
                                panel(
                                    selectInput(inputId = ns("reduction_clustering_method"), label = "Select reduction method: ", choices = c("pca", "ica")),
                                    textInput(inputId = ns("reduction_clustering_count"), label = "Select number of reduction components: ", value = "20"),
                                    selectInput(inputId = ns("algorithm.use"), label = "Select model for scaling: ", choices = c("original Louvain algorithm", "Louvain algorithm with multilevel refinement", "SLM algorithm")),
                                    materialSwitch(inputId = ns("group.singletons"), label = "Group singletons?", value = TRUE),
                                    actionButton(inputId = ns("find_clusters_button"), "Find Clusters")
                                    )
                                   )
                                )
                          )
                    ),
                    style = "primary"),

            bsCollapsePanel("Download Data",
                fluidRow(
                    column(4,
                        fluidRow(
                            column(12,
                                panel(heading = "Seurat",
                                    verbatimTextOutput(outputId = ns("download_seurat_summary")),
                                    downloadButton(outputId = ns("download_seurat_object"), label = "Download Seurat Object")
                                    )
                                   )
                                )
                          ),
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "SCE",
                                            verbatimTextOutput(outputId = ns("download_sce_summary")),
                                            downloadButton(outputId = ns("download_sce_object"), label = "Download SCE Object")
                                            )
                                          )
                                        )
                                  )
                    ),
                    style = "primary")
       ),
    );
}
# ----

# Server for Seurat Workflow ---
Seurat_Server <- function(input, output, session, x) {
    ns <- session$ns

    #seuratObect; single Seurat object which is updated throughout the workflow 
    seuratObject <- reactiveVal()

    #sceObject; single SingleCellExperiment object which is updated throughout the workflow
    sceObject <- reactiveVal()

    #plotObject; contains objects for each of hvg plot, pca plot, ica plot, tsne plot and umap plot (each object accessed through $hvg, $pca, $ica, $tsne and $umap)
    plotObject <- reactiveValues()

    #numberOfReductionComponents; a reactive object to hold number of computed pca & ica components which are used to populate reduction components selectinput (accessed through $pca and $ica)
    numberOfReductionComponents <- reactiveValues()

    #geneNames object; store rownames/featurenames/genenames to avoid issues during sce to seurat conversion and vice versa
    geneNamesSCE <- reactiveVal()
    geneNamesSeurat <- reactiveVal()

    #server events
    observeEvent(input$sce_rds_file, {
        withProgress(message = "Uploading", max = 1, value = 1, {
            seuratObject(SCEtoSeurat(input$sce_rds_file$datapath))
            sceObject(RDStoSCE(input$sce_rds_file$datapath))
            geneNamesSCE(RowNamesSCE(sceObject()))
            geneNamesSeurat(RowNamesSeurat(seuratObject()))
        })
        showNotification("Upload complete")
    })

    output$hvg_output <- renderText({
        GetVariableFeaturesSeurat(seuratObject(), input$hvg_no_features_view)
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
                seuratObject(ScaleDataSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), input$model.use, input$do.scale, input$do.center, input$scale.max))
            })
            sceObject(UpdateSCE(sceObject(), geneNamesSeurat(), seuratObject(), "seuratScaledData", "scale.data"))
            sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
            showNotification("Scale Complete")
        }
    })

    observeEvent(input$run_pca_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Running PCA", max = 1, value = 1, {
                seuratObject(PCASeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), input$pca_no_components))
                numberOfReductionComponents$pca <- dim(seuratObject()[["pca"]])[2]
            })
            sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
            withProgress(message = "Plotting PCA", max = 1, value = 1, {
                plotObject$PCA <- PlotReductionSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), "pca")
            })
            withProgress(message = "Generating Elbow Plot", max = 1, value = 1, {
                plotObject$Elbow <- PlotElbowSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()))
            })
            showNotification("PCA Complete")
        }
    })

    observeEvent(input$run_ica_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Running ICA", max = 1, value = 1, {
                seuratObject(ICASeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), input$ica_no_components))
                numberOfReductionComponents$ica <- dim(seuratObject()[["ica"]])[2]
            })
            sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
            withProgress(message = "Plotting ICA", max = 1, value = 1, {
                plotObject$ICA <- PlotReductionSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), "ica")
            })
            showNotification("ICA Complete")
        }
    })

    observeEvent(input$find_hvg_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
                seuratObject(FindHVGSeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), input$hvg_method, as.numeric(input$hvg_no_features)))
            })
            sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
            withProgress(message = "Plotting HVG", max = 1, value = 1, {
                plotObject$HVG <- PlotHVGSeurat(seuratObject())
            })
            showNotification("Find HVG Complete")
        }
    })

    observeEvent(input$find_clusters_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Finding clusters", max = 1, value = 1, {
                seuratObject(FindClustersSeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), reduction = input$reduction_clustering_method, dims = input$reduction_clustering_count, algorithm = input$algorithm.use, group.singletons = input$group.singletons))
            })
            sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
            showNotification("Find Clusters Complete")
        }
    })

    observeEvent(input$run_tsne_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Running tSNE", max = 1, value = 1, {
                seuratObject(RunTSNESeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), input$reduction_tsne_method, input$reduction_tsne_count))
            })
            sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
            withProgress(message = "Plotting tSNE", max = 1, value = 1, {
                plotObject$TSNE <- PlotReductionSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), "tsne")
            })
            showNotification("tSNE Complete")
        }
    })

    observeEvent(input$run_umap_button, {
        if (!is.null(input$sce_rds_file)) {
            withProgress(message = "Running UMAP", max = 1, value = 1, {
                seuratObject(RunUMAPSeurat(seuratObject = UpdateSeurat(sceObject(), geneNamesSeurat()), input$reduction_umap_method, input$reduction_umap_count))
            })
            sceObject(AddSeuratToMetaDataSCE(sceObject(), seuratObject()))
            withProgress(message = "Plotting UMAP", max = 1, value = 1, {
                plotObject$UMAP <- PlotReductionSeurat(UpdateSeurat(sceObject(), geneNamesSeurat()), "umap")
            })
            showNotification("UMAP Complete")
        }
    })

    output$seurat_metadata_output <- renderText({
        if (!is.null(input$sce_rds_file)) {
            capture.output(names(seuratObject()@meta.data))
        }
    })

    output$seurat_summary_output <- renderText({
        if (!is.null(input$sce_rds_file)) {
            capture.output(seuratObject())
        }
    })

    output$download_seurat_summary <- renderPrint({
        if (!is.null(input$sce_rds_file)) {
            print("Seurat Object Size:")
            paste(as.double(object.size(seuratObject())) / 1000000, "Megabytes")
        }
    })

    output$download_sce_summary <- renderPrint({
        if (!is.null(input$sce_rds_file)) {
            print("SCE Object Size:")
            paste(as.double(object.size(sceObject())) / 1000000, "Megabytes")
        }
    })

    output$download_seurat_object <- downloadHandler(
        filename = function() {
            "seuratObject.rds"
        },
        content = function(file) {
            saveRDS(UpdateSeurat(sceObject(), geneNamesSeurat()), file)
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

    output$plot_hvg <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$HVG
        }
    })

    output$plot_pca <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$PCA
        }
    })

    output$plot_elbow <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$Elbow
        }
    })

    output$plot_ica <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$ICA
        }
    })

    output$plot_tsne <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$TSNE
        }
    })

    output$plot_umap <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$UMAP
        }
    })

    observe({
        if (input$reduction_umap_method == "pca") {
            updateTextInput(session = session, inputId = "reduction_umap_count", value = numberOfReductionComponents$pca)
            }
        else if (input$reduction_umap_method == "ica") {
            updateTextInput(session = session, inputId = "reduction_umap_count", value = numberOfReductionComponents$ica)
            }
        if (input$reduction_clustering_method == "pca") {
            updateTextInput(session = session, inputId = "reduction_clustering_count", value = numberOfReductionComponents$pca)
            }
        else if (input$reduction_clustering_method == "ica") {
            updateTextInput(session = session, inputId = "reduction_clustering_count", value = numberOfReductionComponents$ica)
            }
        if (input$reduction_tsne_method == "pca") {
            updateTextInput(session = session, inputId = "reduction_tsne_count", value = numberOfReductionComponents$pca)
            }
        else if (input$reduction_tsne_method == "ica") {
            updateTextInput(session = session, inputId = "reduction_tsne_count", value = numberOfReductionComponents$ica)
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
#' Normalize the seurat object (data) using the provided parameters and returns an updated seurat object
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


#' ScaleDataSeurat
#' Scales the seurat object (data) using the provided parameters and returns an updated seurat object
#' @param seuratObject
#' @param model.use
#' @param do.scale
#' @param do.center
#' @param scale.max
#'
#' @return seurat object
#' @export
#'
#' @examples
ScaleDataSeurat <- function(seuratObject, model.use, do.scale, do.center, scale.max) {
    return(ScaleData(seuratObject, model.use = model.use, do.scale = do.scale, do.center = do.center, scale.max = as.double(scale.max)))
}

#' PCASeurat
#' Computes the principal component analysis on the provided seurat object (data) and returns the updated seurat object
#' @param seuratObject
#' @param npcs
#'
#' @return seurat object
#' @export
#'
#' @examples
PCASeurat <- function(seuratObject, npcs) {
    return(RunPCA(seuratObject, npcs = as.double(npcs)))
}

#' ICASeurat
#' Computes the independent component analysis on the provided seurat object (data) and returns the updated seurat object
#' @param seuratObject
#' @param nics
#'
#' @return seurat object
#' @export
#'
#' @examples
ICASeurat <- function(seuratObject, nics) {
    return(RunICA(seuratObject, nics = as.double(nics)))
}

#' FindHVGSeurat
#' Find highly variable genes using of the selected methods
#' @param seuratObject
#' @param hvgMethod
#' @param hvgNumber
#'
#' @return seurat object
#' @export
#'
#' @examples
FindHVGSeurat <- function(seuratObject, hvgMethod, hvgNumber) {
    return(FindVariableFeatures(seuratObject, selection.method = hvgMethod, nfeatures = hvgNumber))
}


#' PlotHVGSeurat
#' Plot the highly variable genes
#' @param seuratObject
#'
#' @return plot object
#' @export
#'
#' @examples
PlotHVGSeurat <- function(seuratObject) {
    return(VariableFeaturePlot(seuratObject))
}

#' PlotReductionSeurat
#' Plots the dimensionality reduction algorithms i.e. pca, ica, tsne and umap from the input seurat object
#' @param seuratObject
#' @param reduction
#'
#' @return plot object
#' @export
#'
#' @examples
PlotReductionSeurat <- function(seuratObject, reduction) {

    plot <- DimPlot(seuratObject, reduction = reduction)
    if ("ident" %in% names(plot$data) && "seurat_clusters" %in% names(seuratObject@meta.data)) {
        plot$data$ident <- seuratObject@meta.data$seurat_clusters
    }
    return(plot)
}


#' UpdateSCE
#' Modifies the input sce object to include the updated assays from seurat object
#' @param sce
#' @param geneNames
#' @param seuratObject
#' @param assaySlotSCE
#' @param assaySlotSeurat
#'
#' @return singlecellexperiment object
#' @export
#'
#' @examples
UpdateSCE <- function(sce, geneNames, seuratObject, assaySlotSCE, assaySlotSeurat) {
    assay(sce, assaySlotSCE) <- NULL
    assay(sce, assaySlotSCE) <- slot(seuratObject@assays$RNA, assaySlotSeurat)
    rownames(sce) <- geneNames
    return(sce)
}

#' UpdateSeurat
#' Modifies the seurat object to include updated assays from input sce object
#' @param sce
#' @param geneNames
#'
#' @return seurat object
#' @export
#'
#' @examples
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
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$ica)) {
        seuratObject@reductions$ica <- sce@metadata[["seurat"]]@reductions$ica
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$tsne)) {
        seuratObject@reductions$tsne <- sce@metadata[["seurat"]]@reductions$tsne
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@reductions$umap)) {
        seuratObject@reductions$umap <- sce@metadata[["seurat"]]@reductions$umap
    }
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@meta.data)) {
        seuratObject@meta.data <- sce@metadata[["seurat"]]@meta.data
    }
    return(seuratObject)
}

#' AddSeuratToMetaDataSCE
#' Adds seurat object to the metadata slot of singlecellexperiment object after removing the data matrices so information in the seurat object is not lost
#' @param sce
#' @param seuratObject
#'
#' @return singlecellexperiment object
#' @export
#'
#' @examples
AddSeuratToMetaDataSCE <- function(sce, seuratObject) {
    seuratObject@assays$RNA@counts <- new("dgCMatrix")
    seuratObject@assays$RNA@data <- new("dgCMatrix")
    seuratObject@assays$RNA@scale.data <- matrix()
    sce@metadata[["seurat"]] <- seuratObject
    return(sce)
}

#' RowNamesSeurat
#' Retrieves a list of genenames/rownames/featurenames from seurat object
#' @param seuratObject
#'
#' @return list of rownames
#' @export
#'
#' @examples
RowNamesSeurat <- function(seuratObject) {
    return(rownames(seuratObject))
}

#' RowNamesSCE
#' Retrieves a list of genenames/rownames/featurenames from singlecellexperiment object
#' @param sce
#'
#' @return
#' @export
#'
#' @examples
RowNamesSCE <- function(sce) {
    return(rownames(sce))
}


#' FindClustersSeurat
#' Computes the clusters from the given seurat object
#' @param seuratObject
#' @param reduction
#' @param dims
#' @param algorithm
#' @param group.singletons
#'
#' @return seurat object
#' @export
#'
#' @examples
FindClustersSeurat <- function(seuratObject, reduction, dims, algorithm, group.singletons) {
    seuratObject <- FindNeighbors(seuratObject, reduction = reduction, dims = 1:dims)
    no_algorithm <- 1
    print(algorithm)
    if (algorithm == "original Louvain algorithm") {
        no_algorithm = 1
    } else if (algorithm == "Louvain algorithm with multilevel refinement") {
        no_algorithm = 2
    } else if (algorithm == "SLM algorithm") {
        no_algorithm = 3
    }
    return(FindClusters(seuratObject, algorithm = no_algorithm, group.singletons = group.singletons))
}

#' RunTSNESeurat
#' Computes tSNE from the given seurat object
#' @param seuratObject
#' @param reduction
#' @param dims
#'
#' @return seurat object
#' @export
#'
#' @examples
RunTSNESeurat <- function(seuratObject, reduction, dims) {
    print(seuratObject@reductions)
    seurat <- RunTSNE(seuratObject, reduction = reduction, dims = 1:dims)
    return(seurat)
}

#' RunUMAPSeurat
#' Computes UMAP from the given seurat object
#' @param seuratObject
#' @param reduction
#' @param dims
#'
#' @return seurat object
#' @export
#'
#' @examples
RunUMAPSeurat <- function(seuratObject, reduction, dims) {
    print(seuratObject@reductions)
    seurat <- RunUMAP(seuratObject, reduction = reduction, dims = 1:dims)
    return(seurat)
}

#' GetVariableFeaturesSeurat
#' Retrieves the requested number of variable features (feature names)
#' @param seuratObject
#' @param numberOfFeatures
#'
#' @return
#' @export
#'
#' @examples
GetVariableFeaturesSeurat <- function(seuratObject, numberOfFeatures) {
    if (length(seuratObject@assays$RNA@var.features) > 0) {
        return(print(seuratObject@assays$RNA@var.features[1:numberOfFeatures]))
    }
}

#' PlotElbowSeurat
#' Computes the plot object for elbow plot from the pca slot in the provided seurat object
#' @param seuratObject
#'
#' @return plot object
#' @export
#'
#' @examples
PlotElbowSeurat <- function(seuratObject) {
    return(ElbowPlot(seuratObject))
}

# ----
