library(shiny);
library(shinyWidgets);
library(shinyBS);
library(SingleCellExperiment);
library(Seurat);
library(cowplot);
library(svgPanZoom);
library(shinyjqui);
library(ggplotify);
library(ggplot2);
library(shinyjs);
# User Interface for Seurat Workflow ---
Seurat_UI <- function(id) {
    ns <- NS(id);
    fluidPage(
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary"="border-color:#dddddd;",".panel-primary>.panel-heading+.panel-collapse>.panel-body"="border-color:#dddddd;")),
        bsCollapse(id = ns("SeuratUI"), open = "Data Input",
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
                                            textInput(inputId = ns("pca_no_components"), label = "Select number of components to compute: ", value = "20"),
                                            materialSwitch(inputId = ns("pca_compute_elbow"), label = "Compute ElbowPlot?", value = TRUE),
                                            materialSwitch(inputId = ns("pca_compute_jackstraw"), label = "Compute JackStrawPlot?", value = TRUE),
                                            materialSwitch(inputId = ns("pca_compute_heatmap"), label = "Compute Heatmap?", value = TRUE),
                                            actionButton(inputId = ns("run_pca_button"), "Run PCA")
                                             ),
                                        panel(heading = "Select No. of Components",
                                            h5("Number of components suggested by ElbowPlot: "),
                                            verbatimTextOutput(outputId = ns("pca_significant_pc_output"), placeholder = TRUE),
                                            sliderInput(inputId = ns("pca_significant_pc_slider"), label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10, round = TRUE)
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
                                                    ),
                                            tabPanel(title = "JackStraw Plot",
                                                panel(heading = "JackStraw Plot",
                                                    plotOutput(outputId = ns("plot_jackstraw"))
                                                     )
                                                    ),
                                            tabPanel(title = "Heatmap Plot",
                                                panel(heading = "Heatmap Plot",
                                                    panel(heading = "Plot Options",
                                                        fluidRow(
                                                            column(6,
                                                                pickerInput(inputId = ns("picker_dimheatmap_components_pca"), label = "Select principal components to plot:", choices = c(), options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), multiple = TRUE)
                                                            ),
                                                            column(6,
                                                                sliderInput(inputId = ns("slider_dimheatmap_pca"), label = "Number of columns for the plot: ", min = 1, max = 4, value = 2)
                                                            )
                                                        ),    
                                                        actionButton(inputId = ns("plot_heatmap_pca_button"), "Plot")
                                                        ),
                                                    panel(heading = "Plot",
                                                        jqui_resizable(plotOutput(outputId = ns("plot_heatmap")), options=list(maxWidth=700))
                                                        )
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
                                            textInput(inputId = ns("ica_no_components"), label = "Select number of components to compute: ", value = "20"),
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
            seuratObject(.sceToSeurat(input$sce_rds_file$datapath))
            sceObject(.rdsToSce(input$sce_rds_file$datapath))
            geneNamesSCE(.rowNamesSCE(sceObject()))
            geneNamesSeurat(.rowNamesSeurat(seuratObject()))
        })
        updateCollapse(session = session, "SeuratUI", style = list("Data Input" = "danger"))
        showNotification("Upload complete")
    })

    output$hvg_output <- renderText({
        if (!is.null(sceObject())) {
            if (!is.null(sceObject()@metadata[["seurat"]])) {
                if (length(slot(sceObject()@metadata[["seurat"]], "assays")[["RNA"]]@var.features) > 0) {
                    .seuratGetVariableFeatures(sceObject(), geneNamesSeurat(), input$hvg_no_features_view)
                }
            }
        }
    })

    observeEvent(input$normalize_button, {
        if (!is.null(sceObject())) {
            withProgress(message = "Normalizing", max = 1, value = 1, {
                sceObject(seuratNormalizeData(sce = sceObject(), geneNames = geneNamesSeurat(), input$normalization_method, as.numeric(input$scale_factor)))
            })
            updateCollapse(session = session, "SeuratUI", style = list("Normalize Data" = "danger"))
            showNotification("Normalization Complete")
        }
        else {
            showNotification("Please input dataset (rds file) before normalizing data!", type = "error")
        }
    })

    observeEvent(input$scale_button, {
        if (!is.null(sceObject())) {
            withProgress(message = "Scaling", max = 1, value = 1, {
                sceObject(seuratScaleData(sceObject(), geneNamesSeurat(), input$model.use, input$do.scale, input$do.center, input$scale.max))
            })
            updateCollapse(session = session, "SeuratUI", style = list("Scale Data" = "danger"))
            showNotification("Scale Complete")
        }
        else {
            showNotification("Please input dataset (rds file) before scaling data!", type = "error")
        }
    })

    observeEvent(input$run_pca_button, {
        if (!is.null(sceObject())) {
            if (!is.null(sceObject()@metadata[["seurat"]])) {
                if ((length(slot(sceObject()@metadata[["seurat"]], "assays")[["RNA"]]@scale.data) > 0) && (length(slot(sceObject()@metadata[["seurat"]], "assays")[["RNA"]]@var.features) > 0)) {
                    withProgress(message = "Running PCA", max = 1, value = 1, {
                        sceObject(seuratPCA(sceObject(), geneNamesSeurat(), input$pca_no_components))
                        numberOfReductionComponents$pca <- dim(convertSCEToSeurat(sceObject(), geneNamesSeurat())[["pca"]])[2]
                    })
                    withProgress(message = "Plotting PCA", max = 1, value = 1, {
                        plotObject$PCA <- seuratReductionPlot(sceObject(), geneNamesSeurat(), "pca")
                    })
                    if (input$pca_compute_elbow) {
                        withProgress(message = "Generating Elbow Plot", max = 1, value = 1, {
                            numberOfReductionComponents$significantPC <- .computeSignificantPC(sceObject(), geneNamesSeurat())
                            plotObject$Elbow <- seuratElbowPlot(sceObject(), geneNamesSeurat(), numberOfReductionComponents$significantPC)
                        })
                    }
                    if (input$pca_compute_jackstraw) {
                        withProgress(message = "Generating JackStraw Plot", max = 1, value = 1, {
                            sceObject(seuratComputeJackStraw(sceObject(), geneNamesSeurat(), input$pca_no_components))
                            plotObject$JackStraw <- seuratJackStrawPlot(sceObject(), geneNamesSeurat(), input$pca_no_components)
                        })
                    }
                    if (input$pca_compute_heatmap) {
                        withProgress(message = "Generating Heatmap Plot", max = 1, value = 1, {
                            plotObject$HeatmapCompute <- seuratComputeHeatmap(sceObject(), geneNamesSeurat(), input$pca_no_components)
                            updatePickerInput(session = session, inputId = "picker_dimheatmap_components_pca", choices = .getPCAComponentNames(numberOfReductionComponents$pca))
                        })
                    }
                    addTooltip(session = session, id = ns("reduction_tsne_count"), paste("Maximum components available:", numberOfReductionComponents$pca), placement = "bottom", trigger = "hover")
                    addTooltip(session = session, id = ns("reduction_umap_count"), paste("Maximum components available:", numberOfReductionComponents$pca), placement = "bottom", trigger = "hover")
                    addTooltip(session = session, id = ns("reduction_clustering_count"), paste("Maximum components available:", numberOfReductionComponents$pca), placement = "bottom", trigger = "hover")
                    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "danger"))
                    showNotification("PCA Complete")
                }
                else {
                    showNotification("Please scale data and find highly variable genes before computing pca!", type = "error")
                }
            }
        }
    })

    observeEvent(input$run_ica_button, {
         if (!is.null(sceObject())) {
            if (!is.null(sceObject()@metadata[["seurat"]])) {
                if ((length(slot(sceObject()@metadata[["seurat"]], "assays")[["RNA"]]@scale.data) > 0) && (length(slot(sceObject()@metadata[["seurat"]], "assays")[["RNA"]]@var.features) > 0)) {
                    withProgress(message = "Running ICA", max = 1, value = 1, {
                        sceObject(seuratICA(sceObject(), geneNamesSeurat(), input$ica_no_components))
                        numberOfReductionComponents$ica <- dim(convertSCEToSeurat(sceObject(), geneNamesSeurat())[["ica"]])[2]
                    })
                    withProgress(message = "Plotting ICA", max = 1, value = 1, {
                        plotObject$ICA <- seuratReductionPlot(sceObject(), geneNamesSeurat(), "ica")
                    })
                    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "danger"))
                    showNotification("ICA Complete")
                }
                else {
                    showNotification("Please scale data and find highly variable genes before computing ica!", type = "error")
                }
            }
        }
    })

    observeEvent(input$find_hvg_button, {
        if (!is.null(sceObject())) {
            withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
                sceObject(seuratFindHVG(sceObject(), geneNamesSeurat(), input$hvg_method, as.numeric(input$hvg_no_features)))
            })
            withProgress(message = "Plotting HVG", max = 1, value = 1, {
                plotObject$HVG <- seuratPlotHVG(sceObject(), geneNamesSeurat())
            })
            updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "danger"))
            showNotification("Find HVG Complete")
        }
        else {
            showNotification("Please input dataset (rds file) before computing highly variable genes!", type = "error")
        }
    })

    observeEvent(input$find_clusters_button, {
        if (!is.null(sceObject())) {
            if (!is.null(sceObject()@metadata[["seurat"]])) {
                if (!is.null(slot(sceObject()@metadata[["seurat"]], "reductions")[[input$reduction_clustering_method]])) {
                    withProgress(message = "Finding clusters", max = 1, value = 1, {
                        sceObject(seuratFindClusters(sceObject(), geneNamesSeurat(), reduction = input$reduction_clustering_method, dims = input$reduction_clustering_count, algorithm = input$algorithm.use, group.singletons = input$group.singletons))
                    })
                    updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "danger"))
                    showNotification("Find Clusters Complete")
                }
                else {
                    showNotification("Please compute pca/ica before processing clusters!", type = "error")
                }
            }
        }
    })

    observeEvent(input$run_tsne_button, {
        if (!is.null(sceObject())) {
            if (!is.null(sceObject()@metadata[["seurat"]])) {
                if (!is.null(slot(sceObject()@metadata[["seurat"]], "reductions")[[input$reduction_tsne_method]])) {
                    withProgress(message = "Running tSNE", max = 1, value = 1, {
                        sceObject(seuratRunTSNE(sceObject(), geneNamesSeurat(), input$reduction_tsne_method, input$reduction_tsne_count))
                    })
                    withProgress(message = "Plotting tSNE", max = 1, value = 1, {
                        plotObject$TSNE <- seuratReductionPlot(sceObject(), geneNamesSeurat(), "tsne")
                    })
                    updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "danger"))
                    showNotification("tSNE Complete")
                }
                else {
                    showNotification("Please compute pca/ica before processing tsne!", type = "error")
                }
            }
        }
    })

    observeEvent(input$run_umap_button, {
        if (!is.null(sceObject())) {
            if (!is.null(sceObject()@metadata[["seurat"]])) {
                if (!is.null(slot(sceObject()@metadata[["seurat"]], "reductions")[[input$reduction_umap_method]])) {
                    withProgress(message = "Running UMAP", max = 1, value = 1, {
                        sceObject(seuratRunUMAP(sceObject(), geneNamesSeurat(), input$reduction_umap_method, input$reduction_umap_count))
                    })
                    withProgress(message = "Plotting UMAP", max = 1, value = 1, {
                        plotObject$UMAP <- seuratReductionPlot(sceObject(), geneNamesSeurat(), "umap")
                    })
                    updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "danger"))
                    showNotification("UMAP Complete")
                }
                else {
                    showNotification("Please compute pca/ica before processing umap!", type = "error")
                }
            }
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

    output$pca_significant_pc_output <- renderText({
        if (!is.null(numberOfReductionComponents$significantPC)) {
            numberOfReductionComponents$significantPC
        }
    })

    observe({
        if (!is.null(numberOfReductionComponents$pca)) {
            updateSliderInput(session = session, inputId = "pca_significant_pc_slider", max = numberOfReductionComponents$pca)
        }
    })

    observe({
        if (!is.null(numberOfReductionComponents$significantPC)) {
            updateSliderInput(session = session, inputId = "pca_significant_pc_slider", value = numberOfReductionComponents$significantPC)
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
            saveRDS(convertSCEToSeurat(sceObject(), geneNamesSeurat()), file)
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

    output$plot_jackstraw <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$JackStraw
        }
    })

    output$plot_heatmap <- renderPlot({
        if (!is.null(input$sce_rds_file)) {
            plotObject$Heatmap
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
            updateTextInput(session = session, inputId = "reduction_umap_count", value = input$pca_significant_pc_slider)
        }
        else if (input$reduction_umap_method == "ica") {
            updateTextInput(session = session, inputId = "reduction_umap_count", value = numberOfReductionComponents$ica)
        }
        if (input$reduction_clustering_method == "pca") {
            updateTextInput(session = session, inputId = "reduction_clustering_count", value = input$pca_significant_pc_slider)
        }
        else if (input$reduction_clustering_method == "ica") {
            updateTextInput(session = session, inputId = "reduction_clustering_count", value = numberOfReductionComponents$ica)
        }
        if (input$reduction_tsne_method == "pca") {
            updateTextInput(session = session, inputId = "reduction_tsne_count", value = input$pca_significant_pc_slider)
        }
        else if (input$reduction_tsne_method == "ica") {
            updateTextInput(session = session, inputId = "reduction_tsne_count", value = numberOfReductionComponents$ica)
        }
    })
 
   observeEvent(input$plot_heatmap_pca_button, {
        if (!is.null(input$picker_dimheatmap_components_pca)) {
            plotObject$Heatmap <- seuratHeatmapPlot(plotObject$HeatmapCompute, length(input$picker_dimheatmap_components_pca), input$slider_dimheatmap_pca, input$picker_dimheatmap_components_pca)
        }
   })
}
# ----

# Helper/Wrapper Functions ---

#' .getPCAComponentNames
#' Creates a list of PC components to populate the picker for PC heatmap generation
#' @param maxComponents; number of components to return for the picker
#' @return componentNames; list of component names (appended with "PC")
.getPCAComponentNames <- function(maxComponents) {
    componentNames <- list()
    for (i in 1:maxComponents) {
        componentNames[i] <- paste0("PC",i)
    }
    return(componentNames)
}

#' .rdsToSce
#' Reads rds file (from a local path) and loads into sce object 
#' *Only to be used for first time initialization of the rds file into R environment*
#' @param filePath; path of the rds file to load
#' @return sce object
.rdsToSce <- function(filePath) {
    sce <- readRDS(filePath)
    return(sce)
}

#' .sceToSeurat
#' Converts a sce object to seurat object (using rds filepath)
#' *Only to be used for first time initialization of seurat object*
#' @param filePath; path of the rds file to convert to seurat object
#' @return seurat object
.sceToSeurat <- function(filePath) {
    seuratObject <- CreateSeuratObject(counts = counts(.rdsToSce(filePath)))
    return(seuratObject)
}

#' .addSeuratToMetaDataSCE
#' Adds the input seurat object to the metadata slot of the input sce object (after removing the data matrices)
#' @param sce; sce object to which seurat object should be added in the metadata slot (copy to)
#' @param seuratObject; seurat object which should be added to the metadata slot of sce object (copy from)
#' @return sce; updated sce object which now contains the seurat object in its metadata slot (excluding data matrices)
.addSeuratToMetaDataSCE <- function(sce, seuratObject) {
    seuratObject@assays$RNA@counts <- new("dgCMatrix")
    seuratObject@assays$RNA@data <- new("dgCMatrix")
    seuratObject@assays$RNA@scale.data <- matrix()
    sce@metadata[["seurat"]] <- seuratObject
    return(sce)
}

#' .rowNamesSeurat
#' Retrieves a list of genenames/rownames/featurenames from seurat object
#' @param seuratObject; seurat object from which the genenames/rownames/featurenames should be extracted
#' @return list() of genenames/rownames/featurenames
.rowNamesSeurat <- function(seuratObject) {
    return(rownames(seuratObject))
}

#' .rowNamesSCE
#' Retrieves a list of genenames/rownames/featurenames from sce object
#' @param sce; sce object from which the genenames/rownames/featurenames should be extracted
#' @return list() of genenames/rownames/featurenames
.rowNamesSCE <- function(sce) {
    return(rownames(sce))
}

#' .computeSignificantPC
#' Computes the significant principal components from an input sce object (must containt pca slot) using stdev
#' @param sceObject; sce object with pca computed
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return max_components; a numerical value indicating how many number of components are considered significant
.computeSignificantPC <- function(sceObject, geneNamesSeurat) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    max_components <- 0
    for (i in 1:(length(seuratObject[["pca"]]@stdev) - 1)) {
        if (abs(seuratObject[["pca"]]@stdev[i + 1] - seuratObject[["pca"]]@stdev[i]) > 0.1) {
            max_components <- i
        }
    }
    return(max_components)
}

#' seuratNormalizeData
#' Wrapper for NormalizeData() function from seurat library
#' Normalizes the sce object according to the input parameters 
#' @param sceObject; sce object to normalize
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param normalizationMethod; selected normalization method (default is "LogNormalize")
#' @param scaleFactor; numeric value that represents the scaling factor (default is 10000)
#' @return sceObject; normalized sce object
#' @export
seuratNormalizeData <- function(sceObject, geneNamesSeurat, normalizationMethod, scaleFactor) {
    seuratObject <- NormalizeData(convertSCEToSeurat(sceObject, geneNamesSeurat), normalization.method = normalizationMethod, scale.factor = scaleFactor)
    sceObject <- convertSeuratToSCE(sceObject, geneNamesSeurat, seuratObject, "seuratNormalizedData", "data")
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratScaleData
#' Scales the input sce object according to the input parameters
#' @param sceObject; sce object to scale
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param model.use; selected model to use for scaling data (default is "linear")
#' @param do.scale; boolean if data should be scaled or not (TRUE or FALSE, default is TRUE)
#' @param do.center; boolean if data should be centered or not (TRUE or FALSE, default is TRUE)
#' @param scale.max; maximum numeric value to return for scaled data (default is 10)
#' @return sceObject; scaled sce object
#' @export
seuratScaleData <- function(sceObject, geneNamesSeurat, model.use, do.scale, do.center, scale.max) {
    seuratObject <- ScaleData(convertSCEToSeurat(sceObject, geneNamesSeurat), model.use = model.use, do.scale = do.scale, do.center = do.center, scale.max = as.double(scale.max))
    sceObject <- convertSeuratToSCE(sceObject, geneNamesSeurat, seuratObject, "seuratScaledData", "scale.data")
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratPCA
#' Computes PCA on the input sce object and stores the calculated principal components within the sce object
#' @param sceObject; sce object on which to compute PCA
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param npcs; numeric value of how many components to compute (default is 20)
#' @return sceObject; updated sce object which now contains the computed principal components
#' @export
seuratPCA <- function(sceObject, geneNamesSeurat, npcs) {
    seuratObject <- RunPCA(convertSCEToSeurat(sceObject, geneNamesSeurat), npcs = as.double(npcs))
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratICA
#' Computes ICA on the input sce object and stores the calculated independent components within the sce object
#' @param sceObject; sce object on which to compute ICA
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param nics; numeric value of how many components to compute (default is 20)
#' @return sceObject; updated sce object which now contains the computed independent components
#' @export
seuratICA <- function(sceObject, geneNamesSeurat, nics) {
    seuratObject <- RunICA(convertSCEToSeurat(sceObject, geneNamesSeurat), nics = as.double(nics))
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratComputeJackStraw
#' Compute jackstraw plot and store the computations in the input sce object
#' @param sceObject; sce object on which to compute and store jackstraw plot
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims; numeric value of how many components to use for jackstraw plot (default = number of computed principal components)
#' @return sceObject; updated sce object with jackstraw computations stored in it
#' @export
seuratComputeJackStraw <- function(sceObject, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- JackStraw(seuratObject, dims = as.double(dims))
    seuratObject <- ScoreJackStraw(seuratObject, dims = 1:dims)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratFindHVG
#' Find highly variable genes and store in the input sce object
#' @param sceObject; sce object to compute highly variable genes from and to store back to it
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param hvgMethod; selected method to use for computation of highly variable genes (default is "vst")
#' @param hvgNumber; numeric value of how many genes to select as highly variable (default is 2000)
#' @return sceObject; updated sce object with highly variable genes computation stored
#' @export
seuratFindHVG <- function(sceObject, geneNamesSeurat, hvgMethod, hvgNumber) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- FindVariableFeatures(seuratObject, selection.method = hvgMethod, nfeatures = hvgNumber)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratPlotHVG
#' Plot highly variable genes from input sce object (must have highly variable genes computations stored)
#' @param sceObject; sce object that contains the highly variable genes computations
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return plot object 
#' @export
seuratPlotHVG <- function(sceObject, geneNamesSeurat) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    return(VariableFeaturePlot(seuratObject))
}

#' seuratReductionPlot
#' Plots the selected dimensionality reduction method
#' @param sceObject; sce object which has the selected dimensionality reduction algorithm already computed and stored
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction; one of selected algorithm from pca, ica, tsne and umap
#' @return plot object
#' @export
seuratReductionPlot <- function(sceObject, geneNamesSeurat, reduction) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    plot <- DimPlot(seuratObject, reduction = reduction)
    if ("ident" %in% names(plot$data) && "seurat_clusters" %in% names(seuratObject@meta.data)) {
        plot$data$ident <- seuratObject@meta.data$seurat_clusters
    }
    return(plot)
}


#' convertSeuratToSCE
#' Modifies the input sce object to include the updated assays from seurat object
#' @param sce; outdated sce object in which we wish to include the newly calculated assays from seurat object
#' @param geneNames; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param seuratObject; updated seurat object that contains the newly calculated assays that we wish to store in the sce object
#' @param assaySlotSCE; the relevant assay slot in sce object (copy to)
#' @param assaySlotSeurat; the relevant assay slow in seurat object (copy from)
#' @return sce; sce object that contains the newly added/modified assays from the seurat object
#' @export
convertSeuratToSCE <- function(sce, geneNames, seuratObject, assaySlotSCE, assaySlotSeurat) {
    assay(sce, assaySlotSCE) <- NULL
    assay(sce, assaySlotSCE) <- slot(seuratObject@assays$RNA, assaySlotSeurat)
    rownames(sce) <- geneNames
    return(sce)
}

#' convertSCEToSeurat
#' Converts sce object to seurat while retaining all assays and metadata
#' @param sce; sce object to convert to seurat
#' @param geneNames; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @return seuratObject; updated seurat object that contains all data from the input sce object
#' @export
convertSCEToSeurat <- function(sce, geneNames) {
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
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@assays$RNA@meta.features)) {
        seuratObject@assays$RNA@meta.features <- sce@metadata[["seurat"]]@assays$RNA@meta.features
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
    if (!is.null(sce@metadata[["seurat"]]) && !is.null(sce@metadata[["seurat"]]@commands)) {
        seuratObject@commands <- sce@metadata[["seurat"]]@commands
    }
    return(seuratObject)
}

#' seuratFindClusters
#' Computes the clusters from the input sce object and stores them back in sce object
#' @param sceObject; sce object from which clusters should be computed and stored in
#' @param reduction; selected reduction method to use for computing clusters ("pca" or "ica", default is "pca")
#' @param dims; numeric value of how many components to use for computing clusters (default is 10)
#' @param algorithm; selected algorithm to compute clusters (default is "original Louvain algorithm")
#' @param group.singletons; boolean if singletons should be grouped together or not (TRUE or FALSE, default is TRUE)
#' @return sceObject; updated sce object which now contains the computed clusters
#' @export
seuratFindClusters <- function(sceObject, geneNamesSeurat, reduction, dims, algorithm, group.singletons) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- FindNeighbors(seuratObject, reduction = reduction, dims = 1:dims)
    no_algorithm <- 1
    if (algorithm == "original Louvain algorithm") {
        no_algorithm = 1
    } else if (algorithm == "Louvain algorithm with multilevel refinement") {
        no_algorithm = 2
    } else if (algorithm == "SLM algorithm") {
        no_algorithm = 3
    }
    seuratObject <- FindClusters(seuratObject, algorithm = no_algorithm, group.singletons = group.singletons)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratRunTSNE
#' Computes tSNE from the given sce object and stores the tSNE computations back into the sce object
#' @param sceObject; sce object on which to compute the tSNE
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction; selected reduction algorithm to use for computing tSNE ("pca" or "ica", default is "pca")
#' @param dims; numerical value of how many reduction components to use for tSNE computation (default is 10)
#' @return sceObject; updated sce object with tSNE computations stored
#' @export
seuratRunTSNE <- function(sceObject, geneNamesSeurat, reduction, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- RunTSNE(seuratObject, reduction = reduction, dims = 1:dims)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' seuratRunUMAP
#' Computes UMAP from the given sce object and stores the UMAP computations back into the sce object
#' @param sceObject; sce object on which to compute the UMAP
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param reduction; selected reduction algorithm to use for computing UMAP ("pca" or "ica", default is "pca")
#' @param dims; numerical value of how many reduction components to use for UMAP computation (default is 10)
#' @return sceObject; updated sce object with UMAP computations stored
#' @export
seuratRunUMAP <- function(sceObject, geneNamesSeurat, reduction, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    seuratObject <- RunUMAP(seuratObject, reduction = reduction, dims = 1:dims)
    sceObject <- .addSeuratToMetaDataSCE(sceObject, seuratObject)
    return(sceObject)
}

#' .seuratGetVariableFeatures
#' Retrieves the requested number of variable feature names
#' @param sceObject; sce object from which to extract the variable feature names
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param numberOfFeatures; numerical value indicating how many feature names should be retrieved (default is 100)
#' @return list() of variable feature names
.seuratGetVariableFeatures <- function(sceObject, geneNamesSeurat, numberOfFeatures) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    if (length(seuratObject@assays$RNA@var.features) > 0) {
        return(print(seuratObject@assays$RNA@var.features[1:numberOfFeatures]))
    }
}

#' seuratElbowPlot
#' Computes the plot object for elbow plot from the pca slot in the input sce object
#' @param sceObject; sce object from which to compute the elbow plot (pca should be computed)
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param significantPC; a numerical value indicating the number of significant principal components (used to alter the color of the significant components)
#' @return plot object
#' @export
seuratElbowPlot <- function(sceObject, geneNamesSeurat, significantPC) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    plot <- ElbowPlot(seuratObject)
    plot <- ggplot_build(plot)
    for (i in 1:significantPC) {
        plot$data[[1]]$shape[i] <- 16
        plot$data[[1]]$colour[i] <- "red"
        plot$data[[1]]$size[i] <- 3.5
    }
    plot <- ggplot_gtable(plot)
    plot <- as.ggplot(plot)
    return(plot)
}

#' seuratJackStrawPlot
#' Computes the plot object for jackstraw plot from the pca slot in the input sce object
#' @param sceObject; sce object from which to compute the jackstraw plot (pca should be computed)
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims; a numerical value indicating how many components to use in the computation of jackstraw plot from pca (default is number of pca components computed)
#' @return plot object
#' @export
seuratJackStrawPlot <- function(sceObject, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    return(JackStrawPlot(seuratObject, dims = 1:dims))
}

#' seuratComputeHeatmap
#' Computes the heatmap plot object from the pca slot in the input sce object
#' @param sceObject; sce object from which to compute heatmap (pca should be computed)
#' @param geneNamesSeurat; a list of rowames/genenames/featurenames for seurat object for consistent formatting
#' @param dims; numerical value indicating how many components to use for the computation of heatmap plot object (default is number of pca components computed) 
#' @return plot object
#' @export
seuratComputeHeatmap <- function(sceObject, geneNamesSeurat, dims) {
    seuratObject <- convertSCEToSeurat(sceObject, geneNamesSeurat)
    return(DimHeatmap(seuratObject, dims = 1:dims, fast = FALSE, combine = FALSE))
}

#' seuratHeatmapPlot
#' Modifies the heatmap plot object so it contains specified number of heatmaps in a single plot 
#' @param plotObject; plot object computed from seuratComputeHeatmap() function
#' @param dims; numerical value of how many heatmaps to draw (default is 0)
#' @param ncol; numerical value indicating that in how many columns should the heatmaps be distrbuted (default is 2)
#' @param labels; list() of labels to draw on heatmaps
#' @return modified plot object
#' @export
seuratHeatmapPlot <- function(plotObject, dims, ncol, labels) {
    componentsToPlot <- as.integer(gsub("[^0-9.]", "", labels))
    return(plot_grid(plotlist = plotObject[c(componentsToPlot)], ncol = ncol, labels = labels))
}

# ----
