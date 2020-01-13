library(shiny);
source("Seurat.R");

# Main User Inteface ---
ui <- navbarPage("SingleCellTK", theme = "bootstrap.css", id = "mainPage",
        navbarMenu("Workflows",
            tabPanel("Seurat",
                fluidPage(
                    Seurat_UI(id = "id_1")
                )),
            tabPanel("OSCA")));
# ----

# Main Server ---
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)
server <- function(input, output, session) {
    callModule(module = Seurat_Server, id = "id_1", session)
}
# ----

shinyApp(ui, server);