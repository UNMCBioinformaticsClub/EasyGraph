#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(shinycssloaders)

#sample data
sample_data <- readRDS("de_df_for_volcano.rds")
de <- sample_data[complete.cases(sample_data), ]


# Define UI for application that draws a histogram
ui <- navbarPage(

    title = "EasyGraph",
    
    
    tabPanel("Volcano",
        
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("files", label = "Upload DEG:", accept = ".fasta", multiple = TRUE),
            sliderInput("foldchange", label = "Fold Change", min = -1, max = 1, value = c(-.6,.6), step = .1),
            sliderInput("pvalue", "P-value", min = 0, max = .5, step = .01, value = 0.05)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("volcano") %>% withSpinner()
        )
    )
)
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$volcano <- renderPlot({
        data <- de
        
        top_genes <- de$gene_symbol[1:10]
        
        # add a column of NAs
        de$diffexpressed <- "NO"
        # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
        de$diffexpressed[de$log2FoldChange > input$foldchange[2] & de$pvalue < input$pvalue] <- "UP"
        # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
        de$diffexpressed[de$log2FoldChange < input$foldchange[1] & de$pvalue < input$pvalue] <- "DOWN"
        # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
        de$delabel <- NA
        de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
        
        
        
        de %>% 
            #mutate(Label = ifelse(de$gene_symbol %in% top_genes, de, "")) %>%  
            ggplot(aes(x = log2FoldChange, y = -log10(pvalue), label=delabel, col=diffexpressed)) + 
            geom_point(alpha=0.4) +
            theme_minimal() +
            geom_text_repel() +
            scale_color_manual(values=c("blue", "black", "red")) +
            geom_vline(xintercept=c(input$foldchange[1], input$foldchange[2]), col="red") +
            geom_hline(yintercept=-log10(input$pvalue), col="red")
           
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
