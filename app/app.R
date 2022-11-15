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
library(shinydashboard)
library(shinyWidgets)
library(fontawesome)
library(DT)


#sample rna-seq
sample_data <- readRDS("de_df_for_volcano.rds")
de <- sample_data[complete.cases(sample_data), ]

#sample chip-seq
#chipseq <- read.csv("E6vWT_log2FC.csv")


# Define UI for application that draws a histogram
ui <- navbarPage(

    title = "EasyGraph",
    
    header = tagList(
        useShinydashboard()
    ),
    
    
    
    tabPanel("Volcano",
             
             tabsetPanel(
                 tabPanel("RNA-seq",
                          
                          # Sidebar with a slider input for number of bins 
                          sidebarLayout(
                              sidebarPanel(
                                  fileInput("files", label = "Upload DEG:", accept = ".csv"),
                                  sliderInput("foldchange", label = "Fold Change", min = -1, max = 1, value = c(-.6,.6), step = .1),
                                  sliderInput("pvalue", "P-value", min = 0, max = .5, step = .01, value = 0.05)
                              ),
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                  plotOutput("volcano") %>% withSpinner()
                              )
                          )
                 ),
                 tabPanel("ChIP-seq",
                          
                          # Sidebar with a slider input for number of bins 
                          sidebarLayout(
                              sidebarPanel(
                                  fileInput("files2", label = "Upload mannorm:", accept = ".csv"),
                                  sliderInput("foldchange2", label = "Fold Change", min = -8, max = 8, value = c(-1,1), step = .5),
                                  sliderInput("pvalue2", "P-value", min = 0, max = .5, step = .01, value = 0.05)
                              ),
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                  plotOutput("volcano2") %>% withSpinner(),
                                  infoBoxOutput("peaksLostBox"),
                                  infoBoxOutput("peaksGainedBox"),
                                  dataTableOutput("diffpeaks")
                              )
                          )
                 )
             )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    
    rnaInput <- reactive({
        if(is.null(input$files)){
            rnaData <- de
        }else{
            data <- input$files
            rnaData <- read.csv(data$datapath)
        }
        return(rnaData)
    })

    chipInput <- reactive({
        if(is.null(input$files2)){
            #chipData <- chipseq
            validate(need(input$files2, "Please select manorm output csv"))
        }else{
            data <- input$files2
            chipData <- read.csv(data$datapath)
            chipData$Log2FC <- log2(chipseq[,9]/chipseq[,10])
        }
        return(chipData)
    })



    output$volcano <- renderPlot({
        
        data <- rnaInput()
        
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
    
    
    DEpeaks <- reactive({
        de <- chipInput()
    
        # add a column of NAs
        de$diffpeaks <- "NO"
        # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
        de$diffpeaks[de$Log2FC > input$foldchange2[2] & de$P_value < input$pvalue2 ] <- "UP"
        # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
        de$diffpeaks[de$Log2FC < input$foldchange2[1] & de$P_value < input$pvalue2] <- "DOWN"
        # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
        de$delabel <- NA
        de$peak <- paste(de$chr, de$summit, sep = "_")
        de$delabel[de$diffpeaks != "NO"] <- de$peak[de$diffpeaks != "NO"]
        
        
        up <- sum(de$diffpeaks=="UP")
        down <- sum(de$diffpeaks=="DOWN")
        
        list(de=de, up=up, down=down)
    })
    
    output$diffpeaks <- renderDataTable({
        de <- DEpeaks()$de %>% filter(diffpeaks != "NO")
        de <- de[,0:10]
        
        datatable(de,
                  rownames = FALSE,
                  extensions = 'Buttons',
                  options = list(
                      scrollX = TRUE,
                      searching = TRUE,
                      pageLength = 5,
                      dom = 'Blfrtip',
                      buttons = list(
                          list(extend = 'csv', fieldBoundary = ''),
                          'excel'
                      )
                  )
                 )
        }, server = F)
    
    output$peaksLostBox <- renderInfoBox({
        infoBox(
            "Peaks Lost", DEpeaks()$down, icon = icon("fa-solid fa-down"),
            color = "blue"
        )
    })
    
    
    output$peaksGainedBox <- renderInfoBox({
        infoBox(
            "Peaks Gained", DEpeaks()$up, icon = icon("fa-solid fa-up"),
            color = "red"
        )
    })
    
    
    output$volcano2 <- renderPlot({
        
        de <- DEpeaks()$de
        
        de %>%
            ggplot(aes(x = Log2FC, y = -log10(P_value), label=delabel, col=diffpeaks)) +
            geom_point(alpha=0.4) +
            theme_minimal() +
            scale_color_manual(values=c("blue", "black", "red")) +
            geom_vline(xintercept=c(input$foldchange2[1], input$foldchange2[2]), col="red") +
            geom_hline(yintercept=-log10(input$pvalue2), col="red")
        
        
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
