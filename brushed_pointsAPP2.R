  ## Manually removing some points with Shiny APP
  library(ggplot2)
  library(Cairo)   # For nicer ggplot2 output when deployed on Linux
  library(shiny)
  
  ui <- basicPage(
    plotOutput("plot1", brush = "plot_brush"),
    verbatimTextOutput("info"),
    actionButton("save", "Save")
  )
  
  server <- function(input, output) {
    output$plot1 <- renderPlot({
  p.sd.ADGC.WASHU + geom_point(data = HISPANIC.SUBSET[c("PC1","PC2")], aes(col="NEW_Hispanic")) +
  scale_color_manual(values = c(ADGC = 'black', CEU='red', JPT = 'green', Reported_NHW = 'yellow', Reported_Hispanic = 'grey', Reported_Asian = 'violet', Reported_African = 'orange', YRI = "blue", WashU = "pink", NEW_Hispanic = "gold")) 
    })
    
    data <- reactive({
      brushedPoints(HISPANIC.SUBSET, input$plot_brush, xvar = "PC1", yvar = "PC2")
    })
    
    output$info <- renderPrint({data()})
    observeEvent(input$save, {
    write.table(data(), "brushed_data.csv", sep = "\t", col.names = !file.exists("brushed_data.csv"), append = T)
      
    })
  }
  
  shinyApp(ui, server)
