#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

source("functions.R")

ui <- fluidPage(pageWithSidebar(
  headerPanel = headerPanel("Gene trees and coalescence"),
  
  sidebarPanel(
    sliderInput(
      inputId = "n",
      label = "Population size",
      value = 4,
      min = 2,
      max = 10,
      step = 1
    ),

    actionButton("goButton", "GO"),
    
    helpText(
      a(
        "This Shiny app was adapted from Silas Tittes's source code",
        target = "_blank",
        cex = 0.5,
        href = "https://github.com/silastittes/shiny_popgen"
      )
    )
    
  ),
  
  mainPanel =  mainPanel(plotOutput(outputId = 'viz'))
))


#back end code and response to user input
server <- function(input, output) {
  rand <- eventReactive(input$goButton, {
    return(list(n = input$n))
  })
  
  output$viz <- renderPlot({
    out <- rand()
    for_sim(n = out$n)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
