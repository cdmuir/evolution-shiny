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

    helpText("This is a haploid population, so N is the number of allele copies."),
    
    actionButton("goButton", "GO"),
    
    helpText(
      a(
        "Source code",
        target = "_blank",
        cex = 0.5,
        href = "https://github.com/cdmuir/evolution-shiny"
      )
    ),
    
    helpText(
      a(
        "This Shiny app was adapted from Silas Tittes's source code",
        target = "_blank",
        cex = 0.5,
        href = "https://github.com/silastittes/shiny_popgen"
      )
    )
    
  ),
  
  mainPanel =  mainPanel(
    p("The input to the left allow you to select a different haploid population size to simulate a single gene tree. Blue points are the present-day allele copies and their descendents. Orange lines show relationships between ancesotrs and descendents. Grey points represent lineages that have no descendents in the present-day population. Each time you hit the \"GO\" button you'll get a different random simulation."),
    
    plotOutput(outputId = 'viz'),
    
    h2("Suggested activties:"),
    
    p("1. Simulate a tree with a small population size (e.g. N = 4). On pen and paper, try drawing the phylogenetic tree we would get from this simulation."),
    
    p("2. How many generations does it take for coalesence, where all surviving allele copies share a common ancestor? You'll need to count the number of columns in the plot above. Try several simulations with a small population (N=4) and a larger population (N=10). How does the coalesence time change?")
    )
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
