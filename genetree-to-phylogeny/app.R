#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(BiocManager)
options(repos = BiocManager::repositories())

library(ape)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggtree)
library(magrittr)
library(phangorn)
library(stringr)
library(tibble)
library(tidyr)

source("functions.R")

ui <- fluidPage(pageWithSidebar(
  headerPanel = headerPanel("Coalescent simulations with mutations and alignments"),
  
  sidebarPanel(
    sliderInput(
      inputId = "N",
      label = "Population size",
      value = 10,
      min = 2,
      max = 20,
      step = 1
    ),
    
    helpText("This is a haploid population, so N is the number of allele copies."),
    
    sliderInput(
      inputId = "L",
      label = "Sequence length",
      value = 10,
      min = 1,
      max = 20,
      step = 1
    ),
    
    helpText("The length of the DNA sequence to simulate. Each position starts as blue in the ancestor, then can mutate to orange. Occasionally, you'll see reverse mutations back to blue if you look carefully."),
    
    sliderInput(
      inputId = "u",
      label = "Mutation rate",
      value = 0.1,
      min = 0.01,
      max = 1,
      step = 0.01
    ),
    
    helpText("Mutation rate"),
    
    actionButton("goButton", "GO"),
    
    helpText(
      a(
        "Source code",
        target = "_blank",
        cex = 0.5,
        href = "https://github.com/cdmuir/evolution-shiny"
      )
    )
    
  ),
  
  mainPanel =  mainPanel(
    h5("Gene tree"),
    p("The input to the left allow you to select a different haploid population size to simulate a single gene tree that can accumulate mutations. The phylogeny (gene tree) shows the actual (simulated) relationship among allele copies."),

    h5("Mutations"),
    p("Mutations can occur on any branch at any time at random. Once they occur, they're inherited from parent to offspring. The orange boxes above the branches show where mutations occured and the number indicates the base position of that mutation. You can adjust the mutation rate in the input."),
    
    h5("Alignment"),
    p("On the right side is an alignment of DNA sequences with length L, which you can adjust in the input panel. For simplicity, there are only two states, blue and orange. The ancestor is all blue. Mutations change the base from blue to orange, but reverse mutations can occasionally change them back."),
    
    p("Each time you hit the \"GO\" button you'll get a different random simulation."),
    
    plotOutput(outputId = 'viz1'),
    
    h2("Some things to make sure you understand (ask if you don't!)"),
    
    p("1. Phylogenies can represent the evolutionary relationships among many different entities. What entities (species, alleles, individuals) are being simulated here?"),
    
    p("2. The orange boxes above the branches show where mutations occured in the past. Do you understand what the number means? Given the branch where that mutation occurred, which tips should share that mutation? Did the alignment match your expectations?"),
    
    p("3. The alignment shows a simulated sequence of DNA, one per tip. How many sites are in this alignment? What states can each site have in this simulation?"),
    
    h2("Suggested activity: Compare simulated and inferred phylogenies"),

    helpText("The simulated tree above is the real relationship; the tree below is the inferred set of relationships using the simulated sequence data in the alignment above. How do the actual relationships compare to those inferred using sequence data?"),
    
    helpText("The subtitle displays a standardized distance between real and inferred trees. Low values mean that the inferred tree is very close to the real tree; higher values means that it is farther away. In other words, lower distance indicates a more accurate tree. Try several simulations with different values of N, L, and mutation rate. How do these parameters effect the accuracy of phylogenetic inference?"),
    
    helpText("Does greater N increase or decrease the accuracy?"),
    
    helpText("Does greater L increase or decrease the accuracy?"),
    
    helpText("Does greater mutation rate increase or decrease the accuracy?"),
    
    plotOutput(outputId = 'viz2')
    
)))


#back end code and response to user input
server <- function(input, output) {
  rand <- eventReactive(input$goButton, {
    # return(list(N = input$N, L = input$L, u = input$u))
    simulate_treeseq(N = input$N, L = input$L, u = input$u)
  })
  
  output$viz1 <- renderPlot({
    # out <- rand()
    # s <- simulate_treeseq(N = out$N, L = out$L, u = out$u)
    s <- rand()
    plot_treeseq(s)
  })
  
  output$viz2 <- renderPlot({
    s <- rand()
    fit_treeseq(s)
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
