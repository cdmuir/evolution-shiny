library(dplyr)
library(ggplot2)
library(shiny)
library(cowplot)
library(plotly)
library(tidyr)

rm(list=ls())

## adjustable via sliders
##
initial_p <- 0.01
initial_s <- 0.1
initial_h <- 0.5
initial_n_gen <- 250

## dynamics of viability selection in one population
##
calc_wbar <- function(p, s, h) {
  
  # Equation in class
  # p ^ 2 + 2 * p * (1 - p) * (1 - s * h) + (1 - p) ^ 2 * (1 - s)
  
  # More compact equation (identical)
  1 - 2 * p * (1 - p) * s * h - (1 - p) ^ 2 * s
  
}

calc_deltap <- function(p, s, h, wbar) {
  
  # Equation in class (for h = 0.5)
  # wbar <- calc_wbar(p, s, h)
  # w_11 <- 1 / wbar
  # w_12 <- (1 - h * s) / wbar
  # p ^ 2 * w_11 + p * (1 - p) * w_12 - p
  
  # equation for any h
  p * (1 - p) * (p * h * s + (1 - p) * s * (1 - h)) / wbar
  
}

# Variance in fitness for one-allele, two-locus model
calc_varw <- function(p, s, h, wbar) {
  
  p ^ 2 * (1 - wbar) ^ 2 + 
    2 * p * (1 - p) * ((1 - s * h) - wbar) ^ 2 +
    (1 - p) ^ 2 * ((1 - s) - wbar) ^ 2 
  
}

simulate_selection <- function(p, s, h, n_gen) {
  
  wbar <- calc_wbar(p, s, h)
  deltap <- calc_deltap(p, s, h, wbar)
  varw <- calc_varw(p, s, h, wbar)
  
  ret <- data.frame(
    t = 0:n_gen,
    p = c(p, numeric(n_gen)),
    wbar = c(wbar, numeric(n_gen)),
    deltap = c(deltap, numeric(n_gen)),
    varw = c(varw, numeric(n_gen))
  )
  
  for (i in 1:n_gen) {
    
    p <- p + deltap
    wbar <- calc_wbar(p, s, h)
    deltap <- calc_deltap(p, s, h, wbar)
    varw <- calc_varw(p, s, h, wbar)
    ret$p[i + 1] <- p
    ret$wbar[i + 1] <- wbar
    ret$deltap[i + 1] <- deltap
    ret$varw[i + 1] <- varw
    
  }
  
  ret
  
}

## Define UI
##
ui <- fluidPage(
  titlePanel("Selection at one locus with two alleles"),
  sidebarLayout(
    sidebarPanel(
      numericInput("p",
                  "Initial allele frequency (p)",
                  min = 0.0,
                  max = 1.0,
                  value = initial_p),
      numericInput("s",
                  "Selection coefficient (s)",
                  min = 0.0,
                  max = 1.0,
                  value = initial_s),
      helpText("Higher s means stronger selection."),
      numericInput("h",
                  "Dominance coefficient (h)",
                  min = 0.0,
                  max = 1.0,
                  value = initial_h),
      helpText("When h = 0, the beneficial allele is completely dominant."),
      numericInput("n_gen",
                  "Number of generations",
                  min = 1,
                  max = 10000,
                  value = initial_n_gen),
      actionButton("go", "Simulate Selection!", icon("play")),
    ),
    mainPanel(
      h4("Allele frequency dynamics"),
      helpText("After setting parameters, hit \"Simulate Selection!\" to run the simulation."),
      helpText("Once simulation is done, an empty plot will appear below. Hit \"Play\" to see how the beneficial allele frequency changes over time."),
      plotlyOutput("dynamic_plot1"),
      helpText("You can use slider to check specific timepoint."),
      hr(),
      h4("Evolutionary rate is proportional to the variance in fitness"),
      helpText("Hit \"Play\" to see how the rate of evolution (\"Change in p\") changes along with the variance in fitness (\"Var(W)\")."),
      plotlyOutput("dynamic_plot2"),
      hr(),
      p("This Shiny app was adapted from Kent Holsinger's source code:"),
      uiOutput("github")
    )
  )
)

## Define server logic
##
server <- function(input, output) {
  url_1 <- a("http://darwin.eeb.uconn.edu/eeb348/lecture-notes/selection.pdf",
           href="http://darwin.eeb.uconn.edu/eeb348/lecture-notes/selection.pdf")
  output$darwin <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
           href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })

  simulate <- eventReactive(input$go, {
    df <- simulate_selection(input$p, input$s, input$h, input$n_gen)
    df
  })
  
  output$dynamic_plot1 <- renderPlotly({
    
    ## generate allele frequencies
    ##
    df <- simulate()
    n_gen <- max(df$t)
    
    max_deltap <- max(df$deltap)
    
    ## construct data frame for plot
    ##
    d <- df %>%
      mutate(deltap = deltap / max_deltap, varw = varw / max_deltap) %>%
      rename(`Average fitness` = wbar) %>%
      pivot_longer(-t) %>%
      mutate(size = ifelse(name == "p", 1.2, 1)) %>%
      crossing(frame = 0:n_gen) %>%
      filter(frame >= t) %>%
      filter(name %in% c("p", "Average fitness")) %>%
      mutate(name = factor(name, levels = c("p", "Average fitness")))
    
    if (n_gen > 100) {
      x <- round(seq(0, n_gen, length.out = 100))
    } else {
      x <- 0:n_gen
    }
    
    gp <- ggplot(filter(d, frame %in% x), aes(x = t, y = value, color = name, frame = frame)) +
      geom_line(size = 2) +
      scale_y_continuous(limits = c(0, 1), sec.axis = sec_axis(~ . * max_deltap)) +
      scale_color_brewer("Variable", palette = "Set1") +
      xlab("Time (generations)") +
      ylab("Allele Frequency, p") +
      theme_cowplot()
    
    ## plot it
    ##
    p_plot <- gp %>%
      ggplotly() %>%
      animation_opts(
        frame = 1000 / length(x),
        transition = 0,
        redraw = FALSE
      ) %>%
      animation_slider(
        hide = FALSE
      )
  })

  output$dynamic_plot2 <- renderPlotly({
    
    ## generate allele frequencies
    ##
    df <- simulate()
    n_gen <- max(df$t)
    
    ## construct data frame for plot
    ##
    d <- df %>%
      rename(`Change in p` = deltap, `Var(W)` = varw) %>%
      pivot_longer(-t) %>%
      crossing(frame = 0:n_gen) %>%
      filter(frame >= t) %>%
      filter(name %in% c("Change in p", "Var(W)")) %>%
      mutate(name = factor(name, levels = c("Change in p", "Var(W)")))
    
    if (n_gen > 100) {
      x <- round(seq(0, n_gen, length.out = 100))
    } else {
      x <- 0:n_gen
    }
    
    gp <- ggplot(filter(d, frame %in% x), aes(x = t, y = value, color = name, frame = frame)) +
      geom_line(size = 2) +
      scale_color_brewer("Variable", palette = "Set2") +
      xlab("Time (generations)") +
      ylab("Value") +
      theme_cowplot()
    
    ## plot it
    ##
    p_plot <- gp %>%
      ggplotly() %>%
      animation_opts(
        frame = 1000 / length(x),
        transition = 0,
        redraw = FALSE
      ) %>%
      animation_slider(
        hide = FALSE
      )
  })
  
}

## Run the application
##
shinyApp(ui = ui, server = server)

