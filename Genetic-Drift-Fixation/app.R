library(dplyr)
library(shiny)
library(plotly)
library(stringr)
library(tidyr)

## adjustable via sliders
##
initial_p <- 0.5
initial_N <- 10
initial_n_pops <- 5

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  ## Application title
  ##
  titlePanel("Allele frequency changes with genetic drift"),
  
  ## selection initial allele frequency, population size, and number of
  ## generations for simulation
  ##
  sidebarLayout(
    sidebarPanel(
      numericInput("p",
                  "Initial frequency of A",
                  min = 0.0,
                  max = 1.0,
                  value = initial_p,
                  step = 0.01),
      helpText("Frequency must be between 0 and 1"),
      numericInput("N",
                  "Population size",
                  min = 1,
                  max = 100,
                  value = initial_N),
      helpText("N must be between 1 and 100. The plot will take several seconds to load when N is high."),
      sliderInput("n_pops",
                  "Number of populations",
                  min = 1,
                  max = 9,
                  value = initial_n_pops),
      actionButton("go", "Simulate Genetic Drift!", icon("refresh")),
      helpText("When you click the button above, the figure to the right will update with new simulations using the parameters you entered above."),
    ),
    
    ## Show a plot of the simulated change in allele frequencies
    ##
    mainPanel(
      # p("Notes explaining the principles of genetic drift are available at:",
      #   uiOutput("darwin")),
      p("The inputs to the left allow you to select a different starting initial allele frequency, diploid population size (so the number of gametes is 2N, and number of populations. Each line represents the history of allele frequency change in one population. All populations begin with an identical allele frequency. Each time you change one of the inputs, you'll get a new set of simulation results after you hit the \"Simulate Genetic Drift!\" button. If you hit \"Play\" without changing the input, you'll get a duplicate of the plot you just saw."),
      plotlyOutput("allele_frequency_plot"),
      hr(),
      helpText("To complete the homework, copy-and-paste the lines below into the Google Form:"), uiOutput("google"),
      helpText("Copy this set of numbers into 'Final allele frequency' prompt"),
      verbatimTextOutput("copy1"),
      helpText("Copy this set of numbers into 'Time to fixation' prompt"),
      verbatimTextOutput("copy2"),
      hr(),
      helpText("Here's a nicer looking table describing the output"),
      tableOutput("summary_table"),
      hr(),
      p("This Shiny app was adapted from Kent Holsinger's source code:",
        uiOutput("github"))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  url_1 <- a("https://forms.gle/JXmRExuNmo42HRQm6",
             href="https://forms.gle/JXmRExuNmo42HRQm6")
  output$google <- renderUI({
    tagList("", url_1)
  })
  url_2 <- a("https://kholsinger.github.io/PopGen-Shiny/",
             href="https://kholsinger.github.io/PopGen-Shiny/")
  output$github <- renderUI({
    tagList("", url_2)
  })
  
  simulation <- eventReactive(input$go, {
    ## generate allele frequencies
    ##
    N <- round(input$N)
    n_gen <- 20 * N
    two_n <- 2 * N
    
    df <- data.frame(
      replicate(n = input$n_pops, expr = {
        sims <- rep(NA, n_gen)
        sims[1] <- input$p
        i <- 1
        while (sims[i] > 0 & sims[i] < 1) {
          i <- i + 1
          sims[i] <- rbinom(n = 1, size = two_n, prob = sims[i-1]) /
            two_n
        }
        return(sims)
      })
    ) %>% 
      filter_all(any_vars(!is.na(.))) %>%
      mutate(t = 0:(nrow(.) - 1)) %>% 
      gather(pop, p, -t) 
    
    df <- left_join(df, df %>%
                      filter(!is.na(p)) %>%
                      group_by(pop) %>%
                      summarize(p_final = last(p)), by = "pop"
    ) %>%
      mutate(p = ifelse(is.na(p), p_final, p))
    
    ## construct data frame for plot
    ##
    n_gen <- max(df$t)
    d <- crossing(
      df,
      frame = 0:n_gen
    ) %>%
      filter(frame >= t)
    
    d
    
  })
  
  output$allele_frequency_plot <- renderPlotly({

    d <- simulation()
    n_gen <- max(d$t)
    
    ## plot it
    ##
    if (n_gen > 100) {
      x <- round(seq(0, n_gen, length.out = 100))
    } else {
      x <- 0:n_gen
    }
    
    gp <- ggplot(filter(d, t %in% x), aes(x = t, y = p, color = pop, frame = frame)) +
      geom_line() +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_brewer(palette = "Set1") + 
      xlab("Time (generations)") +
      ylab("Allele Frequency, p") +
      theme_bw() +
      theme(legend.position = "none")
    
    p_plot <- gp %>%
      ggplotly() %>%
      animation_opts(
        frame = 1000 / n_gen,
        transition = 0,
        redraw = FALSE
      ) %>%
      animation_slider(
        hide = TRUE
      )
  })
  
  output$summary_table <- renderTable({
    d <- simulation()
    d %>%
      filter(frame == max(frame)) %>%
      group_by(pop) %>%
      summarize(
        `Initial allele frequency` = first(p),
        `Final allele frequency` = last(p),
        `Time to fixation` = min(which(p %in% c(0, 1)))
      ) %>%
      rename(Population = pop) %>%
      mutate(N = input$N, Population = 1:nrow(.)) %>%
      select(Population, N, `Initial allele frequency`, `Final allele frequency`, `Time to fixation`)
  })

  output$copy1 <- renderPrint({
    
    simulation() %>%
      filter(frame == max(frame)) %>%
      group_by(pop) %>%
      summarize(p = last(p)) %>% 
      pull(p) %>%
      str_c(collapse = ",") %>%
      cat()
    
  })

  output$copy2 <- renderPrint({
    
    simulation() %>%
      filter(frame == max(frame)) %>%
      group_by(pop) %>%
      summarize(t = min(which(p %in% c(0, 1)))) %>% 
      pull(t) %>%
      str_c(collapse = ",") %>%
      cat()
    
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)

