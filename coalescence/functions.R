##Forward time simulation
for_sim <- function(n) {
  
  # n <- 10
  gen <- n * 10

  samps <- 1:n
  sim_forward <- replicate(n = gen - 1, sort(sample(samps, replace = TRUE)))
  
  # Figure out coalescence time ----
  connect_gen1 <- rep(gen:2, each = n)
  connect_gen2 <- connect_gen1 - 1
  
  samps <- sample(1:n, size = n, replace = FALSE)
  t0 <- samps
  c_gen <- gen - 1 
  
  while (c_gen > 0 & length(unique(t0)) > 1) {
    t0 <- sim_forward[t0, c_gen]
    c_gen <- c_gen - 1 # coal is gen - c_gen
  }
  
  # Re-run if no coalesence ----
  while (c_gen == 0) {
    samps <- 1:(n)
    sim_forward <- replicate(n = gen - 1, sort(sample(samps, replace = TRUE)))
    
    # Figure out coalescence time ----
    connect_gen1 <- rep(gen:2, each = n)
    connect_gen2 <- connect_gen1 - 1
    
    samps <- sample(1:n, size = n, replace = FALSE)
    t0 <- samps
    c_gen <- gen - 1 
    
    while (c_gen > 0 & length(unique(t0)) > 1) {
      t0 <- sim_forward[t0, c_gen]
      c_gen <- c_gen - 1 # coal is gen - c_gen
    }
  }
  
  df_balls <- 1:gen %>% 
    map_dfr( ~ {
      data.frame(x = 1:n, y = rep(.x, n))
      }) %>%
    filter(y <= gen - c_gen)
  
  df_segments <- bind_cols(
    data.frame(x = c(sim_forward), y = connect_gen1),
    data.frame(xend = samps, yend = connect_gen2)
  ) %>%
    filter(y <= gen - c_gen)
  
  df_segments1 <- data.frame(
    x = samps, y = 1, 
    xend = sim_forward[samps, gen - 1], yend = 2
  )

  group_gen <- unique(sim_forward[samps, gen - 1])
  df_balls1 <- data.frame(
    x = c(samps, group_gen), y = rep(1:2, c(length(samps), length(group_gen))) 
  )
  
  samps <- sim_forward[samps, gen - 2]

  for(i in 1:(gen - c_gen - 2)){
    
    df_segments1 <- df_segments1 %>%
      bind_rows(
        data.frame(
          x = group_gen, y = i + 1,
          xend = sim_forward[group_gen, gen - (i + 1)], yend = i + 2
        )
      )
    
    df_balls1 <- df_balls1 %>%
      bind_rows(
        data.frame(
          x = c(group_gen, sim_forward[group_gen, gen - (i + 1)]),
          y = rep(
            i + 1:2, 
            c(length(group_gen), length(sim_forward[group_gen, gen - (i + 1)]))
          )
        )
      )
    group_gen <- unique(sim_forward[group_gen, gen - (i + 1)])

  }
  
    
  gp <- ggplot(df_balls, aes(x, y)) +
    geom_segment(
      data = df_segments,
      aes(xend = xend, yend = yend),
      color = "grey"
    ) +
    geom_segment(
      data = df_segments1,
      aes(xend = xend, yend = yend),
      color = "tomato", size = 1.5,
    ) +
    geom_point(color = "grey") +
    geom_point(
      data = data.frame(x = unique(t0), y = gen - c_gen),
      color = "steelblue", size = 5
    ) +
    geom_point(
      data = df_balls1,
      color = "steelblue", size = 2
    ) +
    ylab("Time (generations)") +
    coord_flip(ylim = c(gen - c_gen, 1)) +
    theme_minimal() +
    theme(
      axis.line.x = element_line(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_blank()
    )
  
  gp

}
