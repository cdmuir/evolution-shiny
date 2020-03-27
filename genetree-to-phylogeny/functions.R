##Forward time simulation
simulate_treeseq <- function(N, L, u) {
  
  phy <- rcoal(N)
  
  phy$tip.label <- str_c("t", 1:N)
  data <- simSeq(
    phy, 
    l = L, 
    type = "USER", 
    levels = c("0", "1"), 
    rate = u,
    rootseq = rep("0", L), 
    ancestral = TRUE
  )
  
  list(phy = phy, data = data, N = N, L = L)
  
}

plot_treeseq <- function(.l) {
  
  phy <- .l$phy
  data <- .l$data %>%
    as.character()
  N <- .l$N
  L <- .l$L
  
  x1 <- phy$edge %>%
    as.data.frame() %>%
    set_colnames(c("parent", "node")) %>%
    mutate(
      parent = as.character(parent),
      node = str_c(ifelse(node <= N, "t", ""), node)
    )
  
  x2 <- data %>% 
    as.data.frame() %>%
    rownames_to_column("node") %>%
    mutate_if(is.factor, as.character) %>%
    pivot_longer(-node, names_to = "position") 
  
  x3 <- left_join(select(x1, node = parent), x2, by = "node") %>%
    rename(parent = node, parent_value = value)
  
  x4 <- x2 %>%
    full_join(x1, by = "node") %>%
    full_join(x3, by = c("parent", "position")) %>%
    distinct() %>%
    filter(!is.na(parent_value)) %>%
    mutate(
      label = node,
      node = as.integer(str_remove(node, "^t")),
      is_mutated = value != parent_value,
      mutation_position = ifelse(is_mutated, str_remove(position, "^V"), NA)
    ) %>%
    filter(!is.na(mutation_position)) %>%
    select(parent, node, is_mutated, mutation_position) %>%
    group_by(parent, node) %>%
    summarize(edge_label = str_c(mutation_position, collapse = ","))
  
  p <- ggtree(phy, ladderize = FALSE, size = 2) + 
    geom_rootedge() +
    ggtitle("Gene tree")
    
  p <- p %<+% 
    x4 + 
    geom_label(
      aes(x = branch, label = edge_label), 
      fill = "tomato"
    ) +
    scale_y_continuous(limits = c(0.5, N + 0.5)) +
    theme(
      legend.position = "none",
      text = element_text(size = 18)
    )
  
  x5 <- x2 %>%
    mutate(
      x = as.numeric(str_remove(position, "^V")),
      y = as.numeric(str_remove(node, "^t")),
      text = ifelse(value == 1, x, NA)
    ) %>%
    filter(node %in% phy$tip.label)
  
  q <- ggplot(x5, aes(x, y, fill = value, label = text)) +
    geom_tile(width = 1, height = 1, color = "black") +
    geom_text() +
    scale_fill_manual(values = c("steelblue", "tomato")) +
    scale_y_continuous(limits = c(0.5, N + 0.5)) +
    ggtitle("Sequence alignment") +
    theme_void() +
    theme(
      legend.position = "none"
    )
  
  r <- ggplot(data.frame(x = 1, y = 1:N, az = letters[1:N]), 
              aes(x, y, label = az)) +
    scale_y_continuous(limits = c(0.5, N + 0.5)) +
    geom_text() +
    theme_void()
  plot_grid(p, r, q, ncol = 3, rel_widths = c(1, 0.1, 1), align = "hv")
  
}

fit_treeseq <- function(.l) {
  
  pd <- subset(.l$data, subset = .l$phy$tip.label)
  dm <- dist.hamming(pd)
  phy <- upgma(dm)
  treedist <- RF.dist(.l$phy, phy, normalize = TRUE)
  phy$tip.label <- letters[as.integer(str_remove(phy$tip.label, "^t"))]
  
  ggtree(phy, ladderize = FALSE, size = 2) +
    geom_tiplab(offset = 0.01) +
    ggtitle("Inferred phylogeny", subtitle = bquote(paste("Distance between simulated and inferred tree: ", .(signif(treedist, 3)))))
  
}
