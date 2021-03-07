create_power_graph <- function(dat, alpha_level, beta_level){
  
  dat %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    filter(beta == beta_level) %>%
    mutate(m = paste("m =", m)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_brewer(palette = "Set1") +
    facet_grid(q ~ m, scales = "free_y") + 
    labs(x = "", y = "Power") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
  
}


create_power_graph(dat = power_dat, 
                   alpha_level = ".05",
                   beta_level = .5)


