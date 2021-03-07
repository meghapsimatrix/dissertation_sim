create_power_graph <- function(dat, cov, alpha){
  
  dat %>%
    filter(str_detect(cov_test, cov)) %>%
    filter(alpha == alpha) %>%
    filter(test != "Naive-F") %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    ggplot(aes(x = m, y = rej_rate, fill = test)) + 
    #geom_hline(yintercept = .8, linetype = "solid") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = paste("Power", cov), fill = "") + 
    theme_bw() +
    theme(legend.position = "bottom",
          plot.caption=element_text(hjust = 0, size = 10))
  
}

create_power_graph(dat = power_dat, cov = "X1", alpha = .05)
ggsave("sim_results/graphs/power_cov/power_X1.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_graph(dat = power_dat, cov = "X2", alpha = .05)
ggsave("sim_results/graphs/power_cov/power_X2.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_graph(dat = power_dat, cov = "X3", alpha = .05)
ggsave("sim_results/graphs/power_cov/power_X3.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_graph(dat = power_dat, cov = "X4", alpha = .05)
ggsave("sim_results/graphs/power_cov/power_X4.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_graph(dat = power_dat, cov = "X5", alpha = .05)
ggsave("sim_results/graphs/power_cov/power_X5.png", device = "png", dpi = 500, height = 7, width = 12)




create_power_rat_graph <- function(dat, cov, alpha_level){
  
  
  dat %>%
    filter(str_detect(cov_test, cov)) %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    ggplot(aes(x = m, y = power_ratio, fill = m)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = paste("Power Ratio: HTZ/ CWB", cov)) + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
  
}

create_power_rat_graph(power_ratio, cov = "X1", alpha_level = ".05")
ggsave("sim_results/graphs/power_cov/power_rat_X1.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph(power_ratio, cov = "X2", alpha_level = ".05")
ggsave("sim_results/graphs/power_cov/power_rat_X2.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_rat_graph(power_ratio, cov = "X3", alpha_level = ".05")
ggsave("sim_results/graphs/power_cov/power_rat_X3.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_rat_graph(power_ratio, cov = "X4", alpha_level = ".05")
ggsave("sim_results/graphs/power_cov/power_rat_X4.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph(power_ratio, cov = "X5", alpha_level = ".05")
ggsave("sim_results/graphs/power_cov/power_rat_X5.png", device = "png", dpi = 500, height = 7, width = 12)



