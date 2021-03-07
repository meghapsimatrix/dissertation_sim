cov <- "X1"

make_data <- function(cov){
  
  power_dat %>%
    filter(test != "Naive-F") %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    filter(str_detect(cov_test, cov)) %>%
    mutate(cov = cov)
  
}

covs <- c("X1", "X2", "X3", "X4", "X5")

all_dat <- map_dfr(covs, make_data)

create_power_covs <- function(beta_level, alpha_level, dat = all_dat){
  
  dat %>%
    filter(alpha == alpha_level) %>%
    filter(beta == beta_level) %>%
    ggplot(aes(x = m, y = rej_rate, fill = test)) + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(cov ~ q, scales = "free_y") + 
    labs(x = "Number of Studies", y = "Power", fill = "") +
    ggtitle(bquote(beta == .(beta_level))) + 
    theme_bw() +
    theme(legend.position = "bottom",
          plot.caption=element_text(hjust = 0, size = 10))
  
}


b5 <- create_power_covs(beta_level = .5, alpha_level = ".05")
b1 <- create_power_covs(beta_level = .1, alpha_level = ".05")

b1 + b5 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs/power_05_covs.png", device = "png", dpi = 500, height = 7, width = 12)


b5 <- create_power_covs(beta_level = .5, alpha_level = ".01")
b1 <- create_power_covs(beta_level = .1, alpha_level = ".01")

b1 + b5 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs/power_01_covs.png", device = "png", dpi = 500, height = 7, width = 12)


b5 <- create_power_covs(beta_level = .5, alpha_level = ".10")
b1 <- create_power_covs(beta_level = .1, alpha_level = ".10")

b1 + b5 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs/power_10_covs.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_covs <- function(alpha_level, dat = all_dat){
  
  dat %>%
    filter(alpha == alpha_level) %>%
    filter(q == "q = 1") %>%
    #filter(beta == beta_level) %>%
    ggplot(aes(x = m, y = rej_rate, fill = test)) + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(beta ~ cov, scales = "free_y", labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power", fill = "") +
    #ggtitle(bquote(beta == .(beta_level))) + 
    theme_bw() +
    theme(legend.position = "bottom",
          plot.caption=element_text(hjust = 0, size = 10))
  
}

create_power_covs(alpha_level = ".05")
ggsave("sim_results/graphs/power_05_covs_q1.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_rat_graph_covs <- function(alpha_level, dat = power_ratio_all){
  
  
  dat %>%
    filter(alpha == alpha_level) %>%
    filter(q == "q = 1") %>%
    # filter(beta == beta_level) %>%
    ggplot(aes(x = m, y = power_ratio, fill = m)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(beta ~ cov, scales = "free_y", labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/ CWB") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
  
}

create_power_rat_graph_covs(alpha_level = ".05")
ggsave("sim_results/graphs/power_rat_05_covs_q1.png", device = "png", dpi = 500, height = 7, width = 12)


power_ratio_all <- all_dat %>%
  filter(test %in% c("CWB", "HTZ")) %>%
  select(-c(starts_with("mcse"))) %>%
  spread(test, rej_rate) %>%
  mutate(power_diff = CWB - HTZ,
         power_ratio = HTZ/CWB) %>%
  group_by(m, rho, tau, alpha, q, beta, cov_test, cov) %>%
  summarize_at(vars(power_ratio), mean)

power_ratio %>%
  filter(power_ratio > 1, m == 10) %>%
  View()


create_power_rat_graph_covs <- function(beta_level, alpha_level, dat = power_ratio_all){
  
  
  dat %>%
    filter(alpha == alpha_level) %>%
    filter(beta == beta_level) %>%
    ggplot(aes(x = m, y = power_ratio, fill = m)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(cov ~ q) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/ CWB") + 
    ggtitle(bquote(beta == .(beta_level))) + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
  
}

b5 <- create_power_rat_graph_covs(beta_level = .5, alpha_level = ".05")
b1 <- create_power_rat_graph_covs(beta_level = .1, alpha_level = ".05")

b1 + b5 

ggsave("sim_results/graphs/power_rat_05_covs.png", device = "png", dpi = 500, height = 7, width = 12)


b5 <- create_power_rat_graph_covs(beta_level = .5, alpha_level = ".01")
b1 <- create_power_rat_graph_covs(beta_level = .1, alpha_level = ".01")

b1 + b5 

ggsave("sim_results/graphs/power_rat_01_covs.png", device = "png", dpi = 500, height = 7, width = 12)


b5 <- create_power_rat_graph_covs(beta_level = .5, alpha_level = ".10")
b1 <- create_power_rat_graph_covs(beta_level = .1, alpha_level = ".10")

b1 + b5 

ggsave("sim_results/graphs/power_rat_10_covs.png", device = "png", dpi = 500, height = 7, width = 12)


