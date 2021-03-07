# Power X1 5 --------------------------------------------------------------

create_power_graph <- function(dat, beta, cov, alpha){
  
  dat %>%
    filter(beta_type == beta, alpha == alpha) %>%
    filter(str_detect(cov_test, cov)) %>%
    filter(test != "Naive-F") %>%
    mutate(m = paste("m =", m))  %>%
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_hline(yintercept = .8, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .1)) +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Method", y = paste("Power", beta, cov)) + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
}


# The first binary covariate, X1, is a study level covariate with large imbalance, 
# equaling 1 in 15% of the studies.

create_power_graph(dat = power_dat, beta = "B5", cov = "X1", alpha = .05)


# Power X1 1 --------------------------------------------------------------

create_power_graph(dat = power_dat, beta = "B1", cov = "X1", alpha = .05)


# Power X2 5 --------------------------------------------------------------
# The second binary covariate, X2, is a covariate that varies within studies, 
# equaling 1 in 10% of the effect size estimates overall and in 0 to 20% 
# of the effect size estimates within a study.

create_power_graph(dat = power_dat, beta = "C5", cov = "X2", alpha = .05)


# Power X2 1 --------------------------------------------------------------

create_power_graph(dat = power_dat, beta = "C1", cov = "X2", alpha = .05)


# Power X3 5 --------------------------------------------------------------

# X3 is a normally distributed study level covariate

create_power_graph(dat = power_dat, beta = "D5", cov = "X3", alpha = .05)


# Power X3 1 --------------------------------------------------------------

create_power_graph(dat = power_dat, beta = "D1", cov = "X3", alpha = .05)


# Power X4 5 --------------------------------------------------------------

# X4 is a normally distributed continuous covariate
# that varies within studies

create_power_graph(dat = power_dat, beta = "E5", cov = "X4", alpha = .05)



# Power X4 1 --------------------------------------------------------------

create_power_graph(dat = power_dat, beta = "E1", cov = "X4", alpha = .05)


# Power X5 5 --------------------------------------------------------------

# X5 is a continuous, highly skewed covariate that varies within studies

create_power_graph(dat = power_dat, beta = "F5", cov = "X5", alpha = .05)


# Power X5 1 --------------------------------------------------------------

create_power_graph(dat = power_dat, beta = "F1", cov = "X5", alpha = .05)





