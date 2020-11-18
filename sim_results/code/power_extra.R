# Power X1 5 --------------------------------------------------------------

create_power_graph <- function(dat, beta, cov){
  
  dat %>%
    filter(beta_type == beta) %>%
    filter(test != "Naive-F") %>%
    mutate(m = paste("m =", m),
           q = paste("q =", contrasts))  %>%
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_hline(yintercept = .8, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .1)) +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = paste("Power", beta, cov)) + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
}


# The first binary covariate, X1, is a study level covariate with large imbalance, 
# equaling 1 in 15% of the studies.

create_power_graph(dat = small_res_05, beta = "B5", cov = "X1")
create_power_graph(dat = small_res_01, beta = "B5", cov = "X1")


# Power X1 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "B1", cov = "X1")
create_power_graph(dat = small_res_01, beta = "B1", cov = "X1")


# Power X2 5 --------------------------------------------------------------
# The second binary covariate, X2, is a covariate that varies within studies, 
# equaling 1 in 10% of the effect size estimates overall and in 0 to 20% 
# of the effect size estimates within a study.

create_power_graph(dat = small_res_05, beta = "C5", cov = "X2")
create_power_graph(dat = small_res_01, beta = "C5", cov = "X2")


# Power X2 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "C1", cov = "X2")
create_power_graph(dat = small_res_01, beta = "C1", cov = "X2")


# Power X3 5 --------------------------------------------------------------

# X3 is a normally distributed study level covariate

create_power_graph(dat = small_res_05, beta = "D5", cov = "X3")
create_power_graph(dat = small_res_01, beta = "D5", cov = "X3")


# Power X3 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "D1", cov = "X3")
create_power_graph(dat = small_res_01, beta = "D1", cov = "X3")


# Power X4 5 --------------------------------------------------------------

# X4 is a normally distributed continuous covariate
# that varies within studies

create_power_graph(dat = small_res_05, beta = "E5", cov = "X4")
ggsave("sim_results/graphs/power_X4_beta5.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_graph(dat = small_res_01, beta = "E5", cov = "X4")


# Power X4 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "E1", cov = "X4")
create_power_graph(dat = small_res_01, beta = "E1", cov = "X4")


# Power X5 5 --------------------------------------------------------------

# X5 is a continuous, highly skewed covariate that varies within studies

create_power_graph(dat = small_res_05, beta = "F5", cov = "X5")
create_power_graph(dat = small_res_01, beta = "F5", cov = "X5")


# Power X5 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "F1", cov = "X5")
create_power_graph(dat = small_res_01, beta = "F1", cov = "X5")




t.test(rej_rate ~ tau, data = naive_dat)
t.test(rej_rate ~ rho, data = naive_dat)


t.test(rej_rate ~ tau, data = small_res_05)
t.test(rej_rate ~ rho, data = small_res_05)



