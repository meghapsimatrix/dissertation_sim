library(tidyverse)

load("sim_tacc/data/to_test.RData")


# load results ------------------------------------------------------------

files <- list.files("sim_results/results", full.names = TRUE)

load_res <- function(file){
  
  load(file)
  results
  
}


results <- map_dfr(files, load_res)


results %>%
  filter(K < 50)

table(results$batch)

K_all <- results %>%
  select(batch, iterations) %>%
  distinct(.) %>%
  summarize(iterations = sum(iterations)) %>%
  pull(iterations) 

K_all



results <- left_join(results, 
                     to_test %>% select(cov_test, contrasts),
                     by = "cov_test")

glimpse(results)
  



# Alpha .05 

# Type 1 error ------------------------------------------------------------

small_res_05 <- results %>%
  group_by(m, tau, rho, cov_test, beta_type, contrasts, test) %>%
  summarize(rej_rate = mean(rej_rate_05),
            .groups = "drop") %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))

summary(small_res_05$mcse)

create_type1_graph <- function(dat, intercept){
  
  dat %>%
    filter(beta_type == "A") %>%
    mutate(m = paste("m =", m),
           q = paste("q =", contrasts))  %>%
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_hline(yintercept = intercept, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, .6, .1)) + 
    scale_fill_brewer(palette = "Set1") +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
}



create_type1_graph(dat = small_res_05, intercept = .05)



ggsave("sim_results/graphs/type1.png", device = "png", dpi = 500, height = 7, width = 12)






# Power X1 5 --------------------------------------------------------------

create_power_graph <- function(dat, beta, cov){
  
  dat %>%
    filter(beta_type == beta, str_detect(cov_test, cov)) %>%
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



# Power X1 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "B1", cov = "X1")


# Power X2 5 --------------------------------------------------------------
# The second binary covariate, X2, is a covariate that varies within studies, 
# equaling 1 in 10% of the effect size estimates overall and in 0 to 20% 
# of the effect size estimates within a study.

create_power_graph(dat = small_res_05, beta = "C5", cov = "X2")

# Power X2 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "C1", cov = "X2")

# Power X3 5 --------------------------------------------------------------

# X3 is a normally distributed study level covariate

create_power_graph(dat = small_res_05, beta = "D5", cov = "X3")

# Power X3 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "D1", cov = "X3")

# Power X4 5 --------------------------------------------------------------

# X4 is a normally distributed continuous covariate
# that varies within studies

create_power_graph(dat = small_res_05, beta = "E5", cov = "X4")

ggsave("sim_results/graphs/power_X4_beta5.png", device = "png", dpi = 500, height = 7, width = 12)


# Power X4 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "E1", cov = "X4")

# Power X5 5 --------------------------------------------------------------

# X5 is a continuous, highly skewed covariate that varies within studies

create_power_graph(dat = small_res_05, beta = "F5", cov = "X5")

# Power X5 1 --------------------------------------------------------------

create_power_graph(dat = small_res_05, beta = "F1", cov = "X5")


# Type 1 error .01 --------------------------------------------------------

small_res_01 <- results %>%
  group_by(m, tau, rho, cov_test, beta_type, contrasts, test) %>%
  summarize(rej_rate = mean(rej_rate_01),
            .groups = "drop") %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all)) 

create_type1_graph(dat = small_res_01, intercept = .01)


# Power Difference --------------------------------------------------------

small_res_05_diff <- small_res_05 %>%
  filter(test %in% c("CWB", "HTZ"), beta_type != "A") %>%
  select(-mcse) %>%
  spread(test, rej_rate) %>%
  mutate(power_diff = CWB - HTZ,
         power_ratio = CWB / HTZ) %>%
  group_by(m, rho, tau, contrasts) %>%
  summarize_at(vars(power_diff), mean)

small_res_05_diff %>%
  filter(power_diff < 0) %>%
  View()

small_res_05_diff %>%
  filter(power_ratio < 1) %>%
  View()



