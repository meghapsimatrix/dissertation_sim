library(tidyverse)

load("sim_tacc/data/to_test.RData")


# load results ------------------------------------------------------------

files <- list.files("sim_results/results", full.names = TRUE)

files_james <- list.files("sim_results/results_james", full.names = TRUE)

all_files <- c(files, files_james)

load_res <- function(file){
  
  load(file)
  results
  
}


results <- map_dfr(all_files, load_res)
results <- results %>%
  filter(!(test %in% c("CWB Adjusted fixed tau", "CWB fixed tau")))


# check iterations --------------------------------------------------------

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


# Naive -------------------------------------------------------------------

mcse_01 <- sqrt((.01 * (1 - .01))/ K_all)
mcse_05 <- sqrt((.05 * (1 - .05))/ K_all)
mcse_10 <- sqrt((.10 * (1 - .10))/ K_all)



data_int <- tibble(alpha = c(".01", 
                             ".05", 
                             ".10"),
                   int = c(.01, .05, .10),
                   mcse = c(mcse_01, mcse_05, mcse_10)) %>%
  mutate(error = int + 1.96 * mcse)

naive_dat <- results %>%
  select(test, beta_type, cov_test, rho, tau, m, contrasts, starts_with("rej_rate")) %>%
  filter(test == "Naive-F", beta_type == "A") %>%
  mutate(m = as.character(m),
         q = paste("q =", contrasts)) %>%
  gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
  mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                           alpha == "rej_rate_05" ~ ".05",
                           alpha == "rej_rate_10" ~ ".10"))  %>%
  group_by(m, rho, tau, q, cov_test, alpha) %>%
  summarize_at(vars(rej_rate), mean) %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))

summary(naive_dat$mcse)


naive_dat %>%
  ggplot(aes(x = m, y = rej_rate, fill = m)) +
    geom_hline(data = data_int, aes(yintercept = int), linetype = "solid") + 
    geom_hline(data = data_int, aes(yintercept = error), linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .1)) + 
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(alpha ~ q, scales = "free_y",  labeller = label_bquote(rows = alpha == .(alpha))) + 
    labs(x = "Number of Studies", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  

ggsave("sim_results/graphs_paper/study_1/naivef.png", device = "png", dpi = 500, height = 7, width = 12)


naive_dat %>%
  ungroup() %>%
  group_by(alpha) %>%
  summarize(min = min(mcse),
            max = max(mcse))

# Rej Rate Mean -----------------------------------------------------------


type1_dat <- results %>%
    filter(beta_type == "A") %>%
    mutate(m = paste("m =", m),
           q = paste("q =", contrasts)) %>%
    gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
    mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                             alpha == "rej_rate_05" ~ ".05",
                             alpha == "rej_rate_10" ~ ".10")) %>%
    group_by(m, tau, rho, q, cov_test, test, alpha) %>% 
    summarize(rej_rate = mean(rej_rate),
              .groups = "drop") %>%
    mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))


mcse_type_1 <- type1_dat %>%
  ungroup() %>%
  filter(test != "Naive-F") %>%
  group_by(test, alpha) %>%
  summarize(max = max(mcse)) %>%
  spread(alpha, max) %>%
  mutate_if(is.numeric, round, 3)

#write_csv(mcse_type_1, "sim_results/mcse_type_1.csv")


power_dat <- results %>%
  filter(beta_type != "A") %>%
  mutate(q = paste("q =", contrasts),
         m = as.character(m),
         rho = as.character(rho),
         tau = as.character(tau)) %>%
  gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
  mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                           alpha == "rej_rate_05" ~ ".05",
                           alpha == "rej_rate_10" ~ ".10")) %>%
  group_by(m, tau, rho, q, beta_type, cov_test, test, alpha) %>% 
  summarize(rej_rate = mean(rej_rate),
            .groups = "drop") %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))
  
power_mcse <- power_dat %>%
  ungroup() %>%
  filter(test != "Naive-F") %>%
  group_by(test, alpha) %>%
  summarize(max = max(mcse)) %>%
  spread(alpha, max) %>%
  mutate_if(is.numeric, round, 3)

power_mcse <- power_dat %>%
  ungroup() %>%
  filter(test != "Naive-F") %>%
  group_by(test, alpha, q) %>%
  summarize(max = max(mcse)) %>%
  spread(alpha, max) %>%
  mutate_if(is.numeric, round, 3)



# Alpha .05 

# Type 1 error ------------------------------------------------------------


create_type1_graph <- function(dat, intercept, error){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) +
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    scale_fill_brewer(palette = "Set1") +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
}



create_type1_graph(dat = type1_dat %>% filter(alpha == ".05"), 
                   intercept = .05, 
                   error = data_int %>% filter(int == .05) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/type1_05.png", device = "png", dpi = 500, height = 7, width = 12)

# Type 1 error ------------------------------------------------------------
# 01

create_type1_graph(dat = type1_dat %>% filter(alpha == ".01"),
                   intercept = .01,
                   error = data_int %>% filter(int == .01) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/type1_01.png", device = "png", dpi = 500, height = 7, width = 12)

# Type 1 error ------------------------------------------------------------
# 10

create_type1_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10,
                   error = data_int %>% filter(int == .10) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/type1_10.png", device = "png", dpi = 500, height = 7, width = 12)



# Power -------------------------------------------------------------------

create_power_graph <- function(alpha_level, dat = power_dat){
  
  dat %>%
    filter(alpha == alpha_level) %>%
    filter(test != "Naive-F") %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    ggplot(aes(x = m, y = rej_rate, fill = test)) + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power", fill = "") + 
    theme_bw() +
    theme(legend.position = "bottom",
          plot.caption=element_text(hjust = 0, size = 10))
  
}

glimpse(power_dat)

create_power_graph(alpha_level = ".05")
ggsave("sim_results/graphs_paper/study_1/power_05_abs.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_graph(alpha_level = ".01")
ggsave("sim_results/graphs_paper/study_1/power_01_abs.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_graph(alpha_level = ".10")
ggsave("sim_results/graphs_paper/study_1/power_10_abs.png", device = "png", dpi = 500, height = 7, width = 12)




# Power ratio -------------------------------------------------------------

# the 10 study one throws it off so bad

power_ratio <- power_dat %>%
  filter(test %in% c("CWB", "HTZ")) %>%
  select(-c(starts_with("mcse"))) %>%
  spread(test, rej_rate) %>%
  mutate(power_diff = CWB - HTZ,
         power_ratio = HTZ/CWB) %>%
  group_by(m, rho, tau, alpha, q, beta_type, cov_test) %>%
  summarize_at(vars(power_ratio), mean)



create_power_rat_graph <- function(dat, alpha_level){
  
  
  dat %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    ggplot(aes(x = m, y = power_ratio, fill = m)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/ CWB") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
  
}


create_power_rat_graph(power_ratio, alpha_level = ".05")
ggsave("sim_results/graphs_paper/study_1/power_05.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph(power_ratio, alpha_level = ".01")
ggsave("sim_results/graphs_paper/study_1/power_01.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph(power_ratio, alpha_level = ".10")
ggsave("sim_results/graphs_paper/study_1/power_10.png", device = "png", dpi = 500, height = 7, width = 12)




# Sensitivity Analyses ---------------------------------------------------

# Alpha .05 

# Type 1 error ------------------------------------------------------------

create_type1_tau_graph <- function(dat, intercept, error, br){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    mutate(tau = as.factor(tau)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = tau)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("plum2", "plum4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(tau)) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".05"), intercept = .05, 
                       error = data_int %>% filter(int == .05) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/tau_05.png", device = "png", dpi = 500, height = 7, width = 12)

create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".01"), intercept = .01,
                       error = data_int %>% filter(int == .01) %>% pull(error))
create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10,
                       error = data_int %>% filter(int == .10) %>% pull(error))





create_type1_rho_graph <- function(dat, intercept, error){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    mutate(rho = as.factor(rho)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = rho)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5, position = "dodge") + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("firebrick2", "firebrick4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(rho)) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".05"), intercept = .05,
                       error = data_int %>% filter(int == .05) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/rho_05.png", device = "png", dpi = 500, height = 7, width = 12)

create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".01"), intercept = .01, 
                       error = data_int %>% filter(int == .01) %>% pull(error))
create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10,
                       error = data_int %>% filter(int == .10) %>% pull(error))


# Power sensitivity -------------------------------------------------------



create_power_tau_graph <- function(dat, alpha_level){
  
  
  dat %>%
    mutate(tau = as.factor(tau)) %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    ggplot(aes(x = m, y = power_ratio, fill = tau)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_manual(values = c("plum2", "plum4")) +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/CWB", fill = expression(tau)) + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}

create_power_tau_graph(dat = power_ratio, alpha_level = ".05")
ggsave("sim_results/graphs_paper/study_1/tau_power_05.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_tau_graph(dat = power_ratio, alpha_level = ".01")
create_power_tau_graph(dat = power_ratio, alpha_level = ".10")


create_power_rho_graph <- function(dat, alpha_level){
  
  
  dat %>%
    mutate(rho = as.factor(rho)) %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    ggplot(aes(x = m, y = power_ratio, fill = rho)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_manual(values = c("firebrick2", "firebrick4")) +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/CWB", fill = expression(rho)) + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}

create_power_rho_graph(dat = power_ratio, alpha_level = ".05")
ggsave("sim_results/graphs_paper/study_1/rho_power_05.png", device = "png", dpi = 500, height = 7, width = 12)



# Power single coef -------------------------------------------------------


create_power_covs <- function(alpha_level, dat = power_dat, q_level){
  
  dat <- dat %>%
    mutate(cov_test = str_replace_all(cov_test, "\\+", "\\,"))  %>%
    filter(alpha == alpha_level) %>%
    filter(test != "Naive-F") %>%
    filter(q == q_level) 
  
  if(q_level == "q = 1"){
    
    dat <- dat %>%
      mutate(cov_name = case_when(cov_test == "X1" ~ "Study-level, binary, large imbalance",
                                  cov_test == "X2" ~ "Effect size-level, binary, large imbalance",
                                  cov_test == "X3" ~ "Study-level, continuous, normal",
                                  cov_test == "X4" ~ "Effect size-level, continuous, normal",
                                  cov_test == "X5" ~ "Effect size-level, continuous, skewed"))
    
  }
  
  
  p <- dat %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    #filter(beta == beta_level) %>%
    ggplot(aes(x = m, y = rej_rate, fill = test)) + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_fill_brewer(palette = "Set1") 
  
  if(q_level == "q = 1"){
    p <- p + facet_grid(beta ~ cov_test + cov_name, scales = "free_y", labeller = label_bquote(rows = beta == .(beta)))
  } else {
    p <- p + facet_grid(beta ~ cov_test, scales = "free_y", labeller = label_bquote(rows = beta == .(beta)))

  }
  
  p +
    labs(x = "Number of Studies", y = "Power", fill = "") +
    #ggtitle(bquote(beta == .(beta_level))) + 
    theme_bw() +
    theme(legend.position = "bottom",
          plot.caption=element_text(hjust = 0, size = 10))
  
}



create_power_covs(alpha_level = ".05", q_level = "q = 1")
ggsave("sim_results/graphs_paper/study_1/power_05_q1.png", device = "png", dpi = 500, height = 7, width = 13)



# create_power_covs(alpha_level = ".01", q_level = "q = 1")
# ggsave("sim_results/graphs_paper/study_1/power_01_q1.png", device = "png", dpi = 500, height = 7, width = 12)
# 
# 
# create_power_covs(alpha_level = ".10", q_level = "q = 1")
# ggsave("sim_results/graphs_paper/study_1/power_10_q1.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_covs(alpha_level = ".05", q_level = "q = 2")
ggsave("sim_results/graphs_paper/study_1/power_05_q2.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_covs(alpha_level = ".05", q_level = "q = 3")
ggsave("sim_results/graphs_paper/study_1/power_05_q3.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_covs(alpha_level = ".05", q_level = "q = 4")
ggsave("sim_results/graphs_paper/study_1/power_05_q4.png", device = "png", dpi = 500, height = 7, width = 12)

# create_power_covs(alpha_level = ".05", q_level = "q = 5")
# ggsave("sim_results/graphs_paper/study_1/power_05_q5.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_rat_graph_covs <- function(alpha_level, q_level, dat = power_ratio){
  
  
  dat <- dat %>%
    mutate(cov_test = str_replace_all(cov_test, "\\+", "\\,"))  %>%
    filter(alpha == alpha_level) %>%
    filter(q == q_level) %>%
    mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) 
  # filter(beta == beta_level) %>%
  
  if(q_level == "q = 1"){
    
    dat <- dat %>%
      mutate(cov_name = case_when(cov_test == "X1" ~ "Study-level, binary, large imbalance",
                                  cov_test == "X2" ~ "Effect size-level, binary, large imbalance",
                                  cov_test == "X3" ~ "Study-level, continuous, normal",
                                  cov_test == "X4" ~ "Effect size-level, continuous, normal",
                                  cov_test == "X5" ~ "Effect size-level, continuous, skewed"))
    
  }
  
  p <- dat %>%
    ggplot(aes(x = m, y = power_ratio, fill = m)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    ylim(c(0, NA)) +
    scale_fill_brewer(palette = "Dark2") 
  
  if(q_level == "q = 1"){
    p <- p + facet_grid(beta ~ cov_test + cov_name, scales = "free_y", labeller = label_bquote(rows = beta == .(beta)))
  } else {
    p <- p + facet_grid(beta ~ cov_test, scales = "free_y", labeller = label_bquote(rows = beta == .(beta)))
    
  }
  
  p +
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/ CWB") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
  
}

create_power_rat_graph_covs(alpha_level = ".05", q_level = "q = 1")
ggsave("sim_results/graphs_paper/study_1/power_rat_05_q1.png", device = "png", dpi = 500, height = 7, width = 13)

# create_power_rat_graph_covs(alpha_level = ".01", q_level = "q = 1")
# ggsave("sim_results/graphs_paper/study_1/power_rat_01_q1.png", device = "png", dpi = 500, height = 7, width = 12)
# 
# 
# create_power_rat_graph_covs(alpha_level = ".10", q_level = "q = 1")
# ggsave("sim_results/graphs_paper/study_1/power_rat_10_q1.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_rat_graph_covs(alpha_level = ".05", q_level = "q = 2")
ggsave("sim_results/graphs_paper/study_1/power_rat_05_q2.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph_covs(alpha_level = ".05", q_level = "q = 3")
ggsave("sim_results/graphs_paper/study_1/power_rat_05_q3.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph_covs(alpha_level = ".05", q_level = "q = 4")
ggsave("sim_results/graphs_paper/study_1/power_rat_05_q4.png", device = "png", dpi = 500, height = 7, width = 12)


# create_power_rat_graph_covs(alpha_level = ".05", q_level = "q = 5")
# ggsave("sim_results/graphs_paper/study_1/power_rat_05_q5.png", device = "png", dpi = 500, height = 7, width = 12)



# type 1 covs -------------------------------------------------------------

create_type1_covs <- function(q_level,
                              intercept, error, dat = power_dat){
  
  dat <- dat %>%
    mutate(cov_test = str_replace_all(cov_test, "\\+", "\\,"))  %>%
    #filter(alpha == alpha_level) %>%
    filter(test != "Naive-F") %>%
    filter(q == q_level) 
  
  if(q_level == "q = 1"){
    
    dat <- dat %>%
      mutate(cov_name = case_when(cov_test == "X1" ~ "Study-level, binary, large imbalance",
                                  cov_test == "X2" ~ "Effect size-level, binary, large imbalance",
                                  cov_test == "X3" ~ "Study-level, continuous, normal",
                                  cov_test == "X4" ~ "Effect size-level, continuous, normal",
                                  cov_test == "X5" ~ "Effect size-level, continuous, skewed"))
    
  }
  
  
  p <- dat %>%
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5) +
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    scale_fill_brewer(palette = "Set1") 
  

  if(q_level == "q = 1"){
    p <- p + facet_grid(m ~ cov_test + cov_name, scales = "free_y")
  } else {
    p <- p + facet_grid(m ~ cov_test, scales = "free_y")
    
  }
  
  p +
    labs(x = "Method", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  
}


create_type1_covs(q_level = "q = 1",
                   intercept = .05, 
                   error = data_int %>% filter(int == .05) %>% pull(error),
                   dat = type1_dat %>% filter(alpha == ".05"))


create_power_covs(alpha_level = ".05", q_level = "q = 1")
ggsave("sim_results/graphs_paper/study_1/type1_05_q1.png", device = "png", dpi = 500, height = 7, width = 13)

