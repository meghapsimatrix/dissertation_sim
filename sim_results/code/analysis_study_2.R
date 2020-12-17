library(tidyverse)

load("sim_results/results_study2/sim_test_study_2.RData")

K_all <- results %>%
  select(batch, iterations) %>%
  distinct(.) %>%
  summarize(iterations = sum(iterations)) %>%
  pull(iterations) 

K_all

results <- results %>%
  filter(test != "EDT")


data_int <- tibble(alpha = c(".01", 
                             ".05", 
                             ".10"),
                   int = c(.01, .05, .10))

naive_dat <- results %>%
  select(test, beta_1, rho, tau, m, cov_type, cat_num, starts_with("rej_rate")) %>%
  filter(test == "Naive-F", beta_1 == 0) %>%
  mutate(m = as.character(m),
         q = paste("q =", cat_num)) %>%
  gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
  mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                           alpha == "rej_rate_05" ~ ".05",
                           alpha == "rej_rate_10" ~ ".10"))  %>%
  group_by(m, rho, tau, cov_type, q, alpha) %>%
  summarize_at(vars(rej_rate), mean) %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))

summary(naive_dat$mcse)

naive_check <- naive_dat %>%
  group_by(q, alpha, m, cov_type) %>%
  summarize(min = min(rej_rate),
            max = max(rej_rate),
            median = median(rej_rate))


naive_dat %>%
  ggplot(aes(x = m, y = rej_rate, fill = cov_type)) +
  geom_hline(data = data_int, aes(yintercept = int), linetype = "dashed") + 
  geom_boxplot(alpha = .5) + 
  scale_y_continuous(breaks = seq(0, 1, .1)) + 
  scale_fill_manual(values = c("aquamarine4", "violetred4")) +
  facet_grid(alpha ~ q, scales = "free_y",  labeller = label_bquote(rows = alpha == .(alpha))) + 
  labs(x = "Number of Studies", y = "Type 1 Error Rate", fill = "Covariate Type") + 
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("sim_results/graphs/study_2/naivef_2.png", device = "png", dpi = 500, height = 7, width = 12)


naive_dat %>%
  ungroup() %>%
  group_by(alpha) %>%
  summarize(min = min(mcse),
            max = max(mcse))


# Rej Rate Mean -----------------------------------------------------------


type1_dat <- results %>%
  filter(beta_1 == 0) %>%
  mutate(m = paste("m =", m),
         q = paste("q =", cat_num)) %>%
  gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
  mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                           alpha == "rej_rate_05" ~ ".05",
                           alpha == "rej_rate_10" ~ ".10")) %>%
  group_by(m, tau, rho, q, cov_type, cat_num, test, alpha) %>% 
  summarize(rej_rate = mean(rej_rate),
            .groups = "drop") %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))


type1_dat %>%
  ungroup() %>%
  filter(test != "Naive-F") %>%
  group_by(alpha) %>%
  summarize(min = min(mcse),
            max = max(mcse))

power_dat <- results %>%
  filter(beta_1 != 0) %>%
  mutate(q = paste("q =", cat_num),
         m = as.character(m),
         rho = as.character(rho),
         tau = as.character(tau)) %>%
  gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
  mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                           alpha == "rej_rate_05" ~ ".05",
                           alpha == "rej_rate_10" ~ ".10")) %>%
  group_by(m, tau, rho, q, beta_1, cov_type, cat_num, test, alpha) %>% 
  summarize(rej_rate = mean(rej_rate),
            .groups = "drop") %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))

power_dat %>%
  ungroup() %>%
  group_by(test, beta_1, alpha) %>%
  summarize(min = min(mcse),
            max = max(mcse))


# Alpha .05 

# Type 1 error ------------------------------------------------------------


create_type1_graph <- function(dat, intercept, br){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    ggplot(aes(x = test, y = rej_rate, fill = cov_type)) + 
    geom_hline(yintercept = intercept, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("aquamarine4", "violetred4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = "Covariate Type") + 
    theme_bw() +
    theme(legend.position = "bottom")
}



create_type1_graph(dat = type1_dat %>% filter(alpha == ".05"), intercept = .05, br = .02)

ggsave("sim_results/graphs/study_2/type1_05_2.png", device = "png", dpi = 500, height = 7, width = 12)

# Type 1 error ------------------------------------------------------------
# 01

create_type1_graph(dat = type1_dat %>% filter(alpha == ".01"), intercept = .01, br = .01)

ggsave("sim_results/graphs/study_2/type1_01_2.png", device = "png", dpi = 500, height = 7, width = 12)

# Type 1 error ------------------------------------------------------------
# 10

create_type1_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10, br = .02)

ggsave("sim_results/graphs/study_2/type1_10_2.png", device = "png", dpi = 500, height = 7, width = 12)

# Power ratio -------------------------------------------------------------

# the 10 study one throws it off so bad

power_ratio <- power_dat %>%
  filter(test %in% c("CWB", "HTZ")) %>%
  select(-c(starts_with("mcse"))) %>%
  spread(test, rej_rate) %>%
  mutate(power_diff = CWB - HTZ,
         power_ratio = HTZ/CWB) %>%
  group_by(m, tau, rho, q, beta_1, cov_type, cat_num, alpha) %>%
  summarize_at(vars(power_ratio), mean)


create_power_rat_graph <- function(dat, alpha_level){
  
  
  dat %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = beta_1) %>%
    ggplot(aes(x = m, y = power_ratio, fill = cov_type)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = c("aquamarine4", "violetred4")) +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/ CWB", fill = "Covariate Type") + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}


create_power_rat_graph(power_ratio, alpha_level = ".05")
ggsave("sim_results/graphs/study_2/power_05_2.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph(power_ratio, alpha_level = ".01")
ggsave("sim_results/graphs/study_2/power_01_2.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rat_graph(power_ratio, alpha_level = ".10")
ggsave("sim_results/graphs/study_2/power_10_2.png", device = "png", dpi = 500, height = 7, width = 12)


# Sensitivity Analyses ---------------------------------------------------

# Alpha .05 

# Type 1 error ------------------------------------------------------------

create_type1_tau_graph <- function(dat, intercept, br, cov){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(tau = as.factor(tau)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = tau)) + 
    geom_hline(yintercept = intercept, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("plum3", "plum4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(tau)) + 
    ggtitle(paste0(str_to_title(cov), "-study covariate type")) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".05"), 
                       intercept = .05, 
                       br = .02, 
                       cov = "between")

ggsave("sim_results/graphs/study_2/tau_052_between.png", device = "png", dpi = 500, height = 7, width = 12)


create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".05"), 
                       intercept = .05, 
                       br = .02, 
                       cov = "within")

ggsave("sim_results/graphs/study_2/tau_052_within.png", device = "png", dpi = 500, height = 7, width = 12)

create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".01"), intercept = .01, br = .01)
create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10, br = .02)


create_type1_rho_graph <- function(dat, intercept, br, cov){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(rho = as.factor(rho)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = rho)) + 
    geom_hline(yintercept = intercept, linetype = "dashed") + 
    geom_boxplot(alpha = .5, position = "dodge") + 
    scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("firebrick3", "firebrick4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(rho)) + 
    ggtitle(paste0(str_to_title(cov), "-study covariate type")) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".05"),
                       intercept = .05, 
                       br = .02, 
                       cov = "between")

ggsave("sim_results/graphs/study_2/rho_052_between.png", device = "png", dpi = 500, height = 7, width = 12)

create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".05"),
                       intercept = .05, 
                       br = .02, 
                       cov = "within")

ggsave("sim_results/graphs/study_2/rho_052_within.png", device = "png", dpi = 500, height = 7, width = 12)


create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".01"), intercept = .01, br = .01)
create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10, br = .02)


# Power sensitivity -------------------------------------------------------



create_power_tau_graph <- function(dat, alpha_level, cov){
  
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(tau = as.factor(tau)) %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = beta_1) %>%
    ggplot(aes(x = m, y = power_ratio, fill = tau)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = c("plum3", "plum4")) +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/CWB", fill = expression(tau)) + 
    ggtitle(paste0(str_to_title(cov), "-study covariate type")) + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}

create_power_tau_graph(dat = power_ratio, alpha_level = ".05", cov = "between")
ggsave("sim_results/graphs/study_2/tau_power_052_between.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_tau_graph(dat = power_ratio, alpha_level = ".05", cov = "within")
ggsave("sim_results/graphs/study_2/tau_power_052_within.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_tau_graph(dat = power_ratio, alpha_level = ".01")
create_power_tau_graph(dat = power_ratio, alpha_level = ".10")


create_power_rho_graph <- function(dat, alpha_level, cov){
  
  
  dat %>%
    mutate(cov_type == cov) %>%
    mutate(rho = as.factor(rho)) %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = beta_1) %>%
    ggplot(aes(x = m, y = power_ratio, fill = rho)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "dashed") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_manual(values = c("firebrick3", "firebrick4")) +
    facet_grid(beta ~ q, scales = "free_y",  labeller = label_bquote(rows = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/CWB", fill = expression(rho)) + 
    ggtitle(paste0(str_to_title(cov), "-study covariate type")) + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}

create_power_rho_graph(dat = power_ratio, alpha_level = ".05", cov = "between")
ggsave("sim_results/graphs/study_2/rho_power_052_between.png", device = "png", dpi = 500, height = 7, width = 12)

create_power_rho_graph(dat = power_ratio, alpha_level = ".05", cov = "within")
ggsave("sim_results/graphs/study_2/rho_power_052_within.png", device = "png", dpi = 500, height = 7, width = 12)


create_power_rho_graph(dat = power, alpha_level = ".01")
create_power_rho_graph(dat = power, alpha_level = ".10")

