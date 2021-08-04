library(tidyverse)

load("sim_tacc/data/to_test.RData")


# load results ------------------------------------------------------------

files <- list.files("sim_results/results", full.names = TRUE)

files_james <- list.files("sim_results/results_james", full.names = TRUE)

all_files <- c(files, files_james)

load_res <- function(file) {
  
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

data_int <- 
  tibble(alpha = c(".01", ".05", ".10"),
         int = c(.01, .05, .10)) %>%
  mutate(
    mcse = sqrt(int * (1 - int) / K_all),
    error = int + 1.96 * mcse
  )

naive_dat <- 
  results %>%
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
  ggplot(aes(x = m, y = rej_rate, fill = m, color = m)) +
    geom_hline(data = data_int, aes(yintercept = int), linetype = "solid") + 
    geom_hline(data = data_int, aes(yintercept = error), linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .1)) + 
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(alpha ~ q, scales = "free_y",  labeller = label_bquote(rows = alpha == .(alpha))) + 
    labs(x = "Number of Studies", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  

ggsave("sim_results/graphs_paper/study_1/naivef.png", device = "png", dpi = 500, height = 5, width = 7)

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
    ggplot(aes(x = test, y = rej_rate, fill = test, color = test)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) +
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
}



create_type1_graph(dat = type1_dat %>% filter(alpha == ".05"), 
                   intercept = .05, 
                   error = data_int %>% filter(int == .05) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/type1_05.png", device = "png", dpi = 500, height = 5, width = 7)

# Type 1 error ------------------------------------------------------------
# 01

create_type1_graph(dat = type1_dat %>% filter(alpha == ".01"),
                   intercept = .01,
                   error = data_int %>% filter(int == .01) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/type1_01.png", device = "png", dpi = 500, height = 5, width = 7)

# Type 1 error ------------------------------------------------------------
# 10

create_type1_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10,
                   error = data_int %>% filter(int == .10) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/type1_10.png", device = "png", dpi = 500, height = 5, width = 7)



# Power ratio -------------------------------------------------------------

# the 10 study one throws it off so bad

power_ratio <- 
  power_dat %>%
  select(-c(starts_with("mcse"))) %>%
  spread(test, rej_rate) %>%
  mutate(power_diff = CWB - HTZ,
         power_ratio = HTZ / CWB)



# Power scatterplots ----------------------------------------------------------

power_scatter <- function(data, x, y) {
  
  ggplot(data, aes_string(x, y, color = "m", shape = "m")) + 
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0) + 
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0)) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0)) + 
    facet_wrap(~ q, scales = "free") + 
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c("diamond","circle","triangle","square")) + 
    labs(
      x = paste("Power of", x), 
      y = paste("Power of", y),
      color = "Number of studies (m)", shape = "Number of studies (m)"
    ) + 
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.caption=element_text(hjust = 0, size = 10)
    )  
}

# Power comparison at alpha = .05 (for main text)
power_ratio %>%
  filter(alpha == ".05") %>%
  power_scatter(x = "HTZ", y = "CWB")

ggsave("sim_results/graphs_paper/study_1/power_05_scatter.png", device = "png", dpi = 500, height = 7, width = 12)


# Power comparison at alpha = .01 (for supplementary)  
power_ratio %>%
  filter(alpha == ".01") %>%
  power_scatter(x = "HTZ", y = "CWB")

ggsave("sim_results/graphs_paper/study_1/power_01_scatter.png", device = "png", dpi = 500, height = 7, width = 12)


# Power comparison at alpha = .10 (for supplementary)  
power_ratio %>%
  filter(alpha == ".10") %>%
  power_scatter(x = "HTZ", y = "CWB")

ggsave("sim_results/graphs_paper/study_1/power_10_scatter.png", device = "png", dpi = 500, height = 7, width = 12)


# By covariate combination (all for supplementary)

power_ratio %>%
  filter(alpha == ".05") %>%
  mutate(q_cov = paste0(q, " (", cov_test,")")) %>%
  power_scatter(x = "HTZ", y = "CWB") + 
  facet_wrap(~ q_cov, ncol = 5) 

ggsave("sim_results/graphs_paper/study_1/power_05_scatter_covs.png", device = "png", dpi = 500, height = 7, width = 12)


power_ratio %>%
  filter(alpha == ".01") %>%
  mutate(q_cov = paste0(q, " (", cov_test,")")) %>%
  power_scatter(x = "HTZ", y = "CWB") + 
  facet_wrap(~ q_cov, ncol = 5) 

ggsave("sim_results/graphs_paper/study_1/power_01_scatter_covs.png", device = "png", dpi = 500, height = 7, width = 12)


power_ratio %>%
  filter(alpha == ".10") %>%
  mutate(q_cov = paste0(q, " (", cov_test,")")) %>%
  power_scatter(x = "HTZ", y = "CWB") + 
  facet_wrap(~ q_cov, ncol = 5) 

ggsave("sim_results/graphs_paper/study_1/power_10_scatter_covs.png", device = "png", dpi = 500, height = 7, width = 12)


# CWB versus CWB-adjusted (for supplementary)

power_ratio %>%
  filter(alpha == ".05") %>%
  power_scatter(x = "CWB", y = "`CWB Adjusted`") + 
  labs(y = "Power of CWB Adjusted")

ggsave("sim_results/graphs_paper/study_1/power_05_scatter_cwbs.png", device = "png", dpi = 500, height = 7, width = 12)


# Sensitivity Analyses ---------------------------------------------------

# Alpha .05 

# Type 1 error ------------------------------------------------------------

create_type1_tau_graph <- function(dat, intercept, error, br){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    mutate(tau = as.factor(tau)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = tau, color = tau)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("plum2", "plum4")) +
    scale_color_manual(values = c("plum2", "plum4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(tau), color = expression(tau)) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".05"), intercept = .05, 
                       error = data_int %>% filter(int == .05) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/tau_05.png", device = "png", dpi = 500, height = 5, width = 7)

create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".01"), intercept = .01,
                       error = data_int %>% filter(int == .01) %>% pull(error))
create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10,
                       error = data_int %>% filter(int == .10) %>% pull(error))





create_type1_rho_graph <- function(dat, intercept, error){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    mutate(rho = as.factor(rho)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = rho, color = rho)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5, position = "dodge") + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("firebrick2", "firebrick4")) +
    scale_color_manual(values = c("firebrick2", "firebrick4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(rho), color = expression(rho)) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".05"), intercept = .05,
                       error = data_int %>% filter(int == .05) %>% pull(error))

ggsave("sim_results/graphs_paper/study_1/rho_05.png", device = "png", dpi = 500, height = 5, width = 7)

create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".01"), intercept = .01, 
                       error = data_int %>% filter(int == .01) %>% pull(error))
create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".10"), intercept = .10,
                       error = data_int %>% filter(int == .10) %>% pull(error))


