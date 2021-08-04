library(tidyverse)
library(patchwork)


# load and clean data -----------------------------------------------------

load("sim_results/results_study2/sim_test_study_2.RData")

K_all <- results %>%
  select(batch, iterations) %>%
  distinct(.) %>%
  summarize(iterations = sum(iterations)) %>%
  pull(iterations) 

K_all

results <- results %>%
  mutate(cat_num = cat_num - 1)



# mcse --------------------------------------------------------------------


mcse_01 <- sqrt((.01 * (1 - .01))/ K_all)
mcse_05 <- sqrt((.05 * (1 - .05))/ K_all)
mcse_10 <- sqrt((.10 * (1 - .10))/ K_all)



data_int <- tibble(alpha = c(".01", 
                             ".05", 
                             ".10"),
                   int = c(.01, .05, .10),
                   mcse = c(mcse_01, mcse_05, mcse_10)) %>%
  mutate(error = int + 1.96 * mcse)


# naive f type 1 ----------------------------------------------------------

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

create_naive_graph <- function(level){

  naive_dat %>%
    filter(alpha == level) %>%
    mutate(cov_type = if_else(cov_type == "between", "Study-level",
                              "Effect size-level")) %>%
    ggplot(aes(x = m, y = rej_rate, color = m)) +
    geom_hline(data = data_int %>% filter(alpha == level), aes(yintercept = int), linetype = "solid") + 
    geom_hline(data = data_int %>% filter(alpha == level), aes(yintercept = error), linetype = "dashed") + 
    geom_jitter(alpha = .5, width = 0.1, height =  NULL) + 
    ylim(c(0, NA)) +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(cov_type ~ q, scales = "free_y") + 
    labs(x = "Number of Studies", y = "Type 1 Error Rate") + 
    #ggtitle(title) +
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
}  


create_naive_graph(level = ".05")

ggsave("sim_results/graphs_paper/study_2/naivef_05_2.png", device = "png", dpi = 500, height = 5, width = 7)


create_naive_graph(level = ".01")

ggsave("sim_results/graphs_paper/study_2/naivef_01_2.png", device = "png", dpi = 500, height = 5, width = 7)


create_naive_graph(level = ".10")

ggsave("sim_results/graphs_paper/study_2/naivef_10_2.png", device = "png", dpi = 500, height = 5, width = 7)


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

mcse_type_1 <- type1_dat %>%
  ungroup() %>%
  filter(test != "Naive-F") %>%
  group_by(test, alpha) %>%
  summarize(max = max(mcse)) %>%
  spread(alpha, max) %>%
  mutate_if(is.numeric, round, 3)

#write_csv(mcse_type_1, "sim_results/mcse_type_1_study2.csv")

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

power_mcse <- power_dat %>%
  ungroup() %>%
  filter(test != "Naive-F") %>%
  group_by(test, alpha) %>%
  summarize(max = max(mcse)) %>%
  spread(alpha, max) %>%
  mutate_if(is.numeric, round, 3)

# Alpha .05 

# Type 1 error ------------------------------------------------------------


create_type1_graph <- function(dat, intercept, error, br, cov, title){
  
  dat <- dat %>%
    filter(test != "Naive-F") %>%
    filter(cov_type == cov)
  
  dat %>%
    ggplot(aes(x = test, y = rej_rate, color = test)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_jitter(alpha = .5, width = 0.1, height =  NULL) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    scale_color_brewer(palette = "Set1") +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate") + 
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
}


b <- create_type1_graph(dat = type1_dat %>% filter(alpha == ".05"),
                        intercept = .05,
                        error = data_int %>% filter(int == .05) %>% pull(error), 
                        cov = "between",
                        title = "Study-level covariate type")

w <- create_type1_graph(dat = type1_dat %>% filter(alpha == ".05"),
                        intercept = .05,
                        error = data_int %>% filter(int == .05) %>% pull(error), 
                        cov = "within",
                        title = "Effect size-level covariate type")

w / b

ggsave("sim_results/graphs_paper/study_2/type1_05_2.png", device = "png", 
       dpi = 500, height = 8, width = 7)

# Type 1 error ------------------------------------------------------------
# 01

b <- create_type1_graph(dat = type1_dat %>% filter(alpha == ".01"),
                   intercept = .01,
                   error = data_int %>% filter(int == .01) %>% pull(error), 
                   cov = "between",
                   title = "Study-level covariate type")

w <- create_type1_graph(dat = type1_dat %>% filter(alpha == ".01"),
                        intercept = .01,
                        error = data_int %>% filter(int == .01) %>% pull(error), 
                        cov = "within",
                        title = "Effect size-level covariate type")

w / b

ggsave("sim_results/graphs_paper/study_2/type1_01_2.png", device = "png", 
       dpi = 500, height = 8, width = 7)
# Type 1 error ------------------------------------------------------------
# 10

b <- create_type1_graph(dat = type1_dat %>% filter(alpha == ".10"), 
                   intercept = .10,
                   error = data_int %>% filter(int == .10) %>% pull(error),
                   cov = "between",
                   title = "Study-level covariate type")

w <- create_type1_graph(dat = type1_dat %>% filter(alpha == ".10"), 
                        intercept = .10,
                        error = data_int %>% filter(int == .10) %>% pull(error),
                        cov = "within",
                        title = "Effect size-level covariate type")

w / b

ggsave("sim_results/graphs_paper/study_2/type1_10_2.png", device = "png", dpi = 500, height = 8, width = 7)


# Power scatterplots ----------------------------------------------------------

power_ratio <- 
  power_dat %>%
  select(-c(starts_with("mcse"))) %>%
  spread(test, rej_rate) %>%
  mutate(power_diff = CWB - HTZ,
         power_ratio = HTZ / CWB)

power_scatter <- function(data, x, y) {
  
  data %>%
    mutate(cov_type = if_else(cov_type == "between", "Study-level", "Effect size-level")) %>%
  ggplot(aes_string(x, y, color = "m", shape = "m")) + 
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0) + 
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0)) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0)) + 
    facet_grid(cov_type ~ q, scales = "free") + 
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c("diamond","circle","triangle","square")) + 
    labs(
      x = paste("Power of", x), 
      y = paste("Power of", y),
      color = "Number of studies (m)", shape = "Number of studies (m)"
    ) + 
   # ggtitle(title) +
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


ggsave("sim_results/graphs_paper/study_2/power_05_scatter_2.png", device = "png", dpi = 500, height = 5, width = 7)


# Power comparison at alpha = .01 (for supplementary)  
power_ratio %>%
  filter(alpha == ".01") %>%
  power_scatter(x = "HTZ", y = "CWB")


ggsave("sim_results/graphs_paper/study_2/power_01_scatter_2.png", device = "png", dpi = 500, height = 5, width = 7)

# Power comparison at alpha = .10 (for supplementary)  
power_ratio %>%
  filter(alpha == ".10") %>%
  power_scatter(x = "HTZ", y = "CWB")


ggsave("sim_results/graphs_paper/study_2/power_10_scatter_2.png", device = "png", dpi = 500, height = 5, width = 7)


# CWB versus CWB-adjusted (for supplementary)

power_ratio %>%
  filter(alpha == ".05") %>%
  power_scatter(x = "CWB", y = "`CWB Adjusted`") + 
  labs(y = "Power of CWB-adjusted")


ggsave("sim_results/graphs_paper/study_2/power_05_scatter_cwbs_2.png", device = "png", dpi = 500, height = 7, width = 12)


# Sensitivity Analyses ---------------------------------------------------

# Alpha .05 

# Type 1 error ------------------------------------------------------------

create_type1_tau_graph <- function(dat, intercept, error, br, cov, title){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(tau = as.factor(tau)) %>%
    ggplot(aes(x = test, y = rej_rate, color = tau)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_jitter(alpha = .5, width = 0.1, height =  NULL) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_color_manual(values = c("plum2", "plum4")) +
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", color = expression(tau)) + 
    ggtitle(title) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

b <- create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".05"), 
                       intercept = .05, 
                       error = data_int %>% filter(int == .05) %>% pull(error),
                       cov = "between",
                       title = "Study-level covariate type")



w <- create_type1_tau_graph(dat = type1_dat %>% filter(alpha == ".05"), 
                       intercept = .05, 
                       error = data_int %>% filter(int == .05) %>% pull(error),
                       cov = "within",
                       title = "Effect size-level covariate type")

w / b + plot_layout(guides = "collect") & theme(legend.position = 'bottom')


ggsave("sim_results/graphs_paper/study_2/tau_052.png", device = "png", dpi = 500, height = 8, width = 7)



create_type1_rho_graph <- function(dat, intercept, error, br, cov, title){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(rho = as.factor(rho)) %>%
    ggplot(aes(x = test, y = rej_rate, color = rho)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_jitter(alpha = .5, width = 0.1, height =  NULL) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    scale_color_manual(values = c("firebrick2", "firebrick4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", color = expression(rho)) + 
    ggtitle(title) + 
    theme_bw() +
    theme(legend.position = "bottom")
}

b <- create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".05"),
                       intercept = .05, 
                       error = data_int %>% filter(int == .05) %>% pull(error),
                       br = .02, 
                       cov = "between",
                       title = "Study-level covariate type")


w <- create_type1_rho_graph(dat = type1_dat %>% filter(alpha == ".05"),
                       intercept = .05, 
                       error = data_int %>% filter(int == .05) %>% pull(error),
                       br = .02, 
                       cov = "within",
                       title = "Effect size-level covariate type")

w / b + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs_paper/study_2/rho_052.png", device = "png", dpi = 500, height = 8, width = 7)




