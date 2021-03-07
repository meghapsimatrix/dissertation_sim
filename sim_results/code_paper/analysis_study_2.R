library(tidyverse)
library(patchwork)


load("sim_results/results_study2/sim_test_study_2.RData")

K_all <- results %>%
  select(batch, iterations) %>%
  distinct(.) %>%
  summarize(iterations = sum(iterations)) %>%
  pull(iterations) 

K_all

results <- results %>%
  mutate(cat_num = cat_num - 1)



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

create_naive_graph <- function(cov, title){

  naive_dat %>%
    filter(cov_type == cov) %>%
    ggplot(aes(x = m, y = rej_rate, fill = m)) +
    geom_hline(data = data_int, aes(yintercept = int), linetype = "solid") + 
    geom_hline(data = data_int, aes(yintercept = error), linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, 1, .1)) + 
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(alpha ~ q, scales = "free_y",  labeller = label_bquote(rows = alpha == .(alpha))) + 
    labs(x = "Number of Studies", y = "Type 1 Error Rate") + 
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
}  

b <- create_naive_graph(cov = "between", 
                        title = "Study-level covariate type")

w <- create_naive_graph(cov = "within", 
                        title = "Effect size-level covariate type")


b + w


ggsave("sim_results/graphs_paper/study_2/naivef_2.png", device = "png", dpi = 500, height = 7, width = 12)


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
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    scale_fill_brewer(palette = "Set1") +
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

b + w

ggsave("sim_results/graphs_paper/study_2/type1_05_2.png", device = "png", dpi = 500, height = 7, width = 14)

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

b + w

ggsave("sim_results/graphs_paper/study_2/type1_01_2.png", device = "png", dpi = 500, height = 7, width = 14)

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

b + w

ggsave("sim_results/graphs_paper/study_2/type1_10_2.png", device = "png", dpi = 500, height = 7, width = 14)


# Power -------------------------------------------------------------------

create_power_graph <- function(alpha_level, dat = power_dat, cov, title){
  
  dat %>%
    filter(cov_type == cov) %>%
    filter(alpha == alpha_level) %>%
    filter(test != "Naive-F") %>%
    mutate(beta = beta_1) %>%
    ggplot(aes(x = m, y = rej_rate, fill = test)) + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(q ~ beta, scales = "free_y",  labeller = label_bquote(cols = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power", fill = "") + 
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.caption=element_text(hjust = 0, size = 10))
  
}



b <- create_power_graph(alpha_level = ".05", 
                   cov = "between",
                   title = "Study-level covariate type")

w <- create_power_graph(alpha_level = ".05", 
                        cov = "within",
                        title = "Effect size-level covariate type")

b + w + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs_paper/study_2/power_05_abs_2.png", device = "png", dpi = 500, height = 7, width = 12)

b <- create_power_graph(alpha_level = ".01", 
                        cov = "between",
                        title = "Study-level covariate type")

w <- create_power_graph(alpha_level = ".01", 
                        cov = "within",
                        title = "Effect size-level covariate type")

b + w + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs_paper/study_2/power_01_abs_2.png", device = "png", dpi = 500, height = 7, width = 12)


b <- create_power_graph(alpha_level = ".10", 
                        cov = "between",
                        title = "Study-level covariate type")

w <- create_power_graph(alpha_level = ".10", 
                        cov = "within",
                        title = "Effect size-level covariate type")

b + w + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs_paper/study_2/power_10_abs_2.png", device = "png", dpi = 500, height = 7, width = 12)




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


create_power_rat_graph <- function(dat, alpha_level, cov, title){
  
  
  dat %>%
    filter(alpha == alpha_level) %>%
    filter(cov_type == cov) %>%
    mutate(beta = beta_1) %>%
    ggplot(aes(x = m, y = power_ratio, fill = m)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    ylim(c(0, NA)) +
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(q ~ beta, scales = "free_y",  
               labeller = label_bquote(cols = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/ CWB") + 
    ggtitle(title) +
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))

}

b <- create_power_rat_graph(power_ratio, alpha_level = ".05", 
                            cov = "between",
                            title = "Study-level covariate type")

w <- create_power_rat_graph(power_ratio, alpha_level = ".05", 
                       cov = "within",
                       title = "Effect size-level covariate type")

b + w

ggsave("sim_results/graphs_paper/study_2/power_05_2.png", device = "png", dpi = 500, height = 7, width = 12)

b <- create_power_rat_graph(power_ratio, alpha_level = ".01",
                       cov = "between",
                       title = "Study-level covariate type")

w <- create_power_rat_graph(power_ratio, alpha_level = ".01",
                            cov = "within",
                            title = "Effect size-level covariate type")

b + w

ggsave("sim_results/graphs_paper/study_2/power_01_2.png", device = "png", dpi = 500, height = 7, width = 12)

b <- create_power_rat_graph(power_ratio, alpha_level = ".10",
                       cov = "between",
                       title = "Study-level covariate type")

w <- create_power_rat_graph(power_ratio, alpha_level = ".10",
                            cov = "within",
                            title = "Effect size-level covariate type")

b + w

ggsave("sim_results/graphs_paper/study_2/power_10_2.png", device = "png", dpi = 500, height = 7, width = 12)


# Sensitivity Analyses ---------------------------------------------------

# Alpha .05 

# Type 1 error ------------------------------------------------------------

create_type1_tau_graph <- function(dat, intercept, error, br, cov, title){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(tau = as.factor(tau)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = tau)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_manual(values = c("plum2", "plum4")) +
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(tau)) + 
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

b + w  + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs_paper/study_2/tau_052.png", device = "png", dpi = 500, height = 7, width = 14)



create_type1_rho_graph <- function(dat, intercept, error, br, cov, title){
  
  dat <- dat %>%
    filter(test != "Naive-F")
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(rho = as.factor(rho)) %>%
    ggplot(aes(x = test, y = rej_rate, fill = rho)) + 
    geom_hline(yintercept = intercept, linetype = "solid") + 
    geom_hline(yintercept = error, linetype = "dashed") + 
    geom_boxplot(alpha = .5, position = "dodge") + 
    #scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) + 
    scale_fill_manual(values = c("firebrick2", "firebrick4")) +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate", fill = expression(rho)) + 
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

b + w + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs_paper/study_2/rho_052.png", device = "png", dpi = 500, height = 7, width = 14)




# Power sensitivity -------------------------------------------------------



create_power_tau_graph <- function(dat, alpha_level, cov, title){
  
  
  dat %>%
    filter(cov_type == cov) %>%
    mutate(tau = as.factor(tau)) %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = beta_1) %>%
    ggplot(aes(x = m, y = power_ratio, fill = tau)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    scale_fill_manual(values = c("plum2", "plum4")) +
    ylim(c(0, NA)) +
    facet_grid(q ~ beta, scales = "free_y",  labeller = label_bquote(cols = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/CWB", fill = expression(tau)) + 
    ggtitle(title) + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}

b <- create_power_tau_graph(dat = power_ratio, 
                            alpha_level = ".05", 
                            cov = "between",
                            title = "Study-level covariate type")

w <- create_power_tau_graph(dat = power_ratio, 
                            alpha_level = ".05", 
                            cov = "within",
                            title = "Effect size-level covariate type")

b + w + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs/study_2/tau_power_052.png", device = "png", dpi = 500, height = 7, width = 12)




create_power_rho_graph <- function(dat, alpha_level, cov, title){
  
  
  dat %>%
    mutate(cov_type == cov) %>%
    mutate(rho = as.factor(rho)) %>%
    filter(alpha == alpha_level) %>%
    mutate(beta = beta_1) %>%
    ggplot(aes(x = m, y = power_ratio, fill = rho)) + 
    geom_boxplot(alpha = .5) + 
    geom_hline(yintercept = 1, linetype = "solid") +
    #scale_y_continuous(breaks = seq(0, 1, .05)) + 
    scale_fill_manual(values = c("firebrick2", "firebrick4")) +
    ylim(c(0, NA)) +
    facet_grid(q ~ beta, scales = "free_y",  labeller = label_bquote(cols = beta == .(beta))) + 
    labs(x = "Number of Studies", y = "Power Ratio: HTZ/CWB", fill = expression(rho)) + 
    ggtitle(title) + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}

b <- create_power_rho_graph(dat = power_ratio, 
                            alpha_level = ".05", 
                            cov = "between", 
                            title = "Study-level covariate type")

w <- create_power_rho_graph(dat = power_ratio, 
                            alpha_level = ".05", 
                            cov = "within",
                            title = "Effect size-level covariate type")

b + w + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave("sim_results/graphs_paper/study_2/rho_power_052.png", device = "png", dpi = 500, height = 7, width = 12)


