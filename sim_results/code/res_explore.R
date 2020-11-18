library(tidyverse)

load("sim_tacc/data/to_test.RData")


# load results ------------------------------------------------------------

files <- list.files("sim_results/results", full.names = TRUE)

load_res <- function(file){
  
  load(file)
  results
  
}


results <- map_dfr(files, load_res)



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

data_int <- tibble(alpha = c(".01", 
                             ".05", 
                             ".10"),
                   int = c(.01, .05, .10))

naive_dat <- results %>%
  select(test, beta_type, rho, tau, m, contrasts, starts_with("rej_rate")) %>%
  filter(test == "Naive-F", beta_type == "A") %>%
  mutate(m = as.character(m),
         q = paste("q =", contrasts)) %>%
  gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
  mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                           alpha == "rej_rate_05" ~ ".05",
                           alpha == "rej_rate_10" ~ ".10"))  %>%
  group_by(m, rho, tau, q, alpha) %>%
  summarize_at(vars(rej_rate), mean) %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))

summary(naive_dat$mcse)


naive_dat %>%
  ggplot(aes(x = m, y = rej_rate, fill = m)) +
    geom_hline(data = data_int, aes(yintercept = int), linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, 1, .1)) + 
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(alpha ~ q, scales = "free_y",  labeller = label_bquote(rows = alpha == .(alpha))) + 
    labs(x = "Number of Studies", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
  

ggsave("sim_results/graphs/naivef.png", device = "png", dpi = 500, height = 7, width = 12)



# Rej Rate Mean -----------------------------------------------------------

res_small <- results %>%
  mutate(m = paste("m =", m),
         q = paste("q =", contrasts)) %>%
  gather(alpha, rej_rate, c(rej_rate_01, rej_rate_05, rej_rate_10)) %>%
  mutate(alpha = case_when(alpha == "rej_rate_01" ~ ".01",
                           alpha == "rej_rate_05" ~ ".05",
                           alpha == "rej_rate_10" ~ ".10")) %>%
  group_by(m, tau, rho, beta_type, q, test, alpha) %>%
  summarize(rej_rate = mean(rej_rate),
            .groups = "drop") %>%
  mutate(mcse = sqrt((rej_rate * (1 - rej_rate))/ K_all))

# Alpha .05 

# Type 1 error ------------------------------------------------------------


create_type1_graph <- function(dat, intercept, br){
  
  dat <- dat %>%
    filter(test != "Naive-F",
           beta_type == "A")
  
  dat %>%
    ggplot(aes(x = test, y = rej_rate, fill = test)) + 
    geom_hline(yintercept = intercept, linetype = "dashed") + 
    geom_boxplot(alpha = .5) + 
    scale_y_continuous(breaks = seq(0, .6, br)) + 
    scale_fill_brewer(palette = "Set1") +
    facet_grid(q ~ m) + 
    labs(x = "Method", y = "Type 1 Error Rate") + 
    theme_bw() +
    theme(legend.position = "none",
          plot.caption=element_text(hjust = 0, size = 10))
}



create_type1_graph(dat = res_small %>% filter(alpha == ".05"), intercept = .05, br = .02)



ggsave("sim_results/graphs/type1.png", device = "png", dpi = 500, height = 7, width = 12)

# Type 1 error ------------------------------------------------------------
# 01

create_type1_graph(dat = res_small %>% filter(alpha == ".01"), intercept = .01, br = .01)

ggsave("sim_results/graphs/type1_01.png", device = "png", dpi = 500, height = 7, width = 12)

# Type 1 error ------------------------------------------------------------
# 10


create_type1_graph(dat = res_small %>% filter(alpha == ".10"), intercept = .10, br = .02)

ggsave("sim_results/graphs/type1_10.png", device = "png", dpi = 500, height = 7, width = 12)




# Power Difference --------------------------------------------------------

power <- res_small %>%
  filter(test %in% c("CWB", "HTZ"), beta_type != "A") %>%
  select(-c(starts_with("mcse"))) %>%
  spread(test, rej_rate) %>%
  mutate(power_diff = CWB - HTZ,
         power_ratio = CWB / HTZ) %>%
  group_by(m, rho, tau, alpha, q) %>%
  summarize_at(vars(power_diff), mean)


power %>%
  ggplot(aes(x = m, y = power_diff, fill = m)) + 
  geom_boxplot(alpha = .5) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  #scale_y_continuous(breaks = seq(0, 1, .05)) + 
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(alpha ~ q, scales = "free_y",  labeller = label_bquote(rows = alpha == .(alpha))) + 
  labs(x = "Number of Studies", y = "Difference in Power: CWB - HTZ") + 
  theme_bw() +
  theme(legend.position = "none",
        plot.caption=element_text(hjust = 0, size = 10))


ggsave("sim_results/graphs/power.png", device = "png", dpi = 500, height = 7, width = 12)





