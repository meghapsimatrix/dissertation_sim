library(tidyverse)


load("sim_results/sim_test_1_2_3.RData")
load("sim_tacc/data/to_test.RData")

results %>%
  filter(K < 50)


results <- left_join(results, 
                     to_test %>% select(cov_test, contrasts),
                     by = "cov_test")

small_res <- results %>%
  group_by(m, tau, rho, cov_test, beta_type, contrasts, test) %>%
  summarize(rej_rate_05 = mean(rej_rate_05),
            mcse_05 = mean(mcse_05),
            .groups = "drop")


small_res %>%
  filter(beta_type == "A") %>%
  mutate(m = paste("m =", m),
         q = paste("q =", contrasts))  %>%
  ggplot(aes(x = test, y = rej_rate_05, fill = test)) + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  geom_boxplot(alpha = .5) + 
  facet_grid(q ~ m, scales = "free") + 
  labs(x = "Method", y = "Type 1 Error Rate") + 
  theme_bw() +
  theme(legend.position = "none",
        plot.caption=element_text(hjust = 0, size = 10))



small_res %>%
  filter(beta_type == "B5", str_detect(cov_test, "X1")) %>%
  mutate(m = paste("m =", m),
         q = paste("q =", contrasts))  %>%
  ggplot(aes(x = test, y = rej_rate_05, fill = test)) + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  geom_boxplot(alpha = .5) + 
  facet_grid(q ~ m, scales = "free") + 
  labs(x = "Method", y = "Type 1 Error Rate") + 
  theme_bw() +
  theme(legend.position = "none",
        plot.caption=element_text(hjust = 0, size = 10))



small_res %>%
  filter(beta_type == "C5", str_detect(cov_test, "X2")) %>%
  mutate(m = paste("m =", m),
         q = paste("q =", contrasts))  %>%
  ggplot(aes(x = test, y = rej_rate_05, fill = test)) + 
  geom_hline(yintercept = .8, linetype = "dashed") + 
  geom_boxplot(alpha = .5) + 
  facet_grid(q ~ m, scales = "free") + 
  labs(x = "Method", y = "Power X2 C5") + 
  theme_bw() +
  theme(legend.position = "none",
        plot.caption=element_text(hjust = 0, size = 10))
