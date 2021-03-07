make_data <- function(cov){
  
  type1_dat %>%
    filter(alpha == ".05") %>%
    filter(test != "Naive-F") %>%
    #mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
    filter(str_detect(cov_test, cov)) %>%
    mutate(cov = cov)
  
}

covs <- c("X1", "X2", "X3", "X4", "X5")

all_dat_t1 <- map_dfr(covs, make_data)

glimpse(all_dat_t1)


all_dat_t1 %>%
  filter(alpha == ".05") %>%
  ggplot(aes(x = m, y = rej_rate, fill = test)) + 
  #geom_hline(yintercept = intercept, linetype = "solid") + 
  #geom_hline(yintercept = error, linetype = "dashed") + 
  geom_boxplot(alpha = .5) + 
  #scale_y_continuous(breaks = seq(0, .6, br)) + 
  scale_fill_brewer(palette = "Set1") +
  facet_grid(q ~ cov) + 
  labs(x = "Method", y = "Type 1 Error Rate") + 
  theme_bw() +
  theme(legend.position = "none",
        plot.caption=element_text(hjust = 0, size = 10))
