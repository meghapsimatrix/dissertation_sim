power_ratio %>%
  filter(alpha == ".05") %>%
  mutate(m = as.character(m)) %>%
  mutate(beta = ifelse(str_detect(beta_type, "1"), .1, .5)) %>%
  ggplot(aes(x = m, y = power_ratio, color = cov_test, group = 1)) + 
  geom_line(alpha = .5) + 
  geom_hline(yintercept = 1, linetype = "solid") +
  #scale_y_continuous(breaks = seq(0, 1, .05)) + 
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(~beta, scales = "free_y",  labeller = label_bquote(cols = beta == .(beta))) + 
  labs(x = "Number of Studies", y = "Power Ratio: HTZ/ CWB") + 
  theme_bw()
