load("data/tsl_dat.RData")

tsl_dat <- tsl_dat %>%
  select(study, es_num, delta, v, dv, g2age) %>%
  drop_na() %>%
  mutate(dv = factor(dv),
         dv = fct_recode(dv, bac = "Blood alcohol concentration",
                         comb = "Combined measures (e.g., AUDIT)",
                         fhu = "Frequency of heavy use",
                         fu = "Frequency of use",
                         pc = "Peak consumption",
                         qu = "Quantity of use"))


random_20 <- sample(unique(tsl_dat$study), 20)


tsl_dat <- tsl_dat %>%
  filter(study %in% random_20)

tsl_dat %>%
  group_by(study) %>%
  summarize(es = n_distinct(es_num)) %>%
  ungroup() %>%
  summarize(range(es))

save(tsl_dat, file = "data/tsl_dat_20.RData")
