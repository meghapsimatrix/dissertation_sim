library(tidyverse)

load("../data/sim_test_1104_3.Rdata")

results %>%
  filter(K < 2) %>%
  View()




load("../data/sim_test_1104_2.Rdata")

results %>%
  filter(K < 2) %>%
  View()


load("../data/sim_test_1104.Rdata")

results %>%
  filter(K < 2) %>%
  View()
