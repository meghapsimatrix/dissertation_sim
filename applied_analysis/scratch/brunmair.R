library(tidyverse)
library(janitor)


brunmair_dat <- read_csv("empirical_data/Brunmair and Richter 2019/Brunmair_data.csv") %>%
  clean_names()


num_es <- brunmair_dat %>%
  group_by(study_id) %>%
  count() %>%
  ungroup()

table(num_es$n)

sum(num_es$n)
