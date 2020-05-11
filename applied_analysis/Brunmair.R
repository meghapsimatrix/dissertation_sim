library(tidyverse)
library(janitor)


brunmair_dat <- read_csv("empirical_data/Brunmair and Richter 2019/Brunmair_data.csv") %>%
  clean_names()


num_es <- dat %>%
  group_by(name_year) %>%
  count()
