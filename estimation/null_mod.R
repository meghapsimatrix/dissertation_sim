library(tidyverse)

load("data/meta_data_practice.RData")

vars <- names(meta_data %>% select(starts_with("X")))[2:6]

vars_new <- vars[!str_detect(vars, "X1")]

meta_equation <- paste(vars_new, collapse = " + ")
