library(tidyverse)

# Design matrix -----------------------------------------------------------

design_mat <- read_csv("data/design_matrix.csv") %>%
  filter(!is.na(X5)) %>%
  mutate(X1 = as.numeric(X1), 
         X = 1) %>%
  select(X, everything())

glimpse(design_mat)

save(design_mat, file = "data/design_mat.Rdata")
