library(tidyverse)
library(fastDummies)

# Design matrix -----------------------------------------------------------

design_mat <- read_csv("data/design_matrix.csv") %>%
  filter(!is.na(X5)) %>%
  mutate(X1 = as.numeric(X1), 
         X = 1) %>%
  select(X, everything())


sample(c(rep("A", 9), rep("B", 9), rep("C", 2)))

design_mat$X6 <- rep(sample(c("A", "B", "C"), size = 20, replace = TRUE), each = 10)

design_mat <- dummy_cols(design_mat, select_columns = "X6") %>%
  select(-c(X6_A, X6))

glimpse(design_mat)

save(design_mat, file = "data/design_mat.Rdata")



