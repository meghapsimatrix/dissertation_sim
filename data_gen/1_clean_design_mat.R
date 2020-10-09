library(tidyverse)
library(fastDummies)

# Design matrix -----------------------------------------------------------

design_mat <- read_csv("data/design_matrix.csv") %>%
  filter(!is.na(X5)) %>%
  mutate(X1 = as.numeric(X1), 
         X = 1) %>%
  select(X, everything())

# merge it on study level
#design_mat$X6 <- sample(c(rep("A", 8), rep("B", 8), rep("C", 4)))

#design_mat <- dummy_cols(design_mat, select_columns = "X6") %>%
#  select(-c(X6_A, X6))

glimpse(design_mat)

save(design_mat, file = "data/design_mat.Rdata")

load("data/design_mat.Rdata")



