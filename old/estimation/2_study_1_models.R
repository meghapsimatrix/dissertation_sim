library(tidyverse)
library(robumeta)
library(clubSandwich)

# files -------------------------------------------------------------------
load("data/meta_data_practice.RData")

# combinations ------------------------------------------------------------

covs <- c("X1", "X2", "X3", "X4", "X5")

# terms combinations
comb_terms <- function(m, terms = covs) combn(length(terms), m, simplify = FALSE)

indices <- flatten(map(seq_along(1:5), comb_terms))
indices_test <- map(indices, function(x) x + 1)

contrasts <- map(indices_test, function(x) length(x))

equations <- map_chr(indices, function(x) paste(covs[x], collapse = "+"))

full_mod_indices <- 1:5
terms <- c("X1", "X2", "X3", "X4", "X5")

# null model 
null_indices <- map(indices, function(x) full_mod_indices[-x])
null_terms <- map(null_indices, function(x) covs[x])
null_model <- map_chr(null_terms, function(x) paste(x, collapse = "+"))


# all the null model to fit
to_test <- tibble(equ = equations,
                  type = c(rep("single", 5), rep("mch", 26)),
                  indices = indices) %>%
  mutate(null_indices = null_indices,
         null_terms = null_terms, 
         null_model= null_model, 
         null_model = if_else(null_model == "", "1", null_model),
         null_model = paste("g ~ ", null_model),
         indices_test = indices_test,
         contrasts = unlist(contrasts)) %>%
  select(cov_test = equ, null_model, indices_test, contrasts) 


save(to_test, file = "../sim_tacc/data/to_test.RData")
