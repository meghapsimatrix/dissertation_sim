library(tidyverse)


load("data/meta_data_practice.RData")

vars <- names(meta_data %>% select(starts_with("X")))[2:6]

vars_new <- vars[!str_detect(vars, "X1")]

meta_equation <- paste(vars_new, collapse = " + ")


dat <- tibble(cov = c("X1", "X2", "X3", "X4", "X5"),
       index = c(1:5))


terms <- dat$cov
comb_terms <- function(m, terms = dat$cov) combn(length(terms), m, simplify = FALSE)

indices <- unlist(map(seq_along(1:5), comb_terms), recursive = FALSE)
equations <- map_chr(indices, function(x) paste(terms[x], collapse = "+"))

full_mod_indices <- c(1, 2, 3, 4, 5)

null_indices <- map(indices, function(x) full_mod_indices[-x])
null_terms <- map(null_indices, function(x) terms[x])
null_model <- map_chr(null_terms, function(x) paste(x, collapse = "+"))







to_test <- tibble(equ = equations,
                  type = c(rep("single", 5), rep("mch", 26)),
                  indices = indices) %>%
  mutate(null_indices = null_indices,
         null_terms = null_terms, 
         null_model = null_model)

