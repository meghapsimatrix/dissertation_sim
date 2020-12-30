library(robumeta)
library(clubSandwich)
library(wildmeta)
library(tidyverse)

set.seed(12102020)


full <- robu(d ~ study_type,
             studynum = study,
             var.eff.size = V,
             small = FALSE,
             data = SATcoaching)


cwb(full_model = full,
    indices = 2:3,
    R = 99)

Wald_test(full, constraints = constrain_equal(2:3), vcov = "CR2", test = "HTZ")



load("data/meta_data_practice.RData")

glimpse(meta_data)

full <- robu(g ~ X1 + X2 + X3 + X4 + X5,
             studynum = study,
             var.eff.size = var_g,
             small = FALSE,
             data = meta_data)


cwb(full_model = full,
    indices = 2:6,
    R = 399)

Wald_test(full, constraints = constrain_equal(2:6), vcov = "CR2", test = "HTZ")
