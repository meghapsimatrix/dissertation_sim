library(tidyverse)
library(broom)
library(microbenchmark)
library(robumeta)
library(janitor)

load("data/meta_data_practice.RData")

source("estimation/estimate_tau.R")


f1 <- function(dat) {
  dat %>%
    drop_na() %>%
    group_by(study) %>%
    mutate(k_j = n_distinct(es_num),
           sigma_sq_j = mean(v),
           w_tilde_j = 1/(k_j * sigma_sq_j)) %>%
    ungroup()
}


f2 <- function(dat) {
  
  dat <- drop_na(dat)
  
  k_j <- as.numeric(table(dat$study))
  sigma_sq_j <- tapply(dat$v, dat$study, mean)
  w_tilde <- 1 / (k_j * sigma_sq_j)
  
  dat$w_tilde_j <- w_tilde[dat$study]

  dat  
}

f1(meta_data)
f2(meta_data)

microbenchmark(f1 = f1(meta_data), f2 = f2(meta_data))
