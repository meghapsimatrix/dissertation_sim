library(tidyverse)
library(broom)
library(microbenchmark)
library(robumeta)
library(janitor)

load("data/meta_data_practice.RData")
source("homemade/tau/1_estimate_tau.R")
load("data/tsl_dat.RData")

tsl_dat <- tsl_dat %>%
  select(study, es_num, delta, v, dv, g2age) %>%
  drop_na()

meta_data$v <- meta_data$var_g


system.time(tau_sq_md <- calculate_tau_sq(meta_data, equation_full = "g ~ X1 + X2 + X3 + X4 + X5"))

system.time(calculate_tau_sq(tsl_dat, equation_full = "delta ~ g2age + dv"))

system.time(robu_comp_tsl <- robu(delta ~ g2age + dv, 
                                  studynum = study, 
                                  var.eff.size = v,
                                  small = FALSE,
                                  data = tsl_dat))


system.time(robu_comp_md <- robu(g ~ X1 + X2 + X3 + X4 + X5, 
                              studynum = study, 
                              var.eff.size = v,
                              small = FALSE,
                              data = meta_data))

robu_comp_md$mod_info$tau.sq
tau_sq_md


microbenchmark(
  calculate_tau_sq(meta_data, equation_full = "g ~ X1 + X2 + X3 + X4 + X5"),
  robu(g ~ X1 + X2 + X3 + X4 + X5, 
       studynum = study, 
       var.eff.size = v,
       small = FALSE,
       data = meta_data)
)
