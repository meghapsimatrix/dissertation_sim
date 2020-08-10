library(tidyverse)
library(broom)
library(microbenchmark)
library(robumeta)
library(janitor)

load("data/meta_data_practice.RData")
source("estimation/1_estimate_tau.R")
load("data/tsl_dat.RData")

tsl_dat <- tsl_dat %>%
  select(studyid, esid, delta, v, dv, g2age) %>%
  rename(study = studyid, es_num = esid) %>%
  drop_na()

system.time(calculate_tau_sq(meta_data, equation_full = "smd ~ X1 + X2 + X3 + X4 + X5"))

system.time(calculate_tau_sq(tsl_dat, equation_full = "delta ~ g2age + dv"))

system.time(robu_comp_tsl <- robu(delta ~ g2age + dv, 
                                  studynum = study, 
                                  var.eff.size = v,
                                  small = FALSE,
                                  data = tsl_dat))


system.time(robu_comp_md <- robu(smd ~ X1 + X2 + X3 + X4 + X5, 
                              studynum = study, 
                              var.eff.size = v,
                              small = FALSE,
                              data = meta_data))


microbenchmark(
  calculate_tau_sq(meta_data, equation_full = "smd ~ X1 + X2 + X3 + X4 + X5"),
  robu(smd ~ X1 + X2 + X3 + X4 + X5, 
       studynum = study, 
       var.eff.size = v,
       small = FALSE,
       data = meta_data)
)
