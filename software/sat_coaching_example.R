library(tidyverse)
library(clubSandwich)
library(wildmeta)
library(robumeta)

set.seed(12102020)

glimpse(SATcoaching)

full <- robu(d ~ study_type,
             studynum = study,
             var.eff.size = V,
             small = FALSE,
             data = SATcoaching)

full

Wald_test(full, constraints = constrain_equal(2:3), vcov = "CR2")

null <- robu(d ~ 1,
             studynum = study,
             var.eff.size = V,
             small = FALSE,
             data = SATcoaching)

system.time(res <- cwb(dat = SATcoaching,
    smd = d,
    var = V,
    cluster = study,
    full_model = full,
    null_model = null,
    indices = 2:3,
    R = 99))

res

library(metafor)

rma.mv()

SATcoaching$study_type <- as.factor(SATcoaching$study_type)
SATcoaching$study <- as.factor(SATcoaching$study)

mod_metafor <- rma.mv(yi = d,
                       V = V, 
                       mods = ~ study_type, 
                       data = SATcoaching)

Wald_test(mod_metafor, cluster = SATcoaching$study, 
          constraints = constrain_equal(2:3), vcov = "CR2")

class(mod_metafor)

class(full)


cwb
