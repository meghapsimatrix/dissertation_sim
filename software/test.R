load("/Users/meghajoshi/Box Sync/Dissertation_Joshi/simulation/data/meta_data_practice.RData")
library(wildmeta)
library(robumeta)
library(clubSandwich)

system.time(res <- cwb(dat = meta_data,
                       smd = g,
                       var = var_g,
                       cluster = study,
                       covs_full_form <- "X1 + X2 + X3 + X4 + X5", 
                       indices = 2:6))

system.time(res_adj <- cwb(dat = meta_data,
                       smd = g,
                       var = var_g,
                       cluster = study,
                       covs_full_form <- "X1 + X2 + X3 + X4 + X5", 
                       indices = 2:6,
                       adjust = TRUE))


full <- robu(g ~ X1 + X2 + X3 + X4 + X5,
             studynum = study,
             var.eff.size = var_g,
             small = FALSE,
             data = meta_data)

Wald_test(full, constraints = constrain_zero(2:6), vcov = "CR2")


cwb(dat = SATcoaching,
    smd = d,
    var = V,
    cluster = study,
    covs_full_form = "study_type",
    test_vars = "study_type",
    indices = 2:3,
    R = 99)


