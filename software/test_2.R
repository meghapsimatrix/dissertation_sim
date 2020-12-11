library(clubSandwich)
library(dplyr)
library(robumeta)

full <- robu(g ~ X1 + X2 + X3 + X4 + X5,
             studynum = study,
             var.eff.size = var_g,
             small = FALSE,
             data = meta_data)

null <- robu(g ~ 1,
             studynum = study,
             var.eff.size = var_g,
             small = FALSE,
             data = meta_data)

#
# data <- meta_data
# data$smd <- data$g
# data$var <- data$var_g
# data$cluster <- data$study
# null_model <- null
# full_model <- full
#

system.time(res <- cwb(dat = meta_data,
                       smd = g,
                       var = var_g,
                       cluster = study,
                       full_model = full,
                       null_model = null,
                       indices = 2:6,
                       R = 999,
                       adjust = TRUE))


library(clubSandwich)
library(robumeta)

full <- robu(d ~ study_type,
             studynum = study,
             var.eff.size = V,
             small = FALSE,
             data = SATcoaching)

null <- robu(d ~ 1,
             studynum = study,
             var.eff.size = V,
             small = FALSE,
             data = SATcoaching)

cwb(dat = SATcoaching,
    smd = d,
    var = V,
    cluster = study,
    full_model = full,
    null_model = null,
    indices = 2:3,
    R = 99)


full <- robu(g ~ X1 + X2 + X3 + X4 + X5,
             studynum = study,
             var.eff.size = var_g,
             small = FALSE,
             data = meta_data)

full_model <- full


stats::as.formula(paste("new_y ~ ",
                        stringr::str_split(as.character(full_model$ml), "~")[[3]]))

