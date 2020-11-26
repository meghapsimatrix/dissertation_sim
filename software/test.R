load("/Users/meghajoshi/Box Sync/Dissertation_Joshi/simulation/data/meta_data_practice.RData")
library(wildmeta)


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



system.time(res <- cwb(data = meta_data,
                       smd = g,
                       var = var_g,
                       cluster = study,
                       full_model = full,
                       null_model = null,
                       indices = 2:6,
                       R = 999))
