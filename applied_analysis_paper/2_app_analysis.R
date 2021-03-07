library(tidyverse)
library(broom)
library(robumeta)
library(janitor)
library(clubSandwich)
library(wildmeta)

load("data/tsl_dat_20.RData")


robu_comp_tsl <- robu(delta ~ g2age + dv, 
                      studynum = study, 
                      var.eff.size = v,
                      small = TRUE,
                      data = tsl_dat)


# Single coefficient ------------------------------------------------------

res_age_naive <- Wald_test(robu_comp_tsl, 
          constraints = constrain_zero(2), 
          vcov = "CR1",
          test = "Naive-F")

res_age_htz <- Wald_test(robu_comp_tsl, 
          constraints = constrain_zero(2), 
          vcov = "CR2",
          test = "HTZ")

res_age_edt <- Wald_test(robu_comp_tsl, 
                         constraints = constrain_zero(2), 
                         vcov = "CR2",
                         test = "EDT")

# multiple contrast hypothesis --------------------------------------------

res_mch_naive <- Wald_test(robu_comp_tsl, 
                            constraints = constrain_zero(3:7), 
                            vcov = "CR1", 
                            test = "Naive-F")

res_mch_htz <- Wald_test(robu_comp_tsl, 
                            constraints = constrain_zero(3:7), 
                            vcov = "CR2",
                            test = "HTZ")

res_mch_edt <- Wald_test(robu_comp_tsl, 
                         constraints = constrain_zero(3:7), 
                         vcov = "CR2",
                         test = "EDT")

      

# cluster wild ------------------------------------------------------------

set.seed(01302021)

full_mod <- robu(delta ~ g2age + dv, 
                 studynum = study, 
                 var.eff.size = v,
                 small = TRUE,
                 data = tsl_dat)

cwb_age <- cwb(full_mod, 
               indices = 2)

cwb_adj_age <- cwb(full_mod, 
               indices = 2, 
               adjust = TRUE)

cwb_mch <- cwb(full_mod, 
               indices = 3:7)

cwb_adj_mch <- cwb(full_mod, 
                   indices = 3:7, 
                   adjust = TRUE)


# bind results ------------------------------------------------------------

res_sc <- bind_rows(res_age_naive, res_age_edt, res_age_htz, cwb_age, cwb_adj_age) %>%
  select(Method = test, `F` = Fstat, delta, df_num, df_denom, p = p_val) %>%
  mutate_if(is.numeric, round, 3)

write_csv(res_sc, "applied_analysis_paper/res_sc.csv")


res_mch <- bind_rows(res_mch_naive, res_mch_edt, res_mch_htz, cwb_mch, cwb_adj_mch) %>%
  select(Method = test, `F` = Fstat, delta, df_num, df_denom, p = p_val) %>%
  mutate_if(is.numeric, round, 3)

write_csv(res_mch, "applied_analysis_paper/res_mch.csv")



