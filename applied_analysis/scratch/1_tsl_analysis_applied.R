library(tidyverse)
library(haven)
library(metafor)
library(robumeta)
library(clubSandwich)
library(labelled)
library(fastDummies)

# Read data ---------------------------------------------------------------

tsl_dat <- read_dta("/Users/meghajoshi/Box Sync/Dissertation_Joshi/data_melissa/empirical_data/Tanner_Smith_Lipsey_2015/Data for James.dta") %>%
  select(-referencenum, - refsummary) %>%
  haven::as_factor() %>% 
  mutate_if(is.factor, as.character) %>%
  filter(unit == "1. Individual") 

var_label(tsl_dat)

# res and then one focal moderator at a time
# unusual thing to do 
# don't need to do it
# simple- running looking at couple of moderators one var at a time analysis 
# one one at a time
# one - that has all of the important variable that they use 


# Check num of es per study -----------------------------------------------

es_num <- tsl_dat %>%
  group_by(studyid) %>%
  summarise(es = n_distinct(esid)) %>%
  ungroup()

es_num %>% 
  filter(es %in% c(min(es), max(es)))

summary(es_num$es)

# Check data--------------------------------------------------------------

# check the dv 
table(tsl_dat$dvcat)

# clean dv
tsl_dat <- tsl_dat %>%
  mutate(
    dv = str_sub(dvcat, 4),
    g1gtype = str_sub(g1gtype, 4)
  )

table(tsl_dat$dv)


# check missing data 
map(tsl_dat, ~sum(is.na(.)))


# Effect size and variance ------------------------------------------------

# convert se to sd? 
tsl_dat <- tsl_dat %>%
  mutate(
    sd_tx = replmiss(sd_tx, se_tx * sqrt(n_tx_ob)),
    sd_ct = replmiss(sd_ct, se_ct * sqrt(n_ct_ob))
  )

tsl_dat %>%
  select(esdata:tvaldep, es21) %>%
  View()


# use escalc to calculate delta and v - SMD is Hedges's g right? 
tsl_dat <- escalc(measure="SMD",
                  m1i = mean_tx, sd1i = sd_tx, n1i = n_tx_ob,
                  m2i = mean_ct, sd2i = sd_ct, n2i = n_ct_ob,
                  data = tsl_dat,
                  var.names = c("delta", "v")) 

map(tsl_dat, ~sum(is.na(.)))

tsl_na <- tsl_dat %>%
  filter(is.na(delta)) 

glimpse(tsl_na)


# some don't have sd or se or anything? 
tsl_dat <- tsl_dat %>%
  mutate(delta = replmiss(delta, es21 * -1),
         v = replmiss(v, ((n_tx_ob + n_ct_ob) / (n_tx_ob * n_ct_ob) + 
                            ((delta^2) / (2* (n_tx_ob + n_ct_ob))))))

map(tsl_dat, ~sum(is.na(.)))

tsl_dat_check <- tsl_dat %>%
  select(esdata:tvaldep, es21, delta, v)


save(tsl_dat, file = "data/tsl_dat.RData")

# Pooling ---------------------------------------------------

# Pooled es estimate 

fe_pooled <- rma(yi = delta, vi = v, data = tsl_dat, method = "FE")

summary(fe_pooled)
coef_test(fe_pooled, vcov = "CR2", cluster = tsl_dat$studyid)

re_pooled <- rma(yi = delta, vi = v, data = tsl_dat)

summary(re_pooled)
coef_test(re_pooled, vcov = "CR2", cluster = tsl_dat$studyid)


# rma.mv(yi = delta, vi = v, random = ~ 1 | studyid, data = tsl_dat)


# Analysis one moderator --------------------------------------------------

# correlated
rve_exam <- robu(delta ~ dv, 
                 studynum = studyid, 
                 var.eff.size = v,
                 small = TRUE,
                 data = tsl_dat)

rve_exam


# RVE without correction --------------------------------------------------

# Omnibus test using RVE with no small sample correction
# 2nd to 6th constrained to be 0 
# average for audit is the same as reference level
# average effect is same as ref category effect 
# simpler constrained model - there is an overall intercept 
# Null hypothesis - does not constrain overall average effect 
Wald_test(rve_exam, constraints = 2:6, vcov = "CR2", test = "Naive-F")

# Single coefficinet
Wald_test(rve_exam, constraints = 2, vcov = "CR2", test = "Naive-F")



# AHZ correction ----------------------------------------------------------

# Omnibus test using RVE and the HTZ small sample correction
Wald_test(rve_exam, constraints = 2:6, vcov = "CR2", test = "HTZ")

# Single coefficient
Wald_test(rve_exam, constraints = 2, vcov = "CR2", test = "HTZ")


# Analysis multiple moderator --------------------------------------------------

# two predictors 
# dv and g2perwhite


# correlated
rve_mult <- robu(delta ~ dv + g2perwhite + g2permale + g2age + g1txdur + g1numsessions +  g1delmode + g1txtype, 
                 studynum = studyid, 
                 var.eff.size = v,
                 small = TRUE,
                 data = tsl_dat)

rve_mult


