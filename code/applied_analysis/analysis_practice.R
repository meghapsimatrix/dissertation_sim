library(clubSandwich)
library(robumeta)
library(metafor)
library(tidyverse)

# Need to turn this into a function/ functions
# write this up - empirical examples 
# readers see how methods would be applied
# include the code - readers to template 
# style - make the code as clear and accessible to your audience 

# matrix math version confuse people 
# running the simulation - start by doing the canned routines - using built in functions when you can 
# then you discover that your code takes too long to run - profiling - bottle necks - speed that up
# Wald test might be bottleneck

load("data/meta_data_practice.Rdata")

# correlated
rve_exam <- robu(smd ~ X1 + X2 + X3 + X4 + X5, 
                 studynum = study, 
                 var.eff.size = v,
                 small = TRUE,
                 data = meta_data)
rve_exam

# rma.mv - different working model, rma.mv restricted max likelihood, robu moment estimator (crappy)
# weights exactly inverse var; robu takes shortcuts - can actually fail in werid ways - certain situations
# if you have a covariate varies within study weird results 
# what should be done - rma.mv and clubSandwich on top of that 
# current state - people are using robu

# RVE without correction --------------------------------------------------

# Omnibus test using RVE with no small sample correction
Wald_test(rve_exam, constraints = 2:6, vcov = "CR2", test = "Naive-F")

# Single coefficinet
Wald_test(rve_exam, constraints = 2, vcov = "CR2", test = "Naive-F")



# AHZ correction ----------------------------------------------------------

# Omnibus test using RVE and the HTZ small sample correction
Wald_test(rve_exam, constraints = 2:6, vcov = "CR2", test = "HTZ")

# Single coefficient
Wald_test(rve_exam, constraints = 2, vcov = "CR2", test = "HTZ")


# Cluster Wild Bootstrapping ----------------------------------------------

# how to even start? 

# divide coefficients we care about and coefficients we don't and then do some matrix stuff? 
# set up a generic constrast
# change the contrast to test the different hypothesis 
# tricky 
# generic function 
# bootstrap null - full to null model find a way to do that 
# robu - if x3 to x5 null 
