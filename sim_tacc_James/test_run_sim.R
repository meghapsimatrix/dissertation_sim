library(dplyr)
library(purrr)
library(mvtnorm)
library(clubSandwich)
library(tidyr)
library(stringr)
library(tibble)

# Tipton Pusto design matrix cleaned - clean_design_mat.R
load("data/design_mat.Rdata")
load("data/to_test.RData")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")
source("2_estimation_study_1.R")
source("3_performance_criteria.R")
source("4_run_sim_study.R")

#-----------------------------------------------------------
# Simulation Driver - should return a data.frame or tibble
#-----------------------------------------------------------

iterations <- 2
m <- 20
tau <- 0.2
rho <- 0.4
beta_type <- "A" 
batch <- 1
R <- 39
full_form = "X1 + X2 + X3 + X4 + X5"
test_dat = to_test
design_matrix = design_mat
Ftest_types = c("Naive-F","HTZ","EDT")
CWB_types = c("CWB","CWB-adj","CWB-fix","CWB-adj-fix")
seed = NULL

Ftest_types = c("EDT")
CWB_types = c("CWB-fix","CWB-adj-fix")


system.time(
  x <- run_sim(iterations, m, tau, rho, beta_type, batch, R, 
               Ftest_types = Ftest_types, CWB_types = CWB_types)  
)

x