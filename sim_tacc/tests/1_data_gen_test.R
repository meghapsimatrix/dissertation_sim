library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)


load("data/design_mat.Rdata")

#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_1.R")


set.seed(509985477)

meta_data <- 
  generate_rmeta(m = 10, 
                 tau = 0.1, 
                 rho = 0.8, 
                 covs = design_mat,
                 beta_type = "A")


# generate meta data 
set.seed(342020)


# some of the smds are really high lol
meta_data <- 
  generate_rmeta(m = 80, 
                 tau = 0.4, 
                 rho = 0.8, 
                 covs = design_mat,
                 beta_type = "A")

meta_data <- 
  generate_rmeta(m = 10, 
                 tau = 0.1, 
                 rho = 0.5, 
                 covs = design_mat,
                 beta_type = "A")

# generate meta data 
set.seed(342020)


# some of the smds are really high lol
meta_data_params <- 
  generate_rmeta(m = 80, 
                 tau = 0.4, 
                 rho = 0.8, 
                 covs = design_mat,
                 beta_type = "A", 
                 return_study_params = TRUE)

check <- meta_data %>%
  group_by(study) %>%
  summarize(n = n())

meta_data <- left_join(meta_data, meta_data_params %>% select(k, N) %>% mutate(study = 1:80))



big_meta <- 
  generate_rmeta(m = 1000, 
                 tau = 0.1, 
                 rho = 0.8, 
                 covs = design_mat,
                 beta_type = "B1",
                 return_study_params = FALSE)

mean(big_meta$var_g)
hist(big_meta$var_g)

# Study design features ---------------------------------------------------


# Use return_study_params to get 
# distribution of study design features

study_features <- 
  generate_rmeta(m = 1000, 
                 tau = 0.05, 
                 rho = 0.6, 
                 covs = design_mat,
                 beta_type = "A",
                 return_study_params = TRUE)

# check means of k, N, Psi

hist(study_features$k)
hist(study_features$N)
hist(study_features$Sigma)


sd(study_features$Sigma) # Should be sqrt(rho (1 - rho) / nu)
sqrt(0.6 * 0.4 / 50)



# check distribution of delta

delta_data <- 
  study_features %>%
  mutate(study = 1:n()) %>%
  unnest(cols = delta)


library(lme4)
lmer(delta ~ (1 | study), data = delta_data)
# Estimated SDs should be very close to tau and omega


# True parameter recovery -------------------------------------------------



library(clubSandwich)
library(metafor)

# Set nu = 4000 so that all Psi = rho

meta_data <- 
  generate_rmeta(m = 1000, 
                 tau = 0.05, 
                 rho = 0.6, 
                 covs = design_mat,
                 beta_type = "F1")

save(meta_data, file = "data/meta_data_practice.Rdata")

# Use clubSandwich::impute_covariance_matrix with true rho.

V_mat <- impute_covariance_matrix(vi = meta_data$var_g, 
                                  cluster = meta_data$study, 
                                  r = 0.6)


meta_data %>%
  group_by(study) %>%
 summarize(n = n()) %>% 
  View()

head(V_mat)


