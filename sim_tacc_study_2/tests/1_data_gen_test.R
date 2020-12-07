library(dplyr)
library(purrr)
library(mvtnorm)
library(robumeta)
library(clubSandwich)
library(tidyr)
library(fastDummies)


#-----------------------------------------------------------
# Source the functions
#-----------------------------------------------------------

source("1_data_gen_study_2.R")


set.seed(509985477)

m <- 10
tau <- 0.1
rho <- 0.8
cov_type <- "between"
cat_num <- 3
beta_type <- "A"
beta_type = "A"
k_mean = 4
N_mean = 30
nu = 50


meta_data <- 
  generate_rmeta(m = 10, 
                 tau = 0.1, 
                 rho = 0.8, 
                 cov_type = "between",
                 cat_num = 3,
                 beta_type = "A")

meta_data_w <- 
  generate_rmeta(m = 10, 
                 tau = 0.1, 
                 rho = 0.8, 
                 cat_num = 5,
                 cov_type = "within",
                 beta_type = "A")




# generate meta data 
set.seed(342020)


# some of the smds are really high lol
meta_data <- 
  generate_rmeta(m = 80, 
                 tau = 0.4, 
                 rho = 0.8, 
                 cov_type = "between",
                 cat_num = 5,
                 beta_type = "A")

save(meta_data, file = "../data/meta_data_practice_2.Rdata")



# generate meta data 
set.seed(342020)


# some of the smds are really high lol
meta_data_params <- 
  generate_rmeta(m = 80, 
                 tau = 0.4, 
                 rho = 0.8, 
                 cov_type = "between",
                 cat_num = 5, 
                 beta_type = "B5", 
                 return_study_params = TRUE)

check <- meta_data %>%
  group_by(study) %>%
  summarize(n = n())

meta_data <- left_join(meta_data, 
                       meta_data_params %>% select(k, N) %>% 
                         mutate(study = 1:80))



big_meta <- 
  generate_rmeta(m = 1000, 
                 tau = 0.1, 
                 rho = 0.8, 
                 cov_type = "between",
                 cat_num = 4,
                 beta_type = "B5",
                 return_study_params = FALSE)

save(big_meta, file = "../data/big_meta.Rdata")


mean(big_meta$var_g)
hist(big_meta$var_g)

glimpse(big_meta)

# Study design features ---------------------------------------------------


# Use return_study_params to get 
# distribution of study design features

study_features <- 
  generate_rmeta(m = 1000, 
                 tau = 0.05, 
                 rho = 0.6, 
                 cov_type = "between",
                 cat_num = 5,
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
                 cov_type = "between",
                 cat_num = 5,
                 beta_type = "B5")


# Use clubSandwich::impute_covariance_matrix with true rho.

V_mat <- impute_covariance_matrix(vi = meta_data$var_g, 
                                  cluster = meta_data$study, 
                                  r = 0.6)

head(V_mat)

meta_data %>%
  group_by(study) %>%
  summarize(n = n()) %>% 
  View()


