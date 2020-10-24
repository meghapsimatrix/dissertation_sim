m <- 4
tau <- 2
rho <- 2
K <- 1000
betas_power <- 10
one_run <- 36/60

t_one_error <- (m * tau * rho * one_run * K) 

power <- (m * tau * rho * (one_run/2) * K) 
power_all <- power * betas_power

all <- t_one_error + power_all
# hours on my laptop single core
all

# in years 
all/(24 * 365)

# 500 SU's is 3 years I need like double that + for study 2

# SU's on tacc??
all / (68)  




