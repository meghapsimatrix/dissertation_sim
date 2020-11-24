library(tidyverse)

set.seed(20201123)

# study level 
X1_levels <- LETTERS[1:5]
X1 <- sample(X1_levels, 20, replace=TRUE, prob = c(0.2, 0.1, 0.3, 0.1, .3))

design_mat_bet <- tibble(X1 = rep(X1, each = 10),
                         study = rep(1:20, each = 10),
                         X = 1) %>%
  select(study, X, X1)

design_mat_bet

# es level 
# is it necessary that all studies have all the categories?
# is there a better way to do this one?
X2 <- sample(X1_levels, 200, replace=TRUE, prob = c(0.2, 0.1, 0.3, 0.1, .3))

design_mat_within <- tibble(X2 = X2,
                            study = rep(1:20, each = 10),
                            X = 1) %>%
  select(study, X, X2)

design_mat_within

table(design_mat_within$X2, design_mat_within$study)


# dummy code this, and set intercept to 0 and then 
# set X1_A & X2_A regression coefficient as .5 or  something for power? 
