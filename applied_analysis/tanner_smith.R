library(tidyverse)
library(haven)

# Read data ---------------------------------------------------------------

tsl_dat <- read_dta("empirical_data/Tanner_Smith_Lipsey_2015/Data for James.dta")

# Check num of es per study -----------------------------------------------
es_num <- tsl_dat %>%
  group_by(studyid) %>%
  summarise(es = n_distinct(esid)) %>%
  ungroup()


es_num %>% 
  filter(es %in% c(min(es), max(es)))

summary(es_num$es)



# Check the labels --------------------------------------------------------

vapply(tsl_dat, function(x) attributes(x)[["label"]], character(1))

table(tsl_dat$dvcat)


# Analysis ----------------------------------------------------------------

# correlated
rve_exam <- robu(smd ~ dvcat + design + g1perwhite + g1perblack + g1perhisp + g1peroth + g1permale + g1age + g1txdur + g1numsessions + g1manual + g1gttype, 
                 studynum = studyid, 
                 var.eff.size = v,
                 small = TRUE,
                 data = meta_data)
rve_exam
