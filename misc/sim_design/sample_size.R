library(purrr)
m <- c(10, 20 , 40, 80)

vals <- map(m, function(x) pmin(20+ 2 * rpois(x, 30), 200))

summary(unlist(vals))



big <- pmin(20 + 2 * rpois(100000, 30), 200)

summary(big)

80-72

130-80
80-40


hist(big)
