# Unit tests for Exp weights

setwd("C:/Dropbox/Weighted New/pvalue_weighting_r/R")

# Plots  ----------------
# generate means
J <- 2000
mu <- -abs(rnorm(J))
beta <- 4
q <- 0.5

# find weights
source("exp_weights.R")
w_1 <- exp_weights(mu, beta, q)
plot(mu, w_1)
