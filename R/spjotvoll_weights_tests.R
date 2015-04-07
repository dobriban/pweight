# Unit tests for Spjotvoll weights

setwd("C:/Dropbox/Weighted New/pvalue_weighting_r/R")
source("bayes_weights.R")
source("spjotvoll_weights.R")

# Plots  ----------------
# generate means
J <- 2000
mu <- -abs(rnorm(J))
sigma <- 0.1 * rep(1, J)
alpha <- 1 / J

# find weights
w_1 <- bayes_weights(mu, sigma, alpha)
plot(mu, w_1$w)
w_2 <- spjotvoll_weights(mu, alpha)
points(mu, w_2,col="red", pch=2)
