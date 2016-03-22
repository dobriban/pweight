# Unit tests for Exp weights

source("./R/exp_weights.R")

# Plots  ----------------
# generate means
J <- 2000
mu <- -abs(rnorm(J))
beta <- 4
q <- 5 #change to 5 if desired

# find weights
w_1 <- exp_weights(mu, beta, q)
plot(mu, w_1)
