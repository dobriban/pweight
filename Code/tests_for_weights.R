#This code tests the regularized weighting scheme on synthetic data

# generate means
J = 2000
mu = -abs(rnorm(J))
sigma = 1*rep(1,J)
alpha = 1/J

#find weights
w_1 = bayes_weights(mu,sigma, alpha)

plot(mu,w_1)

#vary alpha
for (j in (1:4)) {
  alpha = 2^j/J
  w_1 = bayes_weights(mu,sigma, alpha)
  plot(mu,w_1)
}
