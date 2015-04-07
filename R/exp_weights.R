# Compute Exponential p-value weights
#
# inputs:
# mu -  the estimated means, a vector of length J 
# beta - weights are proportional to exp(mu*beta)
# q - level at which tests will be performed
#
# Outputs:
# w -the weights
#
exp_weights <- function(mu,  beta,  q) {
  
  J <-  length(mu)
  s <-  rep(0,J)
  
  mu_s <-  sort(mu)
  sort_index <-  order(mu)
  u <-  exp(beta * mu_s)
  c <-  mean(u)
  S  <-  sum(s)
  w <-  u / c * (J - S / q) / (J - S)
  ind <- w > 1 / q
  
  surplus <-  sum(w[ind]) - length(w[ind]) / q
  w[ind] <-  1 / q
  k <-  length(w[ind])
  
  while (surplus > 0) {
    increment <-  min(1 / q - w[J - k], surplus)
    w[J - k] <-  w[J - k]  +  increment
    surplus <-  surplus  -  increment
    k <-  k + 1
  }
  
  v <-  rep(0,J)
  for (i in 1:J){
    v[i] <-  w[sort_index == i]  
  }
  w <-  v
  return(w)
}