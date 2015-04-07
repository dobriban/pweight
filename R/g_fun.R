# compute the dual constraint function: 1 if lambda <= l_prime, and a phi(c_1)
# otherwise
#
# Inputs (vectors of the same length)
# eta, gamma - prior mean and standard error 
# l_prime - crossing points
# lambda - the dual variable (a scalar)
#
# Output
# g_fun - value of the dual constraint (vector of same length)
g_fun <- function(eta, gamma, lambda, l_prime) {
  
  J  <- length(eta)
  g_fun <- rep(0, J)
  
  for (i in 1:J) {
    if (lambda <= l_prime[i]) {
      g_fun[i] <- 1
    } else{ 
      g_fun[i] <- pnorm(Re(c_1(eta[i], gamma[i], lambda)))
    }
  }  
  return(g_fun)
}