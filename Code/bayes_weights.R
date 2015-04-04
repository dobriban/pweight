# Edgar Dobriban
# this code computes the regularized weights for p-value weighting
# of Owen & Dobriban;
# given estimated means and standard errors of test statistics, the weighting scheme optimizes
# the expected number of discoveries at some specific level q
#
# inputs:
# mu - a non-positive vector of length J, the estimated means of test statistics
# sigma - a positive vector of length J, the estimated standard errors of test statistics
# q - level at which the number of discoveries should be maximized; 
#        J*q is the number of expected false rejections, typically a constant
#        even if the problem size grows

regularized_weights = function(mu, sigma, q) {
  
  J = length(mu)
  
  # define auxiliary variables 
  var_plus = sigma^2+1
  var = sigma^2
  mu2 = mu^2
  alpha = mu/var
  
  beta = (mu2 + var*(mu2+var_plus*log(var_plus)))/(var^2)
  gamma = 2*var_plus/var
  
  #define the functions used by Newton's method
  #the goal is to find the zero of the function f below
  f <- function(c) 1/J*sum(pnorm(-alpha - sqrt(beta + gamma*c)))-q;
  df <- function(c) 1/J*sum(-dnorm(-alpha - sqrt(beta + gamma*c))*gamma*(beta + gamma*c)^(-1/2));
  

  x0 = 0;
  tol = 1e-3/J;
  nmax = 100;
  x = rep(0,nmax)
  ex = rep(0,nmax)
  
  #Newton's method iterations
  x[1] = x0 - (f(x0)/df(x0));
  ex[1] = abs(x[1]-x0);
  k = 2;
  
  while ((ex[k-1] >= tol) && (k <= nmax)) {
  x[k] = x[k-1] - (f(x[k-1])/df(x[k-1]));
  ex[k] = abs(x[k]-x[k-1]);
  k = k+1;
  }
  
  c = x[k-1]; #could return this 
  
  w = pnorm(-alpha - sqrt(beta + gamma*c))/q;
  
  return(w)
}
