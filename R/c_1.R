#Find the location of the smaller maximum of the generic term in the dual
#Inputs
# eta, gamma - prior mean and standard error
# lambda - the dual variable
# logscale - (optional) is lambda in log-scale? default = 0

c_1 <- function(eta, gamma, lambda, logscale = 0) {
  
  #not log-scale
  if (logscale == 0) {
    if (lambda < 0) {
      cat("lambda = ",lambda, " is negative")
      #error('c_1:neg_lambda', 'Lambda is not positive')
    }
    x <- log(lambda)    
  } else {    #  log-scale
    x <- lambda
  }
  
  square <- eta^2 + (gamma^2 - 1) * (eta^2 + 2 * gamma^2 * (log(gamma) + x))
  
  # numerical slack
  epsi <- 1e-7
  if (square < -epsi) {
    fprintf("square = ",square, " is negative")
    #error('c_1:neg_square', 'Quadratic expression q is not positive')
  }
  c1  <-  - (eta + sqrt(abs(square))) / (gamma^2 - 1)
}