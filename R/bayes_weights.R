# Compute Bayes p-value weights
#
# inputs:
# mu -  the estimated means, a vector of length J 
# sigma - the standard errors, a positive vector of length J 
# q - level at which tests will be performed
#
# Outputs:
# results - a list containing: 
# w -the optimal weights
# lambda - the dual certificate, normalizing constant produced by solving the optimization problem.
# q_star - true value of q solved for
# q_threshold - maximal value of q for which problem can be solve exactly
#
bayes_weights <- function(mu, sigma, q) {
  
  if (any(sigma<=0)) {warning("Positive variances required")}
  
  sigma_threshold <- 1e-2;
  if (any(sigma <= sigma_threshold)) {
    sigma <- sigma + sigma_threshold
    warning("Some variances too small. Adding a small constant to them")
  }
  
  if ((q <= 0) | (q >= 1)) {stop("Per-comparison error rate must be in (0,1)")}
  
  
  J <- length(mu)
  epsi <- 1e-3
  # define auxiliary variables 
  var_plus <- sigma^2 + 1
  var <- sigma^2
  mu2 <- mu^2
  alpha <- mu / var
  
  beta <- (mu2 + var * (mu2 + var_plus * log(var_plus))) / (var^2)
  gamma <- 2 * var_plus / var
  
  # define the functions used by Newton's method
  # the goal is to find the zero of the function f below
  f <- function(c) 1 / J * sum(pnorm(-alpha - sqrt(beta + gamma * c))) - q
  df <- function(c) 1 / J * sum(-dnorm(-alpha - sqrt(beta + gamma * c)) * gamma * (beta + gamma * c)^(- 1 / 2))
  q_threshold <- f(0)+q   
  
  if (f(0) >=0) { #Newton 
    x0 <- 0
    tol <- 1e-3 / J
    nmax <- 100
    x <- rep(0, nmax)
    ex <- rep(0, nmax)
    
    # Newton's method iterations
    x[1] <- x0 - f(x0) / df(x0)
    ex[1] <- abs(x[1] - x0)
    k <- 2
    
    while ((ex[k - 1] >= tol) && (k <= nmax)) {
      x[k] <- x[k-1] - (f(x[k - 1])/df(x[k - 1]))
      ex[k] <- abs(x[k] - x[k - 1])
      k <- k + 1
    }
    
    c <- x[k - 1] # could return this 
    
    w <- pnorm(-alpha - sqrt(beta + gamma * c)) / q
    q_star <- q
    lambda <- exp(c)
    
  } else { #Brent

    gamma <- sqrt(var_plus)
    overestimate_ind <- 0
    source("find_crossing.R")
    l_prime <- rep(0, J)
    for (i in 1:J) {
      l_prime[i] <- find_crossing(mu[i],gamma[i])
    }
    
    #find the interval where the dual variable lies
    l_sort <- sort(l_prime)
    ind <- order(l_prime)
    
    source("g_fun.R")
    #global dual function
    H <- function(lambda) sum(g_fun(mu, gamma, lambda, l_prime))
    
    l_sort  <- c(l_sort, 1)
    if (H(l_sort[J + 1]) <= J * q) {
      #binary search: l - lower u - upper
      l <- 1; u <- J + 1
      while (u - l > 1) {
        mid <- floor((u + l) / 2)
        if (H(l_sort[mid]) < J * q) {
          u <- mid
        } else {
          l <- mid
        }
      }
      #right interval is l,l+1
      a <- H(l_sort[l])
      c <- H(l_sort[l + 1])
      
      #check a <= J
      #this should never happen, because it corresponds to lambda > 1
      if (a > J) { cat("Error: dual constraint violated at 1")}
      
      #right limit of function at l_sort(l)
      b <- a - 1 + pnorm(c_1(mu[ind[l]],gamma[ind[l]],l_sort[l]))
      
      #should have a>b>c
      if ((a < b) | (b < c)){cat("Error: incorrect ordering of interval endpoints")}
      if ((J * q > a) | (J * q < c)){cat("Error: interval does not contain target")}
      
      #so we know c < b < a, and that c < Jq < a
      #two cases, depending on the ordering of Jq, b
      #(1)  Jq is in (b,a): in this case can't certify the original problem
      if (J * q > b) {
        lambda <- l_sort[l]
        # find closest endpoint
        dist_up <- a - J * q 
        dist_down <- J * q - b
        if (dist_up < dist_down) {
          q_star <- a # instead of Jq
          true_q <- a / J
        } else { 
          q_star <- b # instead of Jq
          true_q <- b/J
          overestimate <- a - b
          overestimate_ind <- ind[l]
        }
        #(2) : Jq is in (c,b): in this case can certify the original problem
      } else {
        # find the zero of the following function f
        f <- function(lambda) H(lambda) - J * q
        # will use standard matlab zero-finding need starting interval
        # ensure starting points have opposite sign: f is decreasing
        epsi <- 1e-7
        l_1 <- l_sort[l] + epsi
        l_2 <- l_sort[l + 1]
        while(f(l_1) < 0) {
          epsi <- epsi/10
          l_1 <- l_sort[l] + epsi
        }
        x0 <- c(l_1, l_2) #  initial interval
        
        lambda <- uniroot(f,x0)
        q_star <- J*q
        true_q <- q
      }
    }
    
    #calculate weights
    #note: dividing by true_q instead of q
    w_temp <-  g_fun(mu, gamma, lambda, l_prime)
    if (overestimate_ind > 0) {
      w_temp[overestimate_ind] <- w_temp[overestimate_ind] - overestimate
    }
    w <-  w_temp / true_q
  }
  
  # error check:
  if (abs(sum(w) - J) > epsi){ cat("Warning: Weights do not average to 1. Mean weight - 1 = ",mean(w) - 1) }
  if (any(w < 0)) {cat("Some weights are negative")}

  results <- list("w" = w, "q_star" = q_star,
"q_threshold" = q_threshold, "lambda" = lambda)
  return(results)
}
