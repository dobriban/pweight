#' informed Genome-Wide Association Study
#'
#'  For more details, see the paper "Optimal Multiple Testing Under a
#'  Gaussian Prior on the Effect Sizes", by Dobriban, Fortney, Kim and Owen,
#'   \url{http://arxiv.org/abs/1504.02935}
#'
#'@param P_current P-values in the current study, a numeric vector of length J, with entries between 0 and 1
#'@param N_current sample size in the current study, a positive integer (or vector of length J)
#'@param P_prior P-values in the prior study, a numeric vector of length J, with entries between 0 and 1
#'@param N_prior sample size in the current study, a positive integer (or vector of length J)
#' @param q (optional) uncorrected level at which tests should be performed. Default \code{q = 0.05}
#' @param phi (optional) dispersion factor used to multiply all standard errors. Default \code{phi = 1}
#'
#'@return A vector of 0-1s indicating the significant tests (1-s)
#'
#'@export
#'
iGWAS <- function(P_current, N_current, P_prior, N_prior, q=0.05, phi = 1) {

  #Error checking: stop if the variables are not in range
  if (any(P_current > 1) | any(P_current < 0) | any(P_prior > 1) | any(P_prior < 0)) {
    stop("P-values must be between 0 and 1")
  }
  if (any(N_current < 1) | any(N_prior < 1)) {
    stop("Sample sizes must be at least 1")
  }
  if ((q <= 0) |
      (q >= 1)) {
    stop("Level at which tests will be performed must be in (0,1)")
  }
  if (phi <= 0) {
    stop("Dispersion parameter must be positive")
  }

  #Define auxiliary variables to compute weights
  J <- length(P_current)
  if (length(N_current)==1) {
    N_current <- N_current *rep(1,J)
  }
  if (length(N_prior)==1) {
    N_prior <- N_prior *rep(1,J)
  }
  T_prior <- qnorm(P_prior/2)
  eta <- T_prior*sqrt(N_current/N_prior)
  sigma <- phi^2*sqrt(N_current/N_prior)

  #Compute weights
  res <- bayes_weights(mu, sigma, q)
  w <- res$w

  #Perform weighted multiple testing
  P_weighted <- P_current/w
  P_w_adjusted <- p.adjust(P_weighted,"bonferroni")
  ind_w <- (P_w_adjusted<q)
  cat(c("Number of significant tests using Weighting: ", sum(ind_w) ))

  return(ind_w)
}
