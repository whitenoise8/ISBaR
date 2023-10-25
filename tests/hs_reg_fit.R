#' hs_reg_fit.R
#' author: Cristian Castiglione
#' creation: 12/09/2023
#' last change: 12/09/2023
#' 
#' description: 
#'   This file implements a Gibbs sampling algorithm for estimating a Bayesian 
#'   linear model predicting the response y given the covariate matrix X.
#'   Assuming that only a small fraction of the columns of X are relevant for 
#'   predicting y, we consider a Horseshoe prior over beta.
#'   The implementation here proposed is just a wrapper of the "horseshoe"
#'   function from the "horseshoe" package.
#'   
#' references:
#'   Bhattacharya A., Chakraborty A., and Mallick B.K (2016), 
#'   Fast sampling with Gaussian scale-mixture priors in high-dimensional regression. 
#'   Biometrika 103(4), 985–991.
#' 
#'   Polson, N.G., Scott, J.G. and Windle, J. (2014) 
#'   The Bayesian Bridge. 
#'   Journal of Royal Statistical Society, B, 76(4), 713-733.
#' 
#'   Rue, H. (2001). 
#'   Fast sampling of Gaussian Markov random fields. 
#'   Journal of the Royal Statistical Society: Series B, 63, 325–338.
#' 
#'   Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010), 
#'   The Horseshoe Estimator for Sparse Signals. 
#'   Biometrika 97(2), 465–480.
#'   
#' model specification:
#' 
#'   y = X beta + e, e ~ N(0, sigma2)
#'   beta_j ~ N(0, sigma2 tau^2 lambda_j^2), (j = 1, ..., p)
#'   lambda_j ~ HC(0, 1), (j = 1, ..., p)
#'   tau ~ HC(0, 1),
#'   sigma2 ~ 1 / sigma2 (Jeffreys prior)
#'   

hs.reg.fit = function (y, X, Z, prior = list(), control = list()) {
  
  require(horseshoe)
  
  time0 = proc.time()
  
  # Sanity check for the control parameters
  check = function (object) {
    flag = FALSE
    if (!is.null(object)) {
      if (is.integer(object)) {
        if (object > 0) {
          flag = TRUE
        }
      }
    }
    return (flag)
  }
  
  default = list(burn = 1000, niter = 5000, thin = 5, verbose = TRUE)
  
  if (check(control$burn)) default$burn = control$burn
  if (check(control$niter)) default$niter = control$niter
  if (check(control$thin)) default$thin = control$thin
  if (is.logical(control$verbose)) default$verbose = control$verbose
  
  # Fit the Horseshoe linear model
  fit = NULL
  if (default$verbose) {
    fit <- horseshoe::horseshoe(y = y, X = cbind(X, Z),
                                method.tau = "halfCauchy",
                                method.sigma = "Jeffreys",
                                burn = default$burn, 
                                nmc = default$niter,
                                thin = default$thin)
  } else {
    invisible(
      capture.output(
        fit <- horseshoe::horseshoe(y = y, X = cbind(X, Z),
                                    method.tau = "halfCauchy",
                                    method.sigma = "Jeffreys",
                                    burn = default$burn, 
                                    nmc = default$niter,
                                    thin = default$thin)
      )
    )
  }
  
  
  # Get the posterior samples
  trace = list(
    beta = t(fit$BetaSamples),
    tau = fit$TauSamples,
    sigma2 = fit$Sigma2Samples)
  
  timef = proc.time()
  
  # Return the posterior model fit
  list(y = y, X = X, Z = Z, 
       prior = prior, control = default, 
       trace = trace, exe.time = (timef - time0)[3])
}
