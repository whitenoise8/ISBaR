# file: hs_reg_fit.R
# author: Cristian Castiglione
# creation: 12/09/2023
# last change: 12/09/2023
# 
# description: 
#   This file implements a Gibbs sampling algorithm for estimating a Bayesian 
#   linear model predicting the response y given the covariate matrix X.
#   Assuming that only a small fraction of the columns of X are relevant for 
#   predicting y, we consider a Horseshoe prior over beta.
#   The implementation here proposed is just a wrapper of the "horseshoe"
#   function from the "horseshoe" package.
#   
# references:
#   Bhattacharya A., Chakraborty A., and Mallick B.K (2016), 
#   Fast sampling with Gaussian scale-mixture priors in high-dimensional regression. 
#   Biometrika 103(4), 985–991.
# 
#   Polson, N.G., Scott, J.G. and Windle, J. (2014) 
#   The Bayesian Bridge. 
#   Journal of Royal Statistical Society, B, 76(4), 713-733.
# 
#   Rue, H. (2001). 
#   Fast sampling of Gaussian Markov random fields. 
#   Journal of the Royal Statistical Society: Series B, 63, 325–338.
# 
#   Carvalho, C. M., Polson, N. G., and Scott, J. G. (2010), 
#   The Horseshoe Estimator for Sparse Signals. 
#   Biometrika 97(2), 465–480.
#   
# model specification:
# 
#   y = X beta + e, e ~ N(0, sigma2)
#   beta_j ~ N(0, sigma2 tau^2 lambda_j^2), (j = 1, ..., p)
#   lambda_j ~ HC(0, 1), (j = 1, ..., p)
#   tau ~ HC(0, 1),
#   sigma2 ~ 1 / sigma2 (Jeffreys prior)
#   


#' @title Fit a Bayesian linear model with Horseshoe prior via Gibbs sampling
#' 
#' @description
#' \code{hs.reg.fit} is a wrapper of the \code{horseshoe} function from the \code{horseshoe} package.
#' 
#' @param y response vector
#' @param X fixed effect design matrix, including all the variables that must not be penalized
#' @param Z random effect design matrix, including all the variable subject to regularization
#' @param prior list of prior parameters. More information are provided in 'Details'.
#' @param control list of control parameter. More information are provided in 'Details'.
#' 
#' @details
#' Additional details...
#' 
#' 
#' @return Returns an object of class "hs", which is a list containing the following elements:
#' \describe{
#'   \item{\code{y}}{response vector}
#'   \item{\code{X}}{fixed effect design matrix}
#'   \item{\code{Z}}{random effect design matrix}
#'   \item{\code{prior}}{list of prior parameters}
#'   \item{\code{control}}{list of control parameters}
#'   \item{\code{burn}}{list containig the sampled chains of all the parameters during the burn-in iterations}
#'   \item{\code{trace}}{list containing the sampled chains of all the parameters after the burn-in period}
#'   \item{\code{exe.time}}{execution time in seconds}
#' }
#' 
#' @export
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
  
  # Check if the data provided are allowed
  if (!is.numeric(y) && !is.vector(y)) stop("'y' must be a numeric vector.")
  if (!is.numeric(X) && !is.matrix(X)) stop("'X' must be a numeric matrix.")
  if (!is.matrix(Z) && !is.list(Z)) stop("'Z' must be a list or a numeric matrix.")
  
  # If Z is a matrix, we cast it into a list
  if (is.matrix(Z))
    Z = list(Z)
  
  # Check if the data dimensions are compatible
  if (length(y) != nrow(X)) stop("'y' and 'X' must have compatible dimensions.")
  for (h in 1:length(Z)) {
    if (length(y) != nrow(Z[[h]])) stop("'y' and 'Z' must have compatible dimensions.")
  }
  
  # Build the completed design matrix
  C = cbind(X, do.call(cbind, Z))
  
  # Fit the Horseshoe linear model
  fit = NULL
  if (default$verbose) {
    fit <- horseshoe::horseshoe(y = y, X = C,
                                method.tau = "halfCauchy",
                                method.sigma = "Jeffreys",
                                burn = default$burn, 
                                nmc = default$niter,
                                thin = default$thin)
  } else {
    invisible(
      capture.output(
        fit <- horseshoe::horseshoe(y = y, X = C,
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
