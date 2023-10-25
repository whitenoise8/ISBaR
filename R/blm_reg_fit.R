# file: blm_reg_fit.R
# author: Cristian Castiglione
# creation: 09/08/2023
# last change: 09/08/2023
# 
# description: 
#   This file implements a Gibbs sampling algorithm for estimating a Bayesian 
#   linear model predicting the response y given the covariate matrix X.
#
# model specification:
# 
#   y = X beta + Z u + e, e ~ N(0, psi^-1)
#   beta ~ N(0, tau^-1 I)
#   tau ~ Gamma(v/2, v/2)
#   psi ~ Gamma(a, b)
#   
#   where a, b, v > 0 are fixed hyperparameters.
# 

#' @title Fit a Bayesian linear model via Gibbs sampling
#' 
#' @description
#' A short description...
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
#' @return Returns an object of class "blm", which is a list containing the following elements:
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
blm.reg.fit = function (y, X, Z, prior = list(), control = list()) {
  
  time0 = proc.time()
  
  # Set the prior and control parameters
  prior = init.blm.prior(prior)
  control = init.blm.control(control)
  totiter = control$burn + control$niter
  
  # Log-prior density function
  get.logprior = function (beta, tau, psi, idx, prior) {
    out = 0
    out = out + dgamma(psi, shape = prior$a, rate = prior$b, log = TRUE)
    for (h in 2:length(idx)) {
      out = out + dgamma(tau[h], shape = 0.5 * prior$v, rate = 0.5 * prior$v, log = TRUE)
      out = out + sum(dnorm(beta[idx[[h]]], mean = 0, sd = 1 / sqrt(tau[h]), log = TRUE))
    }
    return (out)
  }
  
  # Log-likelihood function
  get.loglik = function (y, eta, psi) {
    sum(dnorm(y, mean = eta, sd = 1 / sqrt(psi), log = TRUE))
  }
  
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
  
  # Data dimenstions
  n = nrow(X)
  p = get.block.dim(X, Z)
  H = length(p)
  K = sum(p)
  
  idx = get.block.idx(p)
  
  # Precompute the sufficient statistics
  C = get.matrix(X, Z)
  A = crossprod(C, C)
  b = crossprod(C, y)
  
  # Set the duplication matrix
  L = matrix(0, nrow = K, ncol = H)
  for (h in 2:H) L[idx[[h]],h] = 1

  # Allocate the memory for the Markov chains
  burn = list(beta     = matrix(NA, nrow = control$burn, ncol = K),
              tau      = matrix(NA, nrow = control$burn, ncol = H),
              psi      = rep(NA, length = control$burn),
              logprior = rep(NA, length = control$burn),
              loglik   = rep(NA, length = control$burn),
              logpost  = rep(NA, length = control$burn))
  
  trace = list(beta     = matrix(NA, nrow = control$niter, ncol = K),
               tau      = matrix(NA, nrow = control$niter, ncol = H),
               psi      = rep(NA, length = control$niter),
               logprior = rep(NA, length = control$niter),
               loglik   = rep(NA, length = control$niter),
               logpost  = rep(NA, length = control$niter))
  
  # Initialize the unknown parameters
  D = diag(drop(L %*% (.1 / p)))
  beta = drop(solve(A + D, b))
  eta = drop(C %*% beta)
  
  ap = prior$a + 0.5 * n
  bp = prior$b + 0.5 * sum((y - eta)^2)
  psi = rgamma(1, shape = ap, rate = bp)
  
  tau = rep(0, length = H)
  for (h in 2:H) {
    at = 0.5 * (prior$v + p[h])
    bt = 0.5 * (prior$v + sum(beta[idx[[h]]]^2))
    tau[h] = rgamma(1, shape = at, rate = bt)
  }
  
  # Compute the log prior, likelihood and posterior
  logprior = get.logprior(beta, tau, psi, idx, prior)
  loglik = get.loglik(y, eta, psi)
  logpost = logprior + loglik
  
  # Store the initial values
  burn$beta[1,]    = beta
  burn$tau[1,]     = tau
  burn$psi[1]      = psi
  burn$logprior[1] = logprior
  burn$loglik[1]   = loglik
  burn$logpost[1]  = logpost
  
  if (control$verbose) {
    cat(c(rep("-", 50), "\n"), sep = "")
  }
  
  for (iter in 2:totiter) {
    
    if (control$verbose) {
      if (iter %% control$report == 0) {
        cat(gettextf(" iter: %5d /%5d \t", iter, totiter),
            gettextf(" logp: %.4f \n", logpost))
      }
      
      if (iter == control$burn) {
        cat(c(rep("-", 50), "\n"), sep = "")
      }
    }
    
    # Update: beta
    D = diag(drop(L %*% tau))
    beta = rmnorm(psi * A + D, psi * b)
    eta = drop(C %*% beta)
    
    # Update: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum((y - eta)^2)
    psi = rgamma(1, shape = ap, rate = bp)
    
    # Update: tau
    for (h in 2:H) {
      at = 0.5 * (prior$v + p[h])
      bt = 0.5 * (prior$v + sum(beta[idx[[h]]]^2))
      tau[h] = rgamma(1, shape = at, rate = bt)
    }
    
    # Compute the log prior, likelihood and posterior
    logprior = get.logprior(beta, tau, psi, idx, prior)
    loglik = get.loglik(y, eta, psi)
    logpost  = logprior + loglik
    
    # Store the initial values
    if (iter <= control$burn) {
      t = iter
      burn$beta[t,]    = beta
      burn$tau[t,]     = tau
      burn$psi[t]      = psi
      burn$logprior[t] = logprior
      burn$loglik[t]   = loglik
      burn$logpost[t]  = logpost
    } else {
      t = iter - control$burn
      trace$beta[t,]    = beta
      trace$tau[t,]     = tau
      trace$psi[t]      = psi
      trace$logprior[t] = logprior
      trace$loglik[t]   = loglik
      trace$logpost[t]  = logpost
    }
    
  }
  
  if (control$verbose) {
    cat(c(rep("-", 50), "\n"), sep = "")
  }
  
  if (control$thin > 1) {
    idt = seq(from = 1, to = control$burn, by = control$thin)
    burn$beta     = burn$beta[idt,]
    burn$delta    = burn$delta[idt,]
    burn$phi      = burn$phi[idt,]
    burn$tau      = burn$tau[idt,]
    burn$psi      = burn$psi[idt]
    burn$logprior = burn$logprior[idt]
    burn$loglik   = burn$loglik[idt]
    burn$logpost  = burn$logpost[idt]
    
    idt = seq(from = 1, to = control$niter, by = control$thin)
    trace$beta     = trace$beta[idt,]
    trace$delta    = trace$delta[idt,]
    trace$phi      = trace$phi[idt,]
    trace$tau      = trace$tau[idt,]
    trace$psi      = trace$psi[idt]
    trace$logprior = trace$logprior[idt]
    trace$loglik   = trace$loglik[idt]
    trace$logpost  = trace$logpost[idt]
  }
  
  timef = proc.time()
  
  list(y = y, X = X, Z = Z, idx = idx,
       prior = prior, control = control, 
       burn = burn, trace = trace,
       exe.time = (timef - time0)[3])
}

# Function for initializing the prior parameters of the BLM
init.blm.prior = function (prior) {
  default = list(
    a = 2.01,
    b = 1.01,
    v = 0.1)
  
  if (class(prior) != "list") {
    # Print a warning if `control` is not a list
    warning("The `prior` parameter must be a list \n",
            "Default parameters will be used for the execution.")
  } else {
    if (length(prior) != 0) {
      # Get the default and prior attributes
      nms.default = names(default)
      nms.prior = names(prior)
      nms.undef = nms.prior[!nms.prior %in% nms.default]
      
      # Set the custom hyperparameters
      default[names(prior)] = prior
      
      # Print a warning if some of the prior parameters are not allowed 
      if (length(nms.undef)) {
        warning("Unknown names in prior parameters: ", 
                paste(nms.undef, collapse = ", "))
      }
    }
  }
  
  # Return the control parameters
  return (default)
}

# Function for initializing the control parameters of the algorithms
init.blm.control = function (control) {
  
  # Default control parameters
  default = list(
    niter = 5000,
    burn = 2500,
    thin = 1,
    verbose = TRUE,
    report = 500)
  
  if (class(control) != "list") {
    # Print a warning if `control` is not a list
    warning("The `control` parameter must be a list \n",
            "Default parameters will be used for the execution.")
  } else {
    if (length(control) != 0) {
      # Get the default and control attributes
      nms.default = names(default)
      nms.control = names(control)
      nms.undef = nms.control[!nms.control %in% nms.default]
      
      # Set the custom parameters
      default[names(control)] = control
      
      # Print a warning if some of the control parameters are not allowed 
      if (length(nms.undef)) {
        warning("Unknown names in control: ", 
                paste(nms.undef, collapse = ", "))
      }
    }
  }
  
  # Return the control parameters
  return (default)
}
