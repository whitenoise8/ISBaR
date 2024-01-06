# file: fit_cusp_reg.R
# author: Cristian Castiglione
# creation: 05/01/2024
# last change: 05/01/2024
# 
# description: 
#   This file implements a Gibbs sampling algorithm for estimating a Bayesian 
#   additive model predicting the response y given the covariate matrices 
#   X, Z_1, ..., Z_H. Assuming that the columns of each Z_h (h = 1, ..., H) are 
#   sorted in a decresing order of importance, we elicit for the regreesion 
#   coefficients beta_h an increasing shrinkage prior based on the cumulative
#   shrinkage process by Legramanti, Durante, Dunson (2020).
#   
# references:
#   Legramanti, Durante, Dunson (2020)
#   Bayesian cumulative shrinkage for infinite factorizations
#   Biometrika, 107(3): 745-752
#   
# model specification:
# 
#   y = mu + X * beta + e, e ~ N(0, sigma^-2)
#   sigma^-2 ~ IG(A_s,B_s)
#   beta_g ~ N(0, theta_g I_pg), group g = 1,...,G and each group with pg variables 
#   theta_g = (1 - pi_g) * IG(A_th,B_th) + pi_g * delta_th(inf) 
#   pi_g = om_1 + ... + om_g
#   om_h = nu_h( 1-nu_(h-1) )( ... )(1-nu_1)
#   nu_m ~ Beta(1,A_nu)
#   
#   Integrate out theta_g:
#   beta_g = (1 - pi_g)t(2A_th)(0,B_th/A_th I_pg) + pi_g*N(0,th(inf)I_pg)
#   th(inf) -> 0, th(inf) > 0 (continuous shrinkage)
#   if B_th/A_th > th(inf) we have CUSP
#   
#   where A_s, A_th, A_nu, B_s, B_th, th(inf) > 0 are fixed hyperparameters. 
# 

# Take a look to the following articles for an efficient sampling from
# the full-conditional distribution of beta:
# 
# Bhattacharya, Chakraborty, Mallik (2016)
# Fast sampling with Gaussian scale mixture priors in high-dimensional regression
# Biometrika, 103(4): 985-991
# 
# Rue (2001)
# Fast sampling of Gaussian Markov random fields
# Journal of the Royal Statistical Society, Series B, 63(): 325-338
#


#' @title Fit a Bayesian linear model via Gibbs sampling
#' 
#' @description
#' A short description...
#' 
#' @param y response vector
#' @param X fixed effect design matrix, including all the variables that must not be penalized
#' @param Z list of random effect matrices, including all the variable subject to regularization
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
#'   \item{\code{Z}}{liset of random effect matrices}
#'   \item{\code{prior}}{list of prior parameters}
#'   \item{\code{control}}{list of control parameters}
#'   \item{\code{burn}}{list containig the sampled chains of all the parameters during the burn-in iterations}
#'   \item{\code{trace}}{list containing the sampled chains of all the parameters after the burn-in period}
#'   \item{\code{exe.time}}{execution time in seconds}
#' }
#' 
#' @export
cusp.reg.fit = function (y, X, Z, prior = list(), control = list()) {
  
  # This second implementation adaptively prunes the irrelevant basis along the run
  # y: n x 1 vector of response variables
  # X: n x p1 matrix of covariates
  # Z: list of n x ph basis matrices (h = 2, ..., H)
  
  # Get the initial CPU time
  time0 = proc.time()
  
  # Set the prior and control parameters
  prior = init.cusp.prior(prior)
  control = init.cusp.control(control)
  maxiter = control$burn + control$niter
  
  # Log-prior density function
  get.logprior = function (beta, pi, nu, tau, psi, idx, q, prior) {
    out = 0
    out = out + stats::dgamma(psi, shape = prior$a, rate = prior$b, log = TRUE)
    for (h in 2:length(idx)) {
      qh = q[h]
      ih = idx[[h]][1:qh]
      tauh = tau[h]
      pih = pi[ih]
      nuh = nu[ih]
      betah = beta[ih]
      thetah = pih * prior$t0 + (1 - pih) / tauh 
      out = out + stats::dgamma(tauh, shape = prior$at, rate = prior$bt, log = TRUE)
      out = out + sum(stats::dbeta(nuh, shape1 = 1, shape2 = prior$alpha, log = TRUE))
      out = out + sum(stats::dnorm(betah, mean = 0, sd = sqrt(thetah), log = TRUE))
    }
    return (out)
  }
  
  # Log-likelihood function
  get.loglik = function (y, eta, psi) {
    out = sum(stats::dnorm(y, mean = eta, sd = sqrt(1 / psi), log = TRUE))
    return (out)
  }
  
  # Check if the data provided are allowed
  if (!is.numeric(y) && !is.vector(y)) stop("'y' must be a numeric vector.")
  if (!is.numeric(X) && !is.matrix(X)) stop("'X' must be a numeric matrix.")
  if (!is.matrix(Z) && !is.list(Z)) stop("'Z' must be a list or a numeric matrix.")
  
  # If Z is a matrix, we cast it into a list
  if (is.matrix(Z)) Z = list(Z)
  
  # Check if the data dimensions are compatible
  if (length(y) != nrow(X)) stop("'y' and 'X' must have compatible dimensions.")
  for (h in 1:length(Z)) {
    if (length(y) != nrow(Z[[h]])) stop("'y' and 'Z' must have compatible dimensions.")
  }
  
  # Data dimenstions
  n = length(y)
  p = get.block.dim(X, Z)
  q = p
  H = length(p)
  k = sum(p)
  
  # Block index list
  idx = get.block.idx(p)
  keep = get.flat.idx(idx, p)
  drop = c()
  
  # Precompute the sufficient statistics
  C = get.matrix(X, Z)
  CC = crossprod(C, C)
  Cy = crossprod(C, y)
  
  # Allocate the memory for the Markov chains
  burn = list(beta     = matrix(NA, nrow = control$burn, ncol = k),
              theta    = matrix(NA, nrow = control$burn, ncol = k),
              pi       = matrix(NA, nrow = control$burn, ncol = k),
              omega    = matrix(NA, nrow = control$burn, ncol = k),
              nu       = matrix(NA, nrow = control$burn, ncol = k),
              tau      = matrix(NA, nrow = control$burn, ncol = H),
              psi      = rep(NA, length = control$burn),
              logprior = rep(NA, length = control$burn),
              loglik   = rep(NA, length = control$burn),
              logpost  = rep(NA, length = control$burn),
              npar     = rep(NA, length = control$burn))
  
  trace = list(beta     = matrix(NA, nrow = control$niter, ncol = k),
               theta    = matrix(NA, nrow = control$niter, ncol = k),
               pi       = matrix(NA, nrow = control$niter, ncol = k),
               omega    = matrix(NA, nrow = control$niter, ncol = k),
               nu       = matrix(NA, nrow = control$niter, ncol = k),
               tau      = matrix(NA, nrow = control$niter, ncol = H),
               psi      = rep(NA, length = control$niter),
               logprior = rep(NA, length = control$niter),
               loglik   = rep(NA, length = control$niter),
               logpost  = rep(NA, length = control$niter),
               npar     = rep(NA, length = control$niter))
  
  # Initialize the unknown parameters
  init  = init.cusp.param(y, C, idx, prior)
  beta  = init$beta
  theta = init$theta
  pi    = init$pi
  omega = init$omega
  nu    = init$nu
  tau   = init$tau
  psi   = init$psi
  eta   = C %*% beta
  
  # Compute the log prior, likelihood and posterior
  logprior = get.logprior(beta, pi, nu, tau, psi, idx, q, prior)
  loglik = get.loglik(y, eta, psi)
  logpost = logprior + loglik
  
  # Store the initial values
  burn$beta[1,]    = beta
  burn$theta[1,]   = theta
  burn$pi[1,]      = pi
  burn$omega[1,]   = omega
  burn$nu[1,]      = nu
  burn$tau[1,]     = tau
  burn$psi[1]      = psi
  burn$logprior[1] = logprior
  burn$loglik[1]   = loglik
  burn$logpost[1]  = logpost
  
  if (control$verbose) {
    cat(c(rep("-", 53), "\n"), sep = "")
    print.status(0, maxiter, keep, logpost, control$report, control$burn)
  }
  
  # Gibbs sampling loop
  for (iter in 2:maxiter) {
    
    # Print the MCMC status
    if (control$verbose) {
      print.status(iter, maxiter, keep, logpost, control$report, control$burn)
    }
    
    # Update: beta
    D = diag(1 / theta[keep])
    A = psi * CC[keep,keep] + D
    b = psi * drop(Cy)[keep]
    beta[keep] = rmnorm(A, b)
    beta[drop] = 0
    
    # Update: eta
    eta = drop(C[,keep] %*% beta[keep])
    
    # Update: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum((y - eta)^2)
    psi = stats::rgamma(1, shape = ap, rate = bp)
    
    # Update: z, pi, omega, nu
    for (h in 2:H) {
      # Set the block-specific parameters
      ih = idx[[h]]
      qh = q[h]
      ph = length(ih)
      
      betah = beta[ih]
      omegah = omega[ih]
      nuh = nu[ih]
      tauh = tau[h]
      
      # Keep the relevant variables
      keeph = 1:qh
      droph = setdiff(1:ph, keeph)
      
      # Update tau_h
      at = prior$at + 0.5 * qh
      bt = prior$bt + 0.5 * sum(betah[keeph]^2)
      tauh = stats::rgamma(1, shape = at, rate = bt)
      
      # Update zh
      uph = upper.tri(diag(qh), diag = TRUE)
      loh = lower.tri(diag(qh), diag = FALSE)
      nbh = dnorm(betah[keeph], mean = 0, sd = sqrt(prior$t0))
      tbh = dt(betah[keeph] * (prior$at / prior$bt), df = 2 * prior$at) * (prior$at / prior$bt)
      prh = uph * tcrossprod(omegah, nbh) + loh * tcrossprod(omegah, tbh)
      zh = drop(apply(prh, 2, function (x) sample(1:qh, size = 1, prob = x)))
      
      # Update nu_h
      ah = rep(1, length = qh)
      bh = rep(1, length = qh)
      for (jh in 1:qh) {
        ah[jh] = 1 + sum(zh = jh)
        bh[jh] = prior$alpha + sum(zh > jh)
      }
      nuh = stats::rbeta(qh, shape1 = ah, shape2 = bh)
      nuh[qh] = 1
      
      # Update omega_h and pi_h
      omegah = nuh * cumprod(c(1, 1 - nuh[-1]))
      pih = cumsum(omegah)
      
      # Update theta_h
      thetah = rep(prior$t0, length = qh)
      for (jh in 1:qh) {
        if (zh[jh] > jh) thetah[jh] = 1 / tauh
      }
      
      # Store the updated parameters
      tau[h] = tauh
      theta[ih] = thetah
      pi[ih] = pih
      omega[ih] = omega[ih]
      nu[ih] = nu[h]
    }
    
    # Compute the log prior, likelihood and posterior
    logprior = get.logprior(beta, pi, nu, tau, psi, idx, q, prior)
    loglik = get.loglik(y, eta, psi)
    logpost = logprior + loglik
    
    # Store the initial values
    if (iter <= control$burn) {
      t = iter
      burn$beta[t,]    = beta
      burn$theta[t,]   = theta
      burn$pi[t,]      = pi
      burn$omega[t,]   = omega
      burn$nu[t,]      = nu
      burn$tau[t,]     = tau
      burn$psi[t]      = psi
      burn$logprior[t] = logprior
      burn$loglik[t]   = loglik
      burn$logpost[t]  = logpost
      burn$npar[t]     = length(keep)
    } else {
      t = iter - control$burn
      trace$beta[t,]    = beta
      trace$theta[t,]   = theta
      trace$pi[t,]      = pi
      trace$omega[t,]   = omega
      trace$nu[t,]      = nu
      trace$tau[t,]     = tau
      trace$psi[t]      = psi
      trace$logprior[t] = logprior
      trace$loglik[t]   = loglik
      trace$logpost[t]  = logpost
      trace$npar[t]     = length(keep)
    }
    
    # Adaptation
    if (control$adaptation) {
      a0 = 1; a1 = 0.0005
      pr = exp(- a0 - a1 * iter)
      u = runif(1)
      res = y - eta
      if (u < pr) {
        for (h in 2:H) {
          # Compute the partial explained variance of the m-th basis
          qh = q[h]
          ih = idx[[h]][qh]
          etah = C[,ih] * beta[ih]
          resh = res + etah
          expvar = var(resh) / var(res) - 1
          # If the explained variance is small, we drop the last component
          if (qh > control$minnpar && expvar <= control$tol) {
            q[h] = qh - 1
            beta[ih] = 0
            theta[ih] = 0
          }
          # If the explained variance is high, we add a new component
          # and we sample the associated parameters from the prior
          if (qh < p[h] && expvar > control$tol) {
            q[h] = qh + 1
            beta[ih+1] = mean(C[,ih] * res) / mean(C[,ih]^2)
            theta[ih+1] = theta[ih]
          }
          # Update the linear predictor and the residual vectors
          keep = get.flat.idx(idx, q)
          drop = setdiff(1:k, keep)
          # eta = drop(C[,keep] %*% beta[keep])
          # res = y - eta
        }
      }
    }
  }
  
  if (control$verbose) {
    cat(c(rep("-", 53), "\n"), sep = "")
  }
  
  if (control$thin > 1) {
    idt = seq(from = 1, to = control$burn, by = control$thin)
    burn$beta     = burn$beta[idt,]
    burn$theta    = burn$theta[idt,]
    burn$pi       = burn$pi[idt,]
    burn$omega    = burn$omega[idt,]
    burn$nu       = burn$nu[idt,]
    burn$tau      = burn$tau[idt,]
    burn$psi      = burn$psi[idt]
    burn$logprior = burn$logprior[idt]
    burn$loglik   = burn$loglik[idt]
    burn$logpost  = burn$logpost[idt]
    burn$npar     = burn$npar[idt]
    
    idt = seq(from = 1, to = control$niter, by = control$thin)
    trace$beta     = trace$beta[idt,]
    trace$theta    = trace$theta[idt,]
    trace$pi       = trace$pi[idt,]
    trace$omega    = trace$omega[idt,]
    trace$nu       = trace$nu[idt,]
    trace$tau      = trace$tau[idt,]
    trace$psi      = trace$psi[idt]
    trace$logprior = trace$logprior[idt]
    trace$loglik   = trace$loglik[idt]
    trace$npar     = trace$npar[idt]
  }
  
  # Get the final CPU time
  timef = proc.time()
  
  # Return the estimated models
  list(
    y = y, X = X, Z = Z, idx = idx,
    prior = prior, control = control, 
    burn = burn, trace = trace,
    exe.time = (timef - time0)[3])
}




