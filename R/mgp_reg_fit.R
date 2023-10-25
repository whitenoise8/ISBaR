#' mgp_reg_fit.R
#' author: Cristian Castiglione
#' creation: 04/08/2023
#' last change: 16/10/2023
#' 
#' description: 
#'   This file implements a Gibbs sampling algorithm for estimating a Bayesian 
#'   linear model predicting the response y given the covariate matrix X.
#'   Assuming that the columns of X are sorted in a decresing order of 
#'   importance, we assume for the regreesion coefficients beta an increasing 
#'   shrinkage prior based on the multiplicative gamma process by 
#'   Bhattacharya & Dunson (2011). See also Durante (2017).
#'   
#' references:
#'   Bhattacharya & Dunson (2011)
#'   Sparse Bayesian infinite factor models
#'   Biometrika, 98(2): 291-306
#'   
#'   Durante (2017). 
#'   A note on the multiplicative gamma process
#'   Statistics and Probability Letters, 122: 198-204
#'   
#' model specification:
#' 
#'   y = X beta + Z u + e, e ~ N(0, psi^-1)
#'   beta ~ N(0, eta^-1 I), eta -> inf
#'   u_j ~ N(0, tau^-1 phi_j^-1), (j = 2, ..., p)
#'   phi_j = delta_1 x ... x delta_j, (j = 2, ..., p)
#'   delta_1 ~ Gamma(a1, 1)
#'   delta_j ~ Gamma(a2, 1), (j = 2, ..., p)
#'   tau ~ Gamma(v/2, v/2)
#'   psi ~ Gamma(a, b)
#'   
#'   where a1, a2, a, b, v > 0 are fixed hyperparameters.
#'   To induce an increasing shrinkage on beta, we need to impose a2 > a1.
#'   Moreover, to guarantee for 1 / phi_j (j = 1, ..., n) to have finite
#'   first and second moment, we need a1 > 2 and a3 > 3. 
#' 

#' Take a look to the following articles for an efficient sampling from
#' the full-conditional of beta:
#' Bhattacharya, Chakraborty, Mallik (2016)
#' Fast sampling with Gaussian scale mixture priors in high-dimensional regression
#' Biometrika, 103(4): 985-991
#' 
#' Rue (2001)
#' Fast sampling of Gaussian Markov random fields
#' Journal of the Royal Statistical Society, Series B, 63(): 325-338
#'

mgp.reg.fit = function (y, X, Z, prior = list(), control = list()) {
  
  # Get the initial CPU time
  time0 = proc.time()
  
  # Set the prior and control parameters
  prior = init.mgp.prior(prior)
  control = init.mgp.control(control)
  totiter = control$burn + control$niter
  
  # Log-prior density function
  get.logprior = function (beta, delta, tau, psi, prior) {
    out = 0
    out = out + sum(stats::dnorm(beta, mean = 0, sd = 1 / sqrt(tau * cumprod(delta)), log = TRUE))
    out = out + sum(stats::dgamma(delta[-1], shape = prior$a2, rate = 1, log = TRUE))
    out = out + stats::dgamma(delta[ 1], shape = prior$a1, rate = 1, log = TRUE)
    out = out + stats::dgamma(tau, shape = prior$v / 2, rate = prior$v / 2, log = TRUE)
    out = out + stats::dgamma(psi, shape = prior$a, rate = prior$b, log = TRUE)
    return (out)
  }
  
  # Log-likelihood function
  get.loglik = function (y, eta, psi) {
    out = sum(stats::dnorm(y, mean = eta, sd = sqrt(1/psi), log = TRUE))
    return (out)
  }
  
  # Data dimenstions
  n = nrow(X)
  p = ncol(X)
  q = ncol(Z)
  
  idx = 1:p
  idz = (p+1):(p+q)
  
  # Precompute the sufficient statistics
  XX = crossprod(cbind(X, Z))
  Xy = crossprod(cbind(X, Z), y)
  
  # Allocate the memory for the Markov chains
  burn = list(beta     = matrix(NA, nrow = control$burn, ncol = p+q),
              delta    = matrix(NA, nrow = control$burn, ncol = q),
              phi      = matrix(NA, nrow = control$burn, ncol = q),
              tau      = rep(NA, length = control$burn),
              psi      = rep(NA, length = control$burn),
              logprior = rep(NA, length = control$burn),
              loglik   = rep(NA, length = control$burn),
              logpost  = rep(NA, length = control$burn))
  
  trace = list(beta     = matrix(NA, nrow = control$niter, ncol = p+q),
               delta    = matrix(NA, nrow = control$niter, ncol = q),
               phi      = matrix(NA, nrow = control$niter, ncol = q),
               tau      = rep(NA, length = control$niter),
               psi      = rep(NA, length = control$niter),
               logprior = rep(NA, length = control$niter),
               loglik   = rep(NA, length = control$niter),
               logpost  = rep(NA, length = control$niter))
  
  # Initialize the unknown parameters
  init  = init.mgp.param(y, X, Z, prior)
  mu    = init$mu
  beta  = init$beta
  delta = init$delta
  tau   = init$tau
  psi   = init$psi
  phi   = cumprod(init$delta)
  eta   = X %*% beta[idx] + Z %*% beta[idz]

  # Compute the log prior, likelihood and posterior
  logprior = get.logprior(beta[idz], delta, tau, psi, prior)
  loglik   = get.loglik(y, eta, psi)
  logpost  = logprior + loglik
  
  # Store the initial values
  burn$beta[1,]    = beta
  burn$delta[1,]   = delta
  burn$phi[1,]     = phi
  burn$tau[1]      = tau
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
    D = diag(c(rep(0, p), tau * phi))
    beta = rmnorm(psi * XX + D, psi * Xy)
    eta = drop(X %*% beta[idx] + Z %*% beta[idz])
    
    # Update: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum((y - eta)^2)
    psi = stats::rgamma(1, shape = ap, rate = bp)
    
    # Update: tau
    at = 0.5 * (prior$v + q)
    bt = 0.5 * (prior$v + sum(phi * beta[idz]^2))
    tau = stats::rgamma(1, shape = at, rate = bt)
    
    # Update: delta_1
    ad = prior$a1 + 0.5 * q
    bd = 1 + 0.5 * tau * sum(phi * beta[idz]^2) / delta[1]
    delta[1] = stats::rgamma(1, shape = ad, rate = bd)
    phi = cumprod(delta)
    
    # Update: delta_h
    for (h in 2:q) {
      ad = prior$a2 + 0.5 * (q - h + 1)
      bd = 1 + 0.5 * tau * sum(phi[h:q] * beta[idz[h:q]]^2) / delta[h]
      delta[h] = stats::rgamma(1, shape = ad, rate = bd)
      phi = cumprod(delta)
    }
    
    # Compute the log prior, likelihood and posterior
    logprior = get.logprior(beta[idz], delta, tau, psi, prior)
    loglik   = get.loglik(y, eta, psi)
    logpost  = logprior + loglik
    
    # Store the initial values
    if (iter <= control$burn) {
      t = iter
      burn$beta[t,]    = beta
      burn$delta[t,]   = delta
      burn$phi[t,]     = phi
      burn$tau[t]      = tau
      burn$psi[t]      = psi
      burn$logprior[t] = logprior
      burn$loglik[t]   = loglik
      burn$logpost[t]  = logpost
    } else {
      t = iter - control$burn
      trace$beta[t,]    = beta
      trace$delta[t,]   = delta
      trace$phi[t,]     = phi
      trace$tau[t]      = tau
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
    idx = seq(from = 1, to = control$burn, by = control$thin)
    burn$beta     = burn$beta[idx,]
    burn$delta    = burn$delta[idx,]
    burn$phi      = burn$phi[idx,]
    burn$tau      = burn$tau[idx]
    burn$psi      = burn$psi[idx]
    burn$logprior = burn$logprior[idx]
    burn$loglik   = burn$loglik[idx]
    burn$logpost  = burn$logpost[idx]
    
    idx = seq(from = 1, to = control$niter, by = control$thin)
    trace$beta     = trace$beta[idx,]
    trace$delta    = trace$delta[idx,]
    trace$phi      = trace$phi[idx,]
    trace$tau      = trace$tau[idx]
    trace$psi      = trace$psi[idx]
    trace$logprior = trace$logprior[idx]
    trace$loglik   = trace$loglik[idx]
    trace$logpost  = trace$logpost[idx]
  }
  
  timef = proc.time()
  
  # Return the estimated models
  list(y = y, X = X, Z = Z, 
       prior = prior, control = control, 
       burn = burn, trace = trace,
       exe.time = (timef - time0)[3])
}

mgp.reg.fit2 = function (y, X, Z, prior = list(), control = list()) {
  
  # This second implementation adaptively prunes the irrelevant basis along the run
  
  time0 = proc.time()
  
  # Set the prior and control parameters
  prior = init.mgp.prior(prior)
  control = init.mgp.control(control)
  maxiter = control$burn + control$niter
  
  # Log-prior density function
  get.logprior = function (beta, delta, tau, psi, prior) {
    out = 0
    out = out + sum(stats::dnorm(beta, mean = 0, sd = 1 / sqrt(tau * cumprod(delta)), log = TRUE))
    out = out + sum(stats::dgamma(delta[-1], shape = prior$a2, rate = 1, log = TRUE))
    out = out + stats::dgamma(delta[ 1], shape = prior$a1, rate = 1, log = TRUE)
    out = out + stats::dgamma(tau, shape = prior$v / 2, rate = prior$v / 2, log = TRUE)
    out = out + stats::dgamma(psi, shape = prior$a, rate = prior$b, log = TRUE)
    return (out)
  }
  
  # Log-likelihood function
  get.loglik = function (y, eta, psi) {
    out = sum(stats::dnorm(y, mean = eta, sd = sqrt(1/psi), log = TRUE))
    return (out)
  }
  
  # Data dimenstions
  n = nrow(X)
  p = ncol(X)
  q = ncol(Z)
  m = q
  
  idx = 1:p
  idz = (p+1):(p+q)
  
  # Precompute the sufficient statistics
  XX = crossprod(cbind(X, Z))
  Xy = crossprod(cbind(X, Z), y)
  
  # Allocate the memory for the Markov chains
  burn = list(beta     = matrix(.0, nrow = control$burn, ncol = p+q),
              delta    = matrix(NA, nrow = control$burn, ncol = q),
              phi      = matrix(NA, nrow = control$burn, ncol = q),
              tau      = rep(NA, length = control$burn),
              psi      = rep(NA, length = control$burn),
              logprior = rep(NA, length = control$burn),
              loglik   = rep(NA, length = control$burn),
              logpost  = rep(NA, length = control$burn),
              npar     = rep(NA, length = control$burn))
  
  trace = list(beta     = matrix(.0, nrow = control$niter, ncol = p+q),
               delta    = matrix(NA, nrow = control$niter, ncol = q),
               phi      = matrix(NA, nrow = control$niter, ncol = q),
               tau      = rep(NA, length = control$niter),
               psi      = rep(NA, length = control$niter),
               logprior = rep(NA, length = control$niter),
               loglik   = rep(NA, length = control$niter),
               logpost  = rep(NA, length = control$niter),
               npar     = rep(NA, length = control$niter))
  
  # Initialize the unknown parameters
  init  = init.mgp.param(y, X, Z, prior)
  mu    = init$mu
  beta  = init$beta
  delta = init$delta
  tau   = init$tau
  psi   = init$psi
  phi   = cumprod(init$delta)
  eta   = X %*% beta[idx] + Z %*% beta[idz]
  
  # Compute the log prior, likelihood and posterior
  logprior = get.logprior(beta[idz], delta, tau, psi, prior)
  loglik   = get.loglik(y, eta, psi)
  logpost  = logprior + loglik
  
  # Store the initial values
  burn$beta[1,]    = beta
  burn$delta[1,]   = delta
  burn$phi[1,]     = phi
  burn$tau[1]      = tau
  burn$psi[1]      = psi
  burn$logprior[1] = logprior
  burn$loglik[1]   = loglik
  burn$logpost[1]  = logpost
  
  if (control$verbose) {
    cat(c(rep("-", 53), "\n"), sep = "")
    cat(gettextf(" iter: %5d /%5d  ", 0, maxiter),
        gettextf(" ncomp: %2d  ", length(beta)),
        gettextf(" logp: %.4f \n", logpost))
  }
  
  for (iter in 2:maxiter) {
    
    if (control$verbose) {
      if (iter %% control$report == 0) {
        cat(gettextf(" iter: %5d /%5d  ", iter, maxiter),
            gettextf(" ncomp: %2d  ", length(beta)),
            gettextf(" logp: %.4f \n", logpost))
      }
      
      if (iter == control$burn) {
        cat(c(rep("-", 53), "\n"), sep = "")
      }
    }
    
    # Update: beta
    D = diag(c(rep(0, p), tau * phi))
    beta = rmnorm(psi * XX + D, psi * Xy)
    eta = drop(X %*% beta[idx] + Z[,1:m] %*% beta[idz[1:m]])
    res = y - eta
    
    # Update: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum(res^2)
    psi = stats::rgamma(1, shape = ap, rate = bp)
    
    # Update: tau
    at = 0.5 * (prior$v + q)
    bt = 0.5 * (prior$v + sum(phi * beta[idz[1:m]]^2))
    tau = stats::rgamma(1, shape = at, rate = bt)
    
    # Update: delta_1
    ad = prior$a1 + 0.5 * q
    bd = 1 + 0.5 * tau * sum(phi * beta[idz[1:m]]^2) / delta[1]
    delta[1] = stats::rgamma(1, shape = ad, rate = bd)
    phi = cumprod(delta)
    
    # Update: delta_h
    for (h in 2:m) {
      ad = prior$a2 + 0.5 * (m - h + 1)
      bd = 1 + 0.5 * tau * sum(phi[h:m] * beta[idz[h:m]]^2) / delta[h]
      delta[h] = stats::rgamma(1, shape = ad, rate = bd)
      phi = cumprod(delta)
    }
    
    # Compute the log prior, likelihood and posterior
    logprior = get.logprior(beta[idz[1:m]], delta, tau, psi, prior)
    loglik   = get.loglik(y, eta, psi)
    logpost  = logprior + loglik
    
    # Store the initial values
    if (iter <= control$burn) {
      t = iter
      burn$beta[t,1:(p+m)] = beta
      burn$delta[t,1:m]    = delta
      burn$phi[t,1:m]      = phi
      burn$tau[t]          = tau
      burn$psi[t]          = psi
      burn$logprior[t]     = logprior
      burn$loglik[t]       = loglik
      burn$logpost[t]      = logpost
      burn$npar[t]         = length(beta)
    } else {
      t = iter - control$burn
      trace$beta[t,1:(p+m)] = beta
      trace$delta[t,1:m]    = delta
      trace$phi[t,1:m]      = phi
      trace$tau[t]          = tau
      trace$psi[t]          = psi
      trace$logprior[t]     = logprior
      trace$loglik[t]       = loglik
      trace$logpost[t]      = logpost
      trace$npar[t]         = length(beta)
    }
    
    # Adaptation
    if (control$adaptation) {
      tol = 1e-04
      a0 = 1; a1 = 0.0005
      pr = exp(- a0 - a1 * iter)
      u = runif(1)
      if (u < pr) {
        # Compute the partial explained variance of the m-th basis
        etam = Z[,m] * beta[idz[m]]
        resm = res + etam
        expvar = var(resm) / var(res) - 1
        # If the explained variance is small, we drop the last component
        if (m > control$minnpar && expvar <= control$tol) {
          m = m - 1
          idz = idz[1:m]
          beta = beta[1:(p+m)]
          delta = delta[1:m]
          phi = phi[1:m]
        }
        # If the explained variance is high, we add a new component
        # and we sample the associated parameters from the prior
        if (m < q && expvar > control$tol) {
          m = m + 1
          idz = c(idz, p+m)
          beta = c(beta, mean(Z[,m] * res) / mean(Z[,m]^2))
          delta = c(delta, delta[m-1])
          phi = c(phi, phi[m-1] + delta[m-1])
        }
        # Finally, we update the sufficient statistics
        XX = crossprod(cbind(X, Z[,1:m]))
        Xy = crossprod(cbind(X, Z[,1:m]), y)
      }
    }
  }
  
  if (control$verbose) {
    cat(c(rep("-", 53), "\n"), sep = "")
  }
  
  if (control$thin > 1) {
    idx = seq(from = 1, to = control$burn, by = control$thin)
    burn$beta     = burn$beta[idx,]
    burn$delta    = burn$delta[idx,]
    burn$phi      = burn$phi[idx,]
    burn$tau      = burn$tau[idx]
    burn$psi      = burn$psi[idx]
    burn$logprior = burn$logprior[idx]
    burn$loglik   = burn$loglik[idx]
    burn$logpost  = burn$logpost[idx]
    burn$npar     = burn$npar[idx]
    
    idx = seq(from = 1, to = control$niter, by = control$thin)
    trace$beta     = trace$beta[idx,]
    trace$delta    = trace$delta[idx,]
    trace$phi      = trace$phi[idx,]
    trace$tau      = trace$tau[idx]
    trace$psi      = trace$psi[idx]
    trace$logprior = trace$logprior[idx]
    trace$loglik   = trace$loglik[idx]
    trace$logpost  = trace$logpost[idx]
    trace$npar     = trace$npar[idx]
  }
  
  timef = proc.time()
  
  # Return the estimated models
  list(y = y, X = X, Z = Z, 
       prior = prior, control = control, 
       burn = burn, trace = trace,
       exe.time = (timef - time0)[3])
}

