#' mgp_gam_fit.R
#' author: Cristian Castiglione
#' creation: 15/10/2023
#' last change: 16/10/2023
#' 
#' description: 
#'   This file implements a Gibbs sampling algorithm for estimating a Bayesian 
#'   additive model predicting the response y given the covariate matrices 
#'   X, Z_1, ..., Z_H. Assuming that the columns of each Z_h (h = 1, ..., H) are 
#'   sorted in a decresing order of importance, we elicit for the regreesion 
#'   coefficients beta_h an increasing shrinkage prior based on the multiplicative 
#'   gamma process by Bhattacharya & Dunson (2011). See also Durante (2017).
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
#'   y = X beta_0 + Z1 beta_1 + ... + Z_H beta_H + e, 
#'   e ~ N(0, psi^-1), psi ~ Gamma(a, b),
#'   beta_0 ~ N(0, kappa^-1 I) (kappa -> inf),
#'   beta_jh ~ N(0, tau_h^-1 phi_jh^-1),
#'   phi_jh = delta_1h x ... x delta_jh,
#'   delta_1h ~ Gamma(a1, 1), delta_jh ~ Gamma(a2, 1) (j > 1)
#'   tau_h ~ Gamma(v/2, v/2),
#'   
#'   for j = 1, ..., ph and h = 1, ..., H, where a1, a2, a, b, v > 0 are 
#'   fixed hyperparameters.
#'   To induce an increasing shrinkage on beta, we need to impose a2 > a1.
#'   Moreover, to guarantee for 1 / phi_jh (j = 1, ..., ph) to have finite
#'   first and second moment, we need a1 > 2 and a3 > 3. 
#' 

#' Take a look to the following articles for an efficient sampling from
#' the full-conditional distribution of beta:
#' 
#' Bhattacharya, Chakraborty, Mallik (2016)
#' Fast sampling with Gaussian scale mixture priors in high-dimensional regression
#' Biometrika, 103(4): 985-991
#' 
#' Rue (2001)
#' Fast sampling of Gaussian Markov random fields
#' Journal of the Royal Statistical Society, Series B, 63(): 325-338
#'
mgp.gam.fit = function (y, X, Z, prior = list(), control = list()) {
  
  # y: n x 1 vector of response variables
  # X: n x p1 matrix of covariates
  # Z: list of n x ph basis matrices (h = 2, ..., H)
  
  # Get the initial CPU time
  time0 = proc.time()
  
  # Set the prior and control parameters
  prior = init.mgp.prior(prior)
  control = init.mgp.control(control)
  maxiter = control$burn + control$niter
  
  # Initialize all model parameters
  init.param = function (y, C, idx, prior) {
    n = length(y)
    p = unlist(lapply(idx, length))
    H = length(p)
    k = sum(p)
    
    # Init: beta
    lambda = rep(0, length = k)
    for (h in 2:H) {
      lambda[idx[[h]]] = rep(0.1 / p[h], p[h])
    }
    A = crossprod(C) + diag(lambda)
    b = crossprod(C, y)
    beta = rmnorm(A, b)
    
    # Init: eta
    eta = drop(C %*% beta)
    
    # Init: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum((y - eta)^2)
    psi = stats::rgamma(1, shape = ap, rate = ap)
    
    # Init: tau
    tau = rep(0, length = H)
    for (h in 2:H) {
      ih = idx[[h]]
      at = 0.5 * (prior$v + p[h])
      bt = 0.5 * (prior$v + sum(beta[ih]^2))
      tau[h] = stats::rgamma(1, shape = at, rate = at)
    }
    
    # Init: delta and phi
    delta = rep(0, length = k)
    phi = rep(0, length = k)
    for (h in 2:H) {
      ph = p[h]
      ih = idx[[h]]
      ad = prior$a1 + 0.5 * ph
      bd = 1 + 0.5 * tau[h] * beta[ih]^2
      deltah = rep(NA, length = ph)
      deltah[1] = stats::rgamma(1, shape = ad, rate = ad)
      for (jh in 2:ph) {
        ad = prior$a2 + 0.5 * (ph - jh + 1)
        bd = 1 + 0.5 * tau[h] * beta[ih]^2
        deltah[jh] = stats::rgamma(1, shape = ad, rate = ad)
      }
      delta[ih] = deltah
      phi[ih] = cumprod(deltah)
    }
    
    # output
    list(beta = beta, delta = delta, phi = phi, tau = tau, psi = psi)
  }
  
  # Log-prior density function
  get.logprior = function (beta, delta, tau, psi, idx, prior) {
    out = 0
    out = out + stats::dgamma(psi, shape = prior$a, rate = prior$b, log = TRUE)
    for (h in 2:length(idx)) {
      ih = idx[[h]]
      tauh = tau[h]
      deltah = delta[ih]
      betah = beta[ih]
      lambdah = tauh * cumprod(deltah)
      out = out + stats::dgamma(tauh, shape = prior$v / 2, rate = prior$v / 2, log = TRUE)
      out = out + stats::dgamma(deltah[1], shape = prior$a1, rate = 1, log = TRUE)
      out = out + sum(stats::dgamma(deltah[-1], shape = prior$a2, rate = 1, log = TRUE))
      out = out + sum(stats::dnorm(betah, mean = 0, sd = 1 / sqrt(lambdah), log = TRUE))
    }
    return (out)
  }
  
  # Log-likelihood function
  get.loglik = function (y, eta, psi) {
    out = sum(stats::dnorm(y, mean = eta, sd = sqrt(1 / psi), log = TRUE))
    return (out)
  }
  
  # Print simulation status
  print.status = function (iter, maxiter, logpost, report, burn) {
    if (iter %% report == 0) {
      cat(gettextf(" iter: %5d /%5d \t", iter, maxiter),
          gettextf(" logp: %.4f \n", logpost))
    }
    if (iter == burn) {
      cat(c(rep("-", 50), "\n"), sep = "")
    }
  }
  
  # Data dimenstions
  n = nrow(X)
  p = get.block.dim(X, Z)
  H = length(p)
  k = sum(p)
  
  # Block index list
  idx = get.block.idx(p)
  
  # Precompute the sufficient statistics
  C = get.matrix(X, Z)
  CC = crossprod(C, C)
  Cy = crossprod(C, y)
  
  # Allocate the memory for the Markov chains
  burn = list(beta     = matrix(NA, nrow = control$burn, ncol = k),
              delta    = matrix(NA, nrow = control$burn, ncol = k),
              phi      = matrix(NA, nrow = control$burn, ncol = k),
              tau      = matrix(NA, nrow = control$burn, ncol = H),
              psi      = rep(NA, length = control$burn),
              logprior = rep(NA, length = control$burn),
              loglik   = rep(NA, length = control$burn),
              logpost  = rep(NA, length = control$burn))
  
  trace = list(beta     = matrix(NA, nrow = control$niter, ncol = k),
               delta    = matrix(NA, nrow = control$niter, ncol = k),
               phi      = matrix(NA, nrow = control$niter, ncol = k),
               tau      = matrix(NA, nrow = control$niter, ncol = H),
               psi      = rep(NA, length = control$niter),
               logprior = rep(NA, length = control$niter),
               loglik   = rep(NA, length = control$niter),
               logpost  = rep(NA, length = control$niter))
  
  # Initialize the unknown parameters
  init  = init.param(y, C, idx, prior)
  beta  = init$beta
  delta = init$delta
  tau   = init$tau
  psi   = init$psi
  phi   = init$phi
  eta   = drop(C %*% beta)
  
  # Compute the log prior, likelihood and posterior
  logprior = get.logprior(beta, delta, tau, psi, idx, prior)
  loglik = get.loglik(y, eta, psi)
  logpost = logprior + loglik
  
  # Store the initial values
  burn$beta[1,]    = beta
  burn$delta[1,]   = delta
  burn$phi[1,]     = phi
  burn$tau[1,]     = tau
  burn$psi[1]      = psi
  burn$logprior[1] = logprior
  burn$loglik[1]   = loglik
  burn$logpost[1]  = logpost
  
  if (control$verbose) {
    cat(c(rep("-", 50), "\n"), sep = "")
  }
  
  # Gibbs sampling loop
  for (iter in 2:maxiter) {
    
    # Print the MCMC status
    if (control$verbose) {
      print.status(iter, maxiter, logpost, control$report, control$burn)
    }
    
    # Update: beta
    D = diag(get.penalty(tau, phi, idx))
    beta = rmnorm(psi * CC + D, psi * Cy)
    eta = drop(C %*% beta)
    
    # Update: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum((y - eta)^2)
    psi = stats::rgamma(1, shape = ap, rate = bp)
    
    # Update: tau, delta & phi
    for (h in 2:H) {
      # Block-specific parameters
      ph = p[h]
      ih = idx[[h]]
      betah = beta[ih]
      deltah = delta[ih]
      phih = phi[ih]
      tauh = tau[h]
      
      # Update tau_h
      at = 0.5 * (prior$v + ph)
      bt = 0.5 * (prior$v + sum(phih * betah^2))
      tauh = stats::rgamma(1, shape = at, rate = bt)
      
      # Update delta_1h
      ad = prior$a1 + 0.5 * ph
      bd = 1 + 0.5 * tauh * sum(phih * betah^2) / deltah[1]
      deltah[1] = stats::rgamma(1, shape = ad, rate = bd)
      phih = cumprod(deltah)
      
      # Update delta_jh
      for (j in 2:ph) {
        ad = prior$a2 + 0.5 * (ph - j + 1)
        bd = 1 + 0.5 * tauh * sum(phih[j:ph] * betah[j:ph]^2) / deltah[j]
        deltah[j] = stats::rgamma(1, shape = ad, rate = bd)
        phih = cumprod(deltah)
      }
      
      # Allocation
      tau[h] = tauh
      delta[ih] = deltah
      phi[ih] = phih
    }
    
    # Compute the log prior, likelihood and posterior
    logprior = get.logprior(beta, delta, tau, psi, idx, prior)
    loglik = get.loglik(y, eta, psi)
    logpost = logprior + loglik
    
    # Store the initial values
    if (iter <= control$burn) {
      t = iter
      burn$beta[t,]    = beta
      burn$delta[t,]   = delta
      burn$phi[t,]     = phi
      burn$tau[t,]     = tau
      burn$psi[t]      = psi
      burn$logprior[t] = logprior
      burn$loglik[t]   = loglik
      burn$logpost[t]  = logpost
    } else {
      t = iter - control$burn
      trace$beta[t,]    = beta
      trace$delta[t,]   = delta
      trace$phi[t,]     = phi
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
    idx = seq(from = 1, to = control$burn, by = control$thin)
    burn$beta     = burn$beta[idx,]
    burn$delta    = burn$delta[idx,]
    burn$phi      = burn$phi[idx,]
    burn$tau      = burn$tau[idx,]
    burn$psi      = burn$psi[idx]
    burn$logprior = burn$logprior[idx]
    burn$loglik   = burn$loglik[idx]
    burn$logpost  = burn$logpost[idx]
    
    idx = seq(from = 1, to = control$niter, by = control$thin)
    trace$beta     = trace$beta[idx,]
    trace$delta    = trace$delta[idx,]
    trace$phi      = trace$phi[idx,]
    trace$tau      = trace$tau[idx,]
    trace$psi      = trace$psi[idx]
    trace$logprior = trace$logprior[idx]
    trace$loglik   = trace$loglik[idx]
    trace$logpost  = trace$logpost[idx]
  }
  
  # Get the final CPU time
  timef = proc.time()
  
  # Return the estimated models
  list(y = y, X = X, Z = Z, idx = idx,
       prior = prior, control = control, 
       burn = burn, trace = trace,
       exe.time = (timef - time0)[3])
}

mgp.gam.fit2 = function (y, X, Z, prior = list(), control = list()) {
  
  # This second implementation adaptively prunes the irrelevant basis along the run
  # y: n x 1 vector of response variables
  # X: n x p1 matrix of covariates
  # Z: list of n x ph basis matrices (h = 2, ..., H)
  
  # Get the initial CPU time
  time0 = proc.time()
  
  # Set the prior and control parameters
  prior = init.mgp.prior(prior)
  control = init.mgp.control(control)
  maxiter = control$burn + control$niter
  
  # Initialize all model parameters
  init.param = function (y, C, idx, prior) {
    n = length(y)
    p = unlist(lapply(idx, length))
    H = length(p)
    k = sum(p)
    
    # Init: beta
    lambda = rep(0, length = k)
    for (h in 2:H) {
      lambda[idx[[h]]] = rep(0.1 / p[h], p[h])
    }
    A = crossprod(C) + diag(lambda)
    b = crossprod(C, y)
    beta = rmnorm(A, b)
    
    # Init: eta
    eta = C %*% beta
    
    # Init: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum((y - eta)^2)
    psi = stats::rgamma(1, shape = ap, rate = ap)
    
    # Init: tau
    tau = rep(0, length = H)
    for (h in 2:H) {
      ih = idx[[h]]
      at = 0.5 * (prior$v + p[h])
      bt = 0.5 * (prior$v + sum(beta[ih]^2))
      tau[h] = stats::rgamma(1, shape = at, rate = at)
    }
    
    # Init: delta and phi
    delta = rep(0, length = k)
    phi = rep(0, length = k)
    for (h in 2:H) {
      ph = p[h]
      ih = idx[[h]]
      ad = prior$a1 + 0.5 * ph
      bd = 1 + 0.5 * tau[h] * beta[ih]^2
      deltah = rep(NA, length = ph)
      deltah[1] = stats::rgamma(1, shape = ad, rate = ad)
      for (jh in 2:ph) {
        ad = prior$a2 + 0.5 * (ph - jh + 1)
        bd = 1 + 0.5 * tau[h] * beta[ih]^2
        deltah[jh] = stats::rgamma(1, shape = ad, rate = ad)
      }
      delta[ih] = deltah
      phi[ih] = cumprod(deltah)
    }
    
    # output
    list(beta = beta, delta = delta, phi = phi, tau = tau, psi = psi)
  }
  
  # Log-prior density function
  get.logprior = function (beta, delta, tau, psi, idx, q, prior) {
    out = 0
    out = out + stats::dgamma(psi, shape = prior$a, rate = prior$b, log = TRUE)
    for (h in 2:length(idx)) {
      qh = q[h]
      ih = idx[[h]][1:qh]
      tauh = tau[h]
      deltah = delta[ih]
      betah = beta[ih]
      lambdah = tauh * cumprod(deltah)
      out = out + stats::dgamma(tauh, shape = prior$v / 2, rate = prior$v / 2, log = TRUE)
      out = out + stats::dgamma(deltah[1], shape = prior$a1, rate = 1, log = TRUE)
      out = out + sum(stats::dgamma(deltah[-1], shape = prior$a2, rate = 1, log = TRUE))
      out = out + sum(stats::dnorm(betah, mean = 0, sd = 1 / sqrt(lambdah), log = TRUE))
    }
    return (out)
  }
  
  # Log-likelihood function
  get.loglik = function (y, eta, psi) {
    out = sum(stats::dnorm(y, mean = eta, sd = sqrt(1 / psi), log = TRUE))
    return (out)
  }
  
  # Print simulation status
  print.status = function (iter, maxiter, keep, logpost, report, burn) {
    if (iter %% report == 0) {
      cat(gettextf(" iter: %5d /%5d  ", iter, maxiter),
          gettextf(" ncomp: %2d  ", length(keep)),
          gettextf(" logp: %.4f \n", logpost))
    }
    if (iter == burn) {
      cat(c(rep("-", 53), "\n"), sep = "")
    }
  }
  
  # Data dimenstions
  n = nrow(X)
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
              delta    = matrix(NA, nrow = control$burn, ncol = k),
              phi      = matrix(NA, nrow = control$burn, ncol = k),
              tau      = matrix(NA, nrow = control$burn, ncol = H),
              psi      = rep(NA, length = control$burn),
              logprior = rep(NA, length = control$burn),
              loglik   = rep(NA, length = control$burn),
              logpost  = rep(NA, length = control$burn),
              npar     = rep(NA, length = control$burn))
  
  trace = list(beta     = matrix(NA, nrow = control$niter, ncol = k),
               delta    = matrix(NA, nrow = control$niter, ncol = k),
               phi      = matrix(NA, nrow = control$niter, ncol = k),
               tau      = matrix(NA, nrow = control$niter, ncol = H),
               psi      = rep(NA, length = control$niter),
               logprior = rep(NA, length = control$niter),
               loglik   = rep(NA, length = control$niter),
               logpost  = rep(NA, length = control$niter),
               npar     = rep(NA, length = control$niter))
  
  # Initialize the unknown parameters
  init  = init.param(y, C, idx, prior)
  beta  = init$beta
  delta = init$delta
  tau   = init$tau
  psi   = init$psi
  phi   = init$phi
  eta   = C %*% beta
  
  # Compute the log prior, likelihood and posterior
  logprior = get.logprior(beta, delta, tau, psi, idx, q, prior)
  loglik = get.loglik(y, eta, psi)
  logpost = logprior + loglik
  
  # Store the initial values
  burn$beta[1,]    = beta
  burn$delta[1,]   = delta
  burn$phi[1,]     = phi
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
    D = get.penalty(tau, phi, idx)
    A = psi * CC[keep,keep] + diag(D[keep])
    b = psi * drop(Cy)[keep]
    beta[keep] = rmnorm(A, b)
    beta[drop] = 0
    
    # Update: eta
    eta = drop(C[,keep] %*% beta[keep])
    
    # Update: psi
    ap = prior$a + 0.5 * n
    bp = prior$b + 0.5 * sum((y - eta)^2)
    psi = stats::rgamma(1, shape = ap, rate = bp)
    
    # Update: tau, delta & phi
    for (h in 2:H) {
      # Set the block-specific parameters
      ih = idx[[h]]
      qh = q[h]
      ph = length(ih)
      betah = beta[ih]
      deltah = delta[ih]
      phih = phi[ih]
      tauh = tau[h]
      
      # Keep the relevant variables
      keeph = 1:qh
      droph = setdiff(1:ph, keeph)
      
      # Update tau_h
      at = 0.5 * (prior$v + qh)
      bt = 0.5 * (prior$v + sum(phih[keeph] * betah[keeph]^2))
      tauh = stats::rgamma(1, shape = at, rate = bt)
      
      # Update delta_1h
      ad = prior$a1 + 0.5 * qh
      bd = 1 + 0.5 * tauh * sum(phih[keeph] * betah[keeph]^2) / deltah[1]
      deltah[1] = stats::rgamma(1, shape = ad, rate = bd)
      phih[keeph] = cumprod(deltah[keeph])
      
      # Update delta_jh (j = 2, ..., q)
      for (j in 2:qh) {
        ad = prior$a2 + 0.5 * (qh - j + 1)
        bd = 1 + 0.5 * tauh * sum(phih[j:qh] * betah[j:qh]^2) / deltah[j]
        deltah[j] = stats::rgamma(1, shape = ad, rate = bd)
        phih[keeph] = cumprod(deltah[keeph])
      }
      
      # Fill the remainder components
      deltah[droph] = Inf
      phih[droph] = Inf
      
      # Update delta and phi
      tau[h] = tauh
      delta[ih] = deltah
      phi[ih] = phih
    }
    
    # Compute the log prior, likelihood and posterior
    logprior = get.logprior(beta, delta, tau, psi, idx, q, prior)
    loglik = get.loglik(y, eta, psi)
    logpost = logprior + loglik
    
    # Store the initial values
    if (iter <= control$burn) {
      t = iter
      burn$beta[t,]    = beta
      burn$delta[t,]   = delta
      burn$phi[t,]     = phi
      burn$tau[t,]     = tau
      burn$psi[t]      = psi
      burn$logprior[t] = logprior
      burn$loglik[t]   = loglik
      burn$logpost[t]  = logpost
      burn$npar[t]     = length(keep)
    } else {
      t = iter - control$burn
      trace$beta[t,]    = beta
      trace$delta[t,]   = delta
      trace$phi[t,]     = phi
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
            delta[ih] = Inf
            phi[ih] = Inf
          }
          # If the explained variance is high, we add a new component
          # and we sample the associated parameters from the prior
          if (qh < p[h] && expvar > control$tol) {
            q[h] = qh + 1
            beta[ih+1] = mean(C[,ih] * res) / mean(C[,ih]^2)
            delta[ih+1] = delta[ih]
            phi[ih+1] = phi[ih] * delta[ih]
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
    idx = seq(from = 1, to = control$burn, by = control$thin)
    burn$beta     = burn$beta[idx,]
    burn$delta    = burn$delta[idx,]
    burn$phi      = burn$phi[idx,]
    burn$tau      = burn$tau[idx,]
    burn$psi      = burn$psi[idx]
    burn$logprior = burn$logprior[idx]
    burn$loglik   = burn$loglik[idx]
    burn$logpost  = burn$logpost[idx]
    burn$npar     = burn$npar[idx]
    
    idx = seq(from = 1, to = control$niter, by = control$thin)
    trace$beta     = trace$beta[idx,]
    trace$delta    = trace$delta[idx,]
    trace$phi      = trace$phi[idx,]
    trace$tau      = trace$tau[idx,]
    trace$psi      = trace$psi[idx]
    trace$logprior = trace$logprior[idx]
    trace$loglik   = trace$loglik[idx]
    trace$npar     = trace$npar[idx]
  }
  
  # Get the final CPU time
  timef = proc.time()
  
  # Return the estimated models
  list(y = y, X = X, Z = Z, idx = idx,
       prior = prior, control = control, 
       burn = burn, trace = trace,
       exe.time = (timef - time0)[3])
}




