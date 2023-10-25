#' mgp_reg_utils.R
#' author: Cristian Castiglione
#' creation: 04/08/2023
#' last change: 16/09/2023

# Get the column-dimension of each block
get.block.dim = function (X, Z) {
  c(ncol(X), unlist(lapply(Z, ncol)))
}

# Get the column-indices of each block
get.block.idx = function (p) {
  n = length(p)
  start = cumsum(c(1, p[-n]))
  end = cumsum(p)
  idx = apply(cbind(start, end), 1, function(x) x[1]:x[2])
  return (idx)
}

# Get the column-indices of each block concatenated by row
get.flat.idx = function (idx, q) {
  keep = c()
  for (h in 1:length(q)) {
    keep = c(keep, idx[[h]][1:q[h]])
  }
  return (keep)
}

# Build the model matrix and the sufficient statistics
get.matrix = function (X, Z) {
  C = cbind(X, do.call(cbind, Z))
  return (C)
}

# Get the vector of penalty parameters
get.penalty = function (tau, phi, idx) {
  lambda = rep(0, length = length(phi))
  for (h in 2:length(idx)) {
    ih = idx[[h]]
    lambda[ih] = tau[h] * phi[ih]
  }
  return (lambda)
}

# Function for initializing the prior parameters of the MGP
init.mgp.prior = function (prior) {
  default = list(
    a = 2.1,
    b = 1.1,
    v = 0.1,
    a1 = 2.1,
    a2 = 3.1)
  
  if (class(prior) != "list") {
    # Print a warning if `control` is not a list
    warning("The `prior` parameter must be a list. \n",
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
init.mgp.control = function (control) {
  
  # Default control parameters
  default = list(
    niter = 5000,
    burn = 2500,
    thin = 1,
    adaptation = FALSE,
    tol = 1e-04,
    minnpar = 2,
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

# Function for initializing the unknown parameters of the model
init.mgp.param = function (y, X, Z, prior) {
  n = nrow(X)
  p = ncol(X)
  q = ncol(Z)
  
  idx = 1:p
  idz = (p+1):(p+q)
  
  # Init: beta
  D = diag(c(rep(0, p), rep(0.1/q, q)))
  A = crossprod(cbind(X, Z)) + D
  b = crossprod(cbind(X, Z), y)
  beta = rmnorm(A, b)
  
  # Init: psi
  ap = prior$a + 0.5 * n
  bp = prior$b + 0.5 * sum((y - X %*% beta[idx] - Z %*% beta[idz])^2)
  psi = stats::rgamma(1, shape = ap, rate = ap)
  
  # Init: tau
  at = 0.5 * (prior$v + q)
  bt = 0.5 * (prior$v + sum(beta[idz]^2))
  tau = stats::rgamma(1, shape = at, rate = at)
  
  # Init: delta
  ad = prior$a1 + 0.5 * q
  bd = 1 + 0.5 * tau * beta[idz[1]]^2
  delta = rep(NA, length = q)
  delta[1] = stats::rgamma(1, shape = ad, rate = ad)
  
  for (h in 2:q) {
    ad = prior$a2 + 0.5 * (q - h + 1)
    bd = 1 + 0.5 * tau * beta[idz[h]]^2
    delta[h] = stats::rgamma(1, shape = ad, rate = ad)
  }
  
  # output
  list(beta = beta, delta = delta, tau = tau, psi = psi)
}

# Simulate from the prior distribution
mgp.prior.sim = function (npar, niter, prior) {
  psi = stats::rgamma(niter, shape = prior$a, rate = prior$b)
  tau = stats::rgamma(niter, shape = 0.5 * prior$v, rate = 0.5 * prior$v)
  phi = matrix(NA, nrow = niter, ncol = npar)
  delta = matrix(NA, nrow = niter, ncol = npar)
  delta[,1] = stats::rgamma(niter, shape = prior$a1, rate = 1)
  delta[,2:npar] = stats::rgamma(niter * (npar - 1), shape = prior$a2, rate = 1)
  logphi = t(apply(log(delta), 1, cumsum))
  list(psi = psi, tau = tau, delta = delta, logphi = logphi)
}

# Posterior coefficients
coef.mgp.reg = function(object) {
  beta  = colMeans(object$trace$beta)
  phi   = colMeans(object$trace$phi)
  delta = colMeans(object$trace$delta)
  tau   = mean(object$trace$tau)
  psi   = mean(object$trace$psi)
  
  list(beta = beta, phi = phi, delta = delta, tau = tau, psi = psi)
}

# Posterior variances
var.mgp.reg = function (object) {
  beta  = apply(object$trace$beta, 2, var)
  phi   = apply(object$trace$phi, 2, var)
  delta = apply(object$trace$delta, 2, var)
  tau   = mean(object$trace$tau)
  psi   = mean(object$trace$psi)
  
  list(beta = beta, phi = phi, delta = delta, tau = tau, psi = psi)
}

# Predicted values
predict.mgp.reg = function (object, newdata = NULL) {
  yhat = NULL
  if (is.null(newdata)) {
    yhat = tcrossprod(cbind(object$X, object$Z), object$trace$beta)
  } else {
    yhat = tcrossprod(cbind(newdata$X, newdata$Z), object$trace$beta)
  }
  list(
    mean = drop(rowMeans(yhat)),
    sd  = drop(apply(yhat, 1, sd)),
    q05 = drop(apply(yhat, 1, quantile, probs = 0.05)),
    q25 = drop(apply(yhat, 1, quantile, probs = 0.25)),
    q50 = drop(apply(yhat, 1, quantile, probs = 0.50)),
    q75 = drop(apply(yhat, 1, quantile, probs = 0.75)),
    q95 = drop(apply(yhat, 1, quantile, probs = 0.95)))
}

# Posterior summary statistics
summary.mgp.reg = function (object) {
  
  # beta
  .beta = object$trace$beta
  beta = data.frame(
    mean = c(colMeans(.beta)),
    sd  = c(apply(.beta, 2, sd)),
    q05 = c(apply(.beta, 2, quantile, probs = 0.05)),
    q25 = c(apply(.beta, 2, quantile, probs = 0.25)),
    q50 = c(apply(.beta, 2, quantile, probs = 0.50)),
    q75 = c(apply(.beta, 2, quantile, probs = 0.75)),
    q95 = c(apply(.beta, 2, quantile, probs = 0.95)))
  
  rownames(beta) = paste("beta[", 1:nrow(beta), "]", sep = "")
  
  # phi
  .invphi = 1 / object$trace$phi
  invphi = data.frame(
    mean = c(colMeans(.invphi)),
    q05 = c(apply(.invphi, 2, quantile, probs = 0.05)),
    q25 = c(apply(.invphi, 2, quantile, probs = 0.25)),
    q50 = c(apply(.invphi, 2, quantile, probs = 0.50)),
    q75 = c(apply(.invphi, 2, quantile, probs = 0.75)),
    q95 = c(apply(.invphi, 2, quantile, probs = 0.95)))
  
  rownames(invphi) = paste("1/phi[", 1:nrow(invphi), "]", sep = "")
  
  # tau
  .invtau = 1 / object$trace$tau
  invtau = data.frame(
    mean = c(mean(.invtau)),
    q05 = c(quantile(.invtau, probs = 0.05)),
    q25 = c(quantile(.invtau, probs = 0.25)),
    q50 = c(quantile(.invtau, probs = 0.50)),
    q75 = c(quantile(.invtau, probs = 0.75)),
    q95 = c(quantile(.invtau, probs = 0.95)))
  
  rownames(invtau) = c("1/tau")
  
  # psi
  .invpsi = 1 / object$trace$psi
  invpsi = data.frame(
    mean = c(mean(.invpsi)),
    q05 = c(quantile(.invpsi, probs = 0.05)),
    q25 = c(quantile(.invpsi, probs = 0.25)),
    q50 = c(quantile(.invpsi, probs = 0.50)),
    q75 = c(quantile(.invpsi, probs = 0.75)),
    q95 = c(quantile(.invpsi, probs = 0.95)))
  
  rownames(invpsi) = c("1/psi")
  
  # Print out the result
  cat(c(rep("-", 55), "\n"), sep = "")
  print(round(beta, 3))
  cat(c(rep("-", 55), "\n"), sep = "")
  print(round(invphi, 3))
  # cat(c(rep("-", 55), "\n"), sep = "")
  # print(round(delta, 3))
  cat(c(rep("-", 55), "\n"), sep = "")
  print(round(invtau, 3))
  cat(c(rep("-", 55), "\n"), sep = "")
  print(round(invpsi, 3))
  cat(c(rep("-", 55), "\n"), sep = "")
}

