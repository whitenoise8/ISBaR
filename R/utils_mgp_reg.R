#' file: utils_mgp_reg.R
#' author: Cristian Castiglione
#' creation: 04/08/2023
#' last change: 16/09/2023

#' Function for initializing the prior parameters of the MGP
#' @keywords internal
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
  
  # Safety checks
  if (!all(unlist(lapply(default, is.numeric)))) stop("The prior hyperparameters are not numeric.")
  if (default$a <= 0) stop("The hyperparameter `a` must be a positive real number.")
  if (default$b <= 0) stop("The hyperparameter `b` must be a positive real number.")
  if (default$v <= 0) stop("The hyperparameter `v` must be a positive real number.")
  if (default$a1 <= 0) stop("The hyperparameter `a1` must be a positive real number.")
  if (default$a2 <= 0) stop("The hyperparameter `a2` must be a positive real number.")
  if (default$a1 >= default$a2) stop("The increasing shrinkage condition is not satisfied.")
  
  # Return the control parameters
  return (default)
}

#' Function for initializing the control parameters of the algorithms
#' @keywords internal
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

#' Function for initializing the unknown parameters of the model
#' @keywords internal
init.mgp.param = function (y, C, idx, prior) {
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

#' Simulate from the prior distribution
#' @keywords internal
mgp.prior.sim = function (npar, niter, prior) {
  psi = stats::rgamma(niter, shape = prior$a, rate = prior$b)
  tau = stats::rgamma(niter, shape = 0.5 * prior$v, rate = 0.5 * prior$v)
  phi = matrix(NA, nrow = niter, ncol = npar)
  delta = matrix(NA, nrow = niter, ncol = npar)
  delta[,1] = stats::rgamma(niter, shape = prior$a1, rate = 1)
  delta[,2:npar] = stats::rgamma(niter * (npar - 1), shape = prior$a2, rate = 1)
  logphi = t(apply(log(delta), 1, cumsum))
  logtheta = matrix(NA, nrow = niter, ncol = npar)
  for (h in 1:npar) {logtheta[,h] = - logphi[,h] - log(tau)}
  list(psi = psi, tau = tau, delta = delta, logphi = logphi, logtheta = logtheta)
}


#' @title Compute the MGP prior mean
#' 
#' @description 
#' A short description...
#' 
#' @param t ...
#' @param prior ...
#' 
#' @details
#' Additional details...
#' 
#' 
#' @return
#' A short description...
#' 
#' @export
mgp.prior.mean = function (t, prior) {
  a1 = prior$a1
  a2 = prior$a2
  nu = prior$nu
  mn = nu * (nu-2)^-1 * (a1-1)^-1 * (a2-1)^-(t-1)
  return (mn)
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

