#' file: utils_cusp_reg.R
#' author: Cristian Castiglione
#' creation: 05/01/2024
#' last change: 05/01/2024


#' Function for initializing the prior parameters of the MGP
#' @keywords internal
init.cusp.prior = function (prior) {
  default = list(
    a = 2.1,
    b = 1.1,
    at = 2.1,
    bt = 1.1,
    t0 = 1e-03,
    alpha = 5)
  
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
  if (default$at <= 0) stop("The hyperparameter `at` must be a positive real number.")
  if (default$bt <= 0) stop("The hyperparameter `bt` must be a positive real number.")
  if (default$t0 <= 0) stop("The hyperparameter `t0` must be a positive real number.")
  if (default$alpha <= 0) stop("The hyperparameter `alpha` must be a positive real number.")
  if (default$bt / default$at <= default$t0) stop("The increasing shrinkage condition is not satisfied.")
  
  # Return the control parameters
  return (default)
}

#' Function for initializing the control parameters of the algorithms
#' @keywords internal
init.cusp.control = function (control) {
  
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
init.cusp.param = function (y, C, idx, prior) {
  n = length(y)
  p = unlist(lapply(idx, length))
  H = length(p)
  k = sum(p)
  
  # Init: beta
  lambda = rep(0, length = k)
  for (h in 2:H) {
    ph = p[h]
    ih = idx[[h]]
    lambda[ih] = rep(0.1 / ph, ph)
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
  
  # Init: omega and pi
  theta = rep(0, length = k)
  omega = rep(0, length = k)
  pi = rep(0, length = k)
  nu = rep(0, length = k)
  tau = rep(0, length = H)
  
  theta[idx[[1]]] = Inf
  tau[1] = Inf
  for (h in 2:H) {
    ph = p[h]
    ih = idx[[h]]
    tau[h] = stats::rgamma(1, shape = prior$at + 0.5 * ph, rate = prior$bt + 0.5 * sum(beta[ih]^2))
    nu[ih] = stats::rbeta(ph, shape = 1, shape2 = prior$alpha)
    nu[ih[ph]] = 1
    omega[ih] = nu[ih] * cumprod(c(1, 1 - nu[ih[-1]]))
    pi[ih] = cumsum(omega[ih])
    theta[ih] = pi[ih] * prior$t0 + (1 - pi[ih]) * prior$bt / (prior$at - 1)
  }
  
  # output
  list(beta = beta, theta = theta, pi = pi, 
       omega = omega, nu = nu, tau = tau, psi = psi)
}


#' Simulate from the prior distribution
#' @keywords internal
cusp.prior.sim = function (npar, niter, prior) {
  psi = stats::rgamma(niter, shape = prior$a, rate = prior$b)
  tau = stats::rgamma(niter, shape = prior$at, rate = prior$bt)
  nu = matrix(NA, nrow = niter, ncol = npar)
  omega = matrix(NA, nrow = niter, ncol = npar)
  pi = matrix(NA, nrow = niter, ncol = npar)
  theta = matrix(NA, nrow = niter, ncol = npar)
  for (h in 1:niter) {
    nu[h,] = stats::rbeta(npar, shape1 = 1, shape2 = prior$alpha)
    nu[h,npar] = 1
    omega[h,] = nu[h,] * cumprod(c(1, 1 - nu[h,-1]))
    pi[h,] = cumsum(omega[h,])
    theta[h,] = pi[h,] * prior$t0 + (1 - pi[h,]) / tau[h]
  }
  list(psi = psi, tau = tau, nu = nu, 
       omega = omega, pi = pi, theta = theta)
}

#' @title Compute the CUSP prior mean
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
#' @return
#' A short description...
#' 
#' @export
cusp.prior.mean = function (t, prior, par = "theta") {
  alpha = prior$alpha
  spike = prior$t0
  slab = prior$bt / (prior$at-1)
  if (par == "nu") mn = 1 / (1 + alpha)
  if (par == "omega") mn = alpha^(t-1) / (1 + alpha)^t
  if (par == "pi") mn = 1 - alpha^t / (1 + alpha)^t
  if (par == "theta") mn = spike + alpha^t * (1 + alpha)^-t * (slab - spike)
  return (mn)
} 


