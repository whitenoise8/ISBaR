#' mgp_reg_sim.R
#' author: Cristian Castiglione
#' creation: 15/09/2023
#' last change: 18/09/2023

rm(list = ls())
graphics.off()

## Packages ----
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgcv)

source("spline_basis.R")
source("fourier_basis.R")
source("mgp_reg_utils.R")
source("blm_reg_fit.R")
source("mgp_reg_fit.R")
source("hs_reg_fit.R")

NSIM = 100
NOBS = 250
NBASIS = 30

## Simulation functions ----
data.simulation = function (nsim = 100, nobs = 250, nbasis = 20, fun = NULL) {
  
  # Support limit
  # a = 0
  # b = 2
  
  # True signal
  f = NULL
  if (!is.null(fun) & is.function(fun)) {
    f = fun
  } else {
    # a = 0; b = 2
    # f = function (x) 3 * cos(2 * pi * x) / (2 * pi * x + pi) + 1 + 0.5 * x
    
    # a = -2; b = +2
    # f = function (x) exp(- x**2) * (cos(0.5 * pi * x + 0.3 * pi))**2
    
    a = -2; b = +2
    f = function (x) exp(- 1.25 * x**2) * ((x - 0.25)**2 + 0.3)
  }
  
  # Data storing
  Xn = matrix(NA, nrow = nsim, ncol = nobs)
  Fn = matrix(NA, nrow = nsim, ncol = nobs)
  Yn = matrix(NA, nrow = nsim, ncol = nobs)
  
  Kn = matrix(NA, nrow = nsim, ncol = nbasis)
  Pn = array(NA, dim = c(nsim, nbasis+4, nbasis+4))
  Zn = array(NA, dim = c(nsim, nobs, nbasis+2))
  
  # Simulation loop
  for (sim in 1:nsim) {
    
    # Contaminated signal simulation
    xn = sort(runif(nobs, min = a, max = b))
    en = rnorm(nobs, mean = 0, sd = 0.075)
    fn = f(xn)
    yn = fn + en
    
    # B-spline basis and penalty matrices
    kn = get.bspline.knots(xn, nbasis)
    bn = get.bspline.matrix(xn, kn, a, b)
    pn = get.bspline.penalty(kn, a, b)
    zn = get.ospline.matrix(xn, bn, pn)$z
    
    # Data storing
    Xn[sim, ] = xn
    Fn[sim, ] = fn
    Yn[sim, ] = yn
    Kn[sim, ] = kn
    Pn[sim, , ] = pn
    Zn[sim, , ] = zn
  }
  
  # Returned values
  list(time = Xn,
       signal = Fn,
       data = Yn,
       knots = kn,
       penalty = Pn,
       matrix = Zn)
}

run.simulation = function (nsim = 100, nobs = 250, nbasis = 20, fun = NULL) {
  
  time0 = proc.time()
  
  # Simulation parameters
  alphas = c(.025, .25, .50, .75, .975)
  
  # MCMC parameters
  prior = list(a1 = 2.1, a2 = 3.1, v = 0.1, a = 0.1, b = 0.1)
  control = list(burn = 1000, niter = 5000, thin = 5, verbose = FALSE)
  nmc = floor(control$niter / control$thin)
  
  # Error function
  rmse = function (x, y) sqrt(mean((x - y)^2))
  scale = function (x) (x - mean(x)) / sd(x)
  
  # Data simulation and basis matrix construction
  dat = data.simulation(nsim = nsim, nobs = nobs, nbasis = nbasis, fun = fun)
  
  # Memory allocation
  error = array(NA, dim = c(nsim, 4, nmc))
  
  # Initialize the simulation bar
  progress.bar(0, 0, 50)

  # Simulation loop
  for (sim in 1:nsim) {
    
    # Simulated data
    y = dat$data[sim, ]
    x = dat$time[sim, ]
    f = dat$signal[sim, ]
    X = cbind(1, scale(x))
    Z = dat$matrix[sim, , ]
    
    # Model fit
    blm = blm.reg.fit(y, X, Z, prior = list(), control = control)
    hs  =  hs.reg.fit(y, X, Z, prior = list(), control = control)
    mgp = mgp.reg.fit(y, X, Z, prior =  prior, control = control)
    rml = mgcv::gam(y ~ s(x, bs = "cr", k = nbasis, m = 2), 
                    data = data.frame(x = x, y = y),
                    family = gaussian(), method = "REML")
    
    # Predictions
    f.hat.blm = tcrossprod(cbind(X, Z), blm$trace$beta)
    f.hat.hs  = tcrossprod(cbind(X, Z),  hs$trace$beta)
    f.hat.mgp = tcrossprod(cbind(X, Z), mgp$trace$beta)
    f.hat.rml = mgcv::predict.gam(rml, se.fit = FALSE)
    
    # Reconstruction error
    error[sim, 1, ] = apply(f.hat.blm, 2, function (x) rmse(x, f)**2)
    error[sim, 2, ] = apply(f.hat.hs , 2, function (x) rmse(x, f)**2)
    error[sim, 3, ] = apply(f.hat.mgp, 2, function (x) rmse(x, f)**2)
    error[sim, 4, ] = rmse(f.hat.rml, f)**2
    
    # Plot the posterior distribution of the RMSE at each iteration
    df = data.frame(model = rep(c("BLM", "HS", "MGP"), each = nmc), 
                    RMSE = sqrt(c(error[sim, 1, ], error[sim, 2, ], error[sim, 3, ])))
    
    plt = ggplot(data = df, mapping = aes(x = RMSE, fill = model, color = model)) +
      geom_density(alpha = .2) +
      geom_vline(xintercept = mean(sqrt(error[sim, 1, ])), lwd = 0.8, lty = 2, color = 2) + # BLM
      geom_vline(xintercept = mean(sqrt(error[sim, 2, ])), lwd = 0.8, lty = 2, color = 3) + # HS
      geom_vline(xintercept = mean(sqrt(error[sim, 3, ])), lwd = 0.8, lty = 2, color = 4) + # MGP
      geom_vline(xintercept = mean(sqrt(error[sim, 4, ])), lwd = 0.8, lty = 2, color = 6) + # REML
      ggtitle(paste("Posterior distribution of the RMSE (", sim, ")", sep = "")) + 
      theme_bw()

    print(plt)
    
    # Print the simulation bar
    time = (proc.time() - time0)[3]
    progress = floor(50 * sim / nsim)
    progress.bar(progress, time, 50)
  }
  
  # Return the signal reconstruction error
  return (error)
}

progress.bar = function (progress, time, length = 50) {
  if (time > 60) {
    time = paste(floor(time / 60), "m", sep = "")
  } else {
    time = paste(floor(time), "s", sep = "")
  }
  if (progress < length) {
    bar = c(rep("=", progress), ">", rep(".", length - progress - 1))
    cat("\r |", bar, "| ", 2 * progress, "%, ", time, "  ", sep = "")
  } else {
    cat("\r |", rep("=", progress), "| 100%, ", time, "\n", sep = "")
  }
}

## Main execution ----
main = function () {
  
  nsim = NSIM
  nobs = NOBS
  nbasis = NBASIS
  
  # Run the simulation and collect the reconstruction error
  error = run.simulation(nsim = nsim, nobs = nobs, nbasis = nbasis, fun = NULL)
  error.matrix = sqrt(rowMeans(error, dims = 2))
  
  # Plot the distribution of the reconstruction error
  df = data.frame(
    model = rep(c("BLM", "HS", "MGP", "REML"), each = nsim),
    error = c(error.matrix) # c(error.matrix[,1], error.matrix[,2], error.matrix[,3], error.matrix[,4])
  )
  
  plt = ggplot(data = df, mapping = aes(x = model, y = error, col = model, fill = model)) + 
    geom_boxplot(alpha = 0.4) + theme_bw() + 
    labs(x = "Model", y = "RMSE", title = "Average RMSE a posteriori")
  
  print(plt)

  # Return the plot object
  return (plt)
}

main()

## End of file ----
