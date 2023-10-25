#' mgp_reg_test.R
#' author: Cristian Castiglione
#' creation: 05/08/2023
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
source("fast_rmnorm.R")
source("mgp_reg_utils.R")
source("blm_reg_fit.R")
source("mgp_reg_fit.R")
source("hs_reg_fit.R")

## Data simulation ----
n = 250; p = 20; a = 0; b = 2
x = sort(runif(n, min = a, max = b))
f = 3 * cos(2 * pi * x) / (2 * pi * x + pi) + 1 + 0.5 * x
e = rnorm(n, mean = 0, sd = 0.1)
y = f + e

plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)

n = 250; p = 20; a = -2; b = +2
x = sort(runif(n, min = a, max = b))
f = exp(-x**2) * (cos(0.5 * pi * x + 0.3 * pi))**2
e = rnorm(n, mean = 0, sd = 0.1)
y = f + e

plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)

n = 250; p = 40; a = -2; b = +2
x = sort(runif(n, min = a, max = b))
f = exp(-1.25*x**2) * ((x-0.25)**2 + 0.3)
e = rnorm(n, mean = 0, sd = 0.075)
y = f + e

plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)


## Spline regression ----

# Compute the basis and penalty matrices
k = get.bspline.knots(x, p)
B = get.bspline.matrix(x, k, a, b)
P = get.bspline.penalty(k, a, b)
Z = get.ospline.matrix(x, B, P)$z
X = cbind(1, x-mean(x))
C = cbind(X, Z)

# Fit a Bayesian linear model with standard prior using MCMC
blm.prior = list(v = 0.1, a = 0.1, b = 0.1)
blm.control = list(burn = 1000, niter = 5000, report = 500, verbose = TRUE, thin = 5)
blm.fit = blm.reg.fit(y, X, Z, blm.prior, blm.control)
blm.pred = tcrossprod(C, blm.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Fit a Bayesian linear models with Horseshoe prior using MCMC
hs.control = list(burn = 1000, niter = 5000, thin = 5)
hs.fit = hs.reg.fit(y, X, Z, prior = list(), control = hs.control)
hs.pred = tcrossprod(C, hs.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Fit a Bayesian linear model with MGP prior using MCMC
mgp.prior = list(a1 = 2.1, a2 = 3.1, v = 0.1, a = 0.1, b = 0.1)
mgp.control = list(burn = 1000, niter = 5000, tol = 1e-05, adaptation = TRUE, report = 500, verbose = TRUE, thin = 5)
mgp.fit = mgp.reg.fit2(y, X, Z, mgp.prior, mgp.control)
mgp.pred = tcrossprod(C, mgp.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Fit a cubic regression spline using PIRLS + REML
mgcv.fit = mgcv::gam(y ~ s(x, bs = "cr", k = p, m = 2), 
                     data = data.frame(x = x, y = y),
                     family = gaussian(), method = "REML")
mgcv.pred = mgcv::predict.gam(mgcv.fit, se.fit = FALSE)
mgcv.se = mgcv::predict.gam(mgcv.fit, se.fit = TRUE)$se.fit
mgcv.pred = cbind(mgcv.pred - 1.96 * mgcv.se, mgcv.pred, mgcv.pred + 1.96 * mgcv.se)

# Compare the estimated curves (posterior mean)
{
  plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)
  lines(x, blm.pred[,1], col = 3, lty = 2, lwd = 1)
  lines(x, blm.pred[,2], col = 3, lty = 1, lwd = 2)
  lines(x, blm.pred[,3], col = 3, lty = 2, lwd = 1)
  lines(x, mgp.pred[,1], col = 4, lty = 2, lwd = 1)
  lines(x, mgp.pred[,2], col = 4, lty = 1, lwd = 2)
  lines(x, mgp.pred[,3], col = 4, lty = 2, lwd = 1)
  lines(x, hs.pred[,1], col = 7, lty = 2, lwd = 1)
  lines(x, hs.pred[,2], col = 7, lty = 1, lwd = 2)
  lines(x, hs.pred[,3], col = 7, lty = 2, lwd = 1)
  lines(x, mgcv.pred[,1], col = 5, lty = 2, lwd = 1)
  lines(x, mgcv.pred[,2], col = 5, lty = 1, lwd = 2)
  lines(x, mgcv.pred[,3], col = 5, lty = 2, lwd = 1)
  legend("topright", col = c(2, 3, 4, 7, 5), lty = 1, lwd = 2,
         legend = c("TRUE", "BLM", "MGP", "HS", "MGCV"))
}

{
  df = data.frame(
    x = rep(x, times = 5),
    y = c(f, blm.pred[,2], hs.pred[,2], mgp.pred[,2], mgcv.pred[,2]),
    ylo = c(f, blm.pred[,1], hs.pred[,1], mgp.pred[,1], mgcv.pred[,1]),
    yup = c(f, blm.pred[,3], hs.pred[,3], mgp.pred[,3], mgcv.pred[,3]),
    method = rep(c("TRUE", "BLM", "HS", "MGP", "MGCV"), each = n))
  
  df$method = factor(df$method, levels = c("TRUE", "BLM", "HS", "MGP", "MGCV"))
  
  plt = ggplot(data = df) + theme_bw() + labs(x = "", y = "") +
    geom_point(data = data.frame(x = x, y = y), mapping = aes(x = x, y = y), alpha = 0.4) +
    geom_line(mapping = aes(x = x, y = y, color = method), lwd = 0.8) + 
    geom_ribbon(mapping = aes(x = x, ymax = ylo, ymin = yup, fill = method), alpha = 0.3) +
    facet_wrap(~ method)
  
  print(plt)
}

# Let's have a look to the posterior of the (invese) shrinkage factors
{
  sigma = 1 / sqrt(mgp.fit$trace$tau * mgp.fit$trace$phi)
  sigma = t(apply(sigma, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE))
  matplot(sigma, type = "o", col = 2:6, lty = 1, pch = 20, 
          xlab = "j", ylab = expression(sigma[j]), 
          main = expression(paste(
            "Posterior distribution of ", 
            sigma[j], " = ", (tau * phi[j])^"-1/2")))
  abline(h = median(1 / sqrt(blm.fit$trace$tau)), col = 1, lty = 2)
  legend("topright", col = 2:6, lty = 1, pch = 20, 
         legend = c("5%", "25%", "50%", "75%", "95%"))
}

# Let's have a look to the posterior of the regression coefficients
{
  beta.blm = t(apply(blm.fit$trace$beta, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  beta.mgp = t(apply(mgp.fit$trace$beta, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  beta.hs = t(apply(hs.fit$trace$beta, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  plot(0, 0, type = "n", xlim = c(1, p+4), ylim = c(-1,+1) * 2.5,
       xlab = "j", ylab = expression(beta[j]), 
       main = expression(paste(
         "Posterior credibility intervals of ", beta[j], " (j = 1, ..., p)")))
  abline(h = 0, lty = 2, col = 8)
  # Bayesian linear model
  points(1:(p+4), beta.blm[,2], col = 4, pch = 19)
  arrows(1:(p+4), beta.blm[,1], 1:(p+4), beta.blm[,3], 
         angle = 90, code = 3, length = 0.02, col = 4)
  # Multiplicative Gamma process
  points(1:(p+4)+0.2, beta.mgp[,2], col = 2, pch = 19)
  arrows(1:(p+4)+0.2, beta.mgp[,1], 1:(p+4)+0.2, beta.mgp[,3], 
         angle = 90, code = 3, length = 0.02, col = 2)
  # Horseshoe prior
  points(1:(p+4)+0.4, beta.hs[,2], col = 3, pch = 19)
  arrows(1:(p+4)+0.4, beta.hs[,1], 1:(p+4)+0.4, beta.hs[,3], 
         angle = 90, code = 3, length = 0.02, col = 3)
  legend("bottomleft", col = c(4, 2, 3), lty = 1, lwd = 2,
         legend = c("BLM", "MGP", "HS"))
}

# Let's have a look to the posterior correlation matrix of beta
corrplot::corrplot(cor(blm.fit$trace$beta), method = "color")
corrplot::corrplot(cor( hs.fit$trace$beta), method = "color")
corrplot::corrplot(cor(mgp.fit$trace$beta), method = "color")

# Signal recontruction error
rmse = function (x, y = 0) sqrt(mean((x - y)^2))

rmse.blm = tcrossprod(C, blm.fit$trace$beta) |> apply(2, function(x) rmse(x, f))
rmse.mgp = tcrossprod(C, mgp.fit$trace$beta) |> apply(2, function(x) rmse(x, f))
rmse.hs = tcrossprod(C, hs.fit$trace$beta) |> apply(2, function(x) rmse(x, f))
rmse.mgcv = rmse(mgcv.pred[,2], f)

df = data.frame(model = rep(c("BLM", "MGP", "HS"), each = 1000), 
                RMSE = c(rmse.blm, rmse.mgp, rmse.hs))

ggplot(data = df, mapping = aes(x = RMSE, fill = model, color = model)) +
  geom_density(alpha = .2) +
  geom_vline(xintercept = sqrt(mean(rmse.blm**2)), lwd = 0.8, lty = 2, color = 2) +
  geom_vline(xintercept = sqrt(mean(rmse.hs **2)), lwd = 0.8, lty = 2, color = 3) +
  geom_vline(xintercept = sqrt(mean(rmse.mgp**2)), lwd = 0.8, lty = 2, color = 4) +
  geom_vline(xintercept = rmse.mgcv, lwd = 0.8, lty = 2, color = 6) +
  ggtitle("Posterior distribution of the RMSE") + 
  theme_bw()



## Fourier regression ----

# Compute the basis and penalty matrices
Z = get.fourier.matrix(x, p)
X = cbind(1, x-1)
C = cbind(X, Z)

# Fit a Bayesian linear model with standard prior using MCMC
blm.prior = list(v = 0.1, a = 0.1, b = 0.1)
blm.control = list(burn = 1000, niter = 5000, report = 500, verbose = TRUE, thin = 5)
blm.fit = blm.reg.fit(y, X, Z, blm.prior, blm.control)
blm.pred = tcrossprod(C, blm.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Fit a Bayesian linear model with MGP prior using MCMC
mgp.prior = list(a1 = 2.1, a2 = 3.1, v = 0.1, a = 0.1, b = 0.1)
mgp.control = list(burn = 1000, niter = 5000, report = 500, verbose = TRUE, thin = 5)
mgp.fit = mgp.reg.fit(y, X, Z, mgp.prior, mgp.control)
mgp.pred = tcrossprod(C, mgp.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Compare the estimated curves (posterior mean)
{
  plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)
  lines(x, blm.pred[,1], col = 3, lty = 2, lwd = 1)
  lines(x, blm.pred[,2], col = 3, lty = 1, lwd = 2)
  lines(x, blm.pred[,3], col = 3, lty = 2, lwd = 1)
  lines(x, mgp.pred[,1], col = 4, lty = 2, lwd = 1)
  lines(x, mgp.pred[,2], col = 4, lty = 1, lwd = 2)
  lines(x, mgp.pred[,3], col = 4, lty = 2, lwd = 1)
  legend("bottomright", col = c(2, 3, 4), lty = 1, lwd = 2,
         legend = c("TRUE", "BLM", "MGP"))
  }

# Let's have a look to the posterior of the (invese) shrinkage factors
{
  sigma = 1 / sqrt(mgp.fit$trace$tau * mgp.fit$trace$phi)
  sigma = t(apply(sigma, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
  matplot(sigma, type = "o", col = 2:6, lty = 1, pch = 20, 
          xlab = "j", ylab = expression(sigma[j]), 
          main = expression(paste(
            "Posterior distribution of ", 
            sigma[j], " = ", (tau * phi[j])^"-1/2")))
  abline(h = median(1 / sqrt(blm.fit$trace$tau)), col = 1, lty = 2)
  legend("topright", col = 2:6, lty = 1, pch = 20, 
         legend = c("5%", "25%", "50%", "75%", "95%"))
}

# Let's have a look to the posterior of the regression coefficients
{
  beta.blm = t(apply(blm.fit$trace$beta, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  beta.mgp = t(apply(mgp.fit$trace$beta, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  plot(0, 0, type = "n", xlim = c(1, p+2), ylim = c(-1,+1) * 30,
       xlab = "j", ylab = expression(beta[j]), 
       main = expression(paste(
         "Posterior distribution of ", beta[j], " (j = 1, ..., p)")))
  abline(h = 0, lty = 2, col = 8)
  points(1:(p+2), beta.blm[,2], col = 4, pch = 19)
  arrows(1:(p+2), beta.blm[,1], 1:(p+2), beta.blm[,3], 
         angle = 90, code = 3, length = 0.03, col = 4)
  points(1:(p+2)+0.4, beta.mgp[,2], col = 2, pch = 19)
  arrows(1:(p+2)+0.4, beta.mgp[,1], 1:(p+2)+0.4, beta.mgp[,3], 
         angle = 90, code = 3, length = 0.03, col = 2)
}

# Let's have a look to the posterior correlation matrix of beta
corrplot::corrplot(cor(mgp.fit$trace$beta), method = "color")








