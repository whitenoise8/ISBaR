#' mgp_gam_test.R
#' author: Cristian Castiglione
#' creation: 16/10/2023
#' last change: 16/10/2023

rm(list = ls())
graphics.off()

## Packages ----
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgcv)

devtools::load_all()

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
x = runif(n, min = a, max = b)
y = runif(n, min = a, max = b)
fx = exp(-1.25*x**2) * ((x-0.25)**2 + 0.3)
fy = .5 * cos(4 * pi/exp(y-a)) / (y-a+1)
e = rnorm(n, mean = 0, sd = 0.05)
u = fx + fy + e

par(mfrow = c(2,2))
plot(x, fx, type = "p")
plot(y, fy, type = "p")
plot(x, u, type = "p")
plot(y, u, type = "p")
par(mfrow = c(1,1))

## Spline regression ----

# Compute the basis matrices
X  = cbind(1, x-(b-a)/2, y-(b-a)/2)
Zx = get.ospline(x, p, a, b)$Z
Zy = get.ospline(y, p, a, b)$Z
Z  = list(Zx, Zy)
C  = cbind(X, Zx, Zy)

# Fit a Bayesian linear model with standard prior using MCMC
blm.prior = list(v = 0.1, a = 0.1, b = 0.1)
blm.control = list(burn = 1000, niter = 5000, report = 500, verbose = TRUE, thin = 5)
blm.fit = blm.reg.fit(u, X, Z, blm.prior, blm.control)
blm.pred = tcrossprod(C, blm.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Fit a Bayesian linear models with Horseshoe prior using MCMC
hs.control = list(burn = 1000, niter = 5000, thin = 5)
hs.fit = hs.reg.fit(u, X, Z, prior = list(), control = hs.control)
hs.pred = tcrossprod(C, hs.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Fit a Bayesian linear model with MGP prior using MCMC
mgp.prior = list(a1 = 2.1, a2 = 3.1, v = 0.1, a = 0.1, b = 0.1)
mgp.control = list(burn = 1000, niter = 5000, tol = 1e-05, adaptation = TRUE, report = 1000, verbose = TRUE, thin = 5)
mgp.fit = mgp.reg.fit(u, X, Z, mgp.prior, mgp.control)
mgp.pred = tcrossprod(C, mgp.fit$trace$beta) %>% 
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()

# Fit a cubic regression spline using PIRLS + REML
mgcv.fit = mgcv::gam(u ~ s(x, bs = "cr", k = p, m = 2) + s(y, bs = "cr", k = p, m = 2), 
                     data = data.frame(x = x, y = y, u = u),
                     family = gaussian(), method = "REML")
mgcv.pred = mgcv::predict.gam(mgcv.fit, se.fit = FALSE)
mgcv.se = mgcv::predict.gam(mgcv.fit, se.fit = TRUE)$se.fit
mgcv.pred = cbind(mgcv.pred - 1.96 * mgcv.se, mgcv.pred, mgcv.pred + 1.96 * mgcv.se)

# Let's have a look to the posterior of the (invese) shrinkage factors
{
  idx1 = mgp.fit$idx[[2]]
  idx2 = mgp.fit$idx[[3]]
  sigma1 = 1 / sqrt(mgp.fit$trace$tau[,2] + mgp.fit$trace$phi[,idx1])
  sigma1 = t(apply(sigma1, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE))
  sigma2 = 1 / sqrt(mgp.fit$trace$tau[,3] + mgp.fit$trace$phi[,idx2])
  sigma2 = t(apply(sigma2, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE))
  par(mfrow = c(1,2))
  matplot(sigma1, type = "o", col = 2:6, lty = 1, pch = 20, 
          xlab = "j", ylab = expression(sigma[j1]), 
          main = expression(paste(
            "Posterior distribution of ", 
            sigma[j1], " = ", (tau[1] * phi[j1])^"-1/2")))
  abline(h = median(1 / sqrt(blm.fit$trace$tau[,2])), col = 1, lty = 2)
  matplot(sigma2, type = "o", col = 2:6, lty = 1, pch = 20, 
          xlab = "j", ylab = expression(sigma[j2]), 
          main = expression(paste(
            "Posterior distribution of ", 
            sigma[j2], " = ", (tau[2] * phi[j2])^"-1/2")))
  abline(h = median(1 / sqrt(blm.fit$trace$tau[,3])), col = 1, lty = 2)
  legend("topright", col = 2:6, lty = 1, pch = 20, 
         legend = c("5%", "25%", "50%", "75%", "95%"))
  par(mfrow = c(1,1))
}

# Let's have a look to the posterior of the regression coefficients
{
  idx1 = mgp.fit$idx[[2]]
  idx2 = mgp.fit$idx[[3]]
  
  n1 = length(idx1)
  n2 = length(idx2)
  
  prb = c(0.025, 0.5, 0.975)
  prb = c(0.1, 0.5, 0.9)
  
  beta.blm = t(apply(blm.fit$trace$beta, 2, quantile, probs = prb))
  beta.mgp = t(apply(mgp.fit$trace$beta, 2, quantile, probs = prb))
  beta.hs = t(apply(hs.fit$trace$beta, 2, quantile, probs = prb))
  
  par(mfrow = c(1,2))
  plot(0, 0, type = "n", xlim = c(1, n1), ylim = c(-1,+1) * 2.5,
       xlab = "j", ylab = expression(beta[j1]), 
       main = expression(paste(
         "Posterior credibility intervals of ", beta[j1], " (j = 1, ..., p)")))
  abline(h = 0, lty = 2, col = 8)
  # Bayesian linear model
  points(1:n1, beta.blm[idx1,2], col = 4, pch = 19)
  arrows(1:n1, beta.blm[idx1,1], 1:n1, beta.blm[idx1,3], 
         angle = 90, code = 3, length = 0.02, col = 4)
  # Multiplicative Gamma process
  points(1:n1+0.2, beta.mgp[idx1,2], col = 2, pch = 19)
  arrows(1:n1+0.2, beta.mgp[idx1,1], 1:n1+0.2, beta.mgp[idx1,3], 
         angle = 90, code = 3, length = 0.02, col = 2)
  # Horseshoe prior
  points(1:n1+0.4, beta.hs[idx1,2], col = 3, pch = 19)
  arrows(1:n1+0.4, beta.hs[idx1,1], 1:n1+0.4, beta.hs[idx1,3], 
         angle = 90, code = 3, length = 0.02, col = 3)
  
  plot(0, 0, type = "n", xlim = c(1, n2), ylim = c(-1,+1) * 9,
       xlab = "j", ylab = expression(beta[j2]), 
       main = expression(paste(
         "Posterior credibility intervals of ", beta[j2], " (j = 1, ..., p)")))
  abline(h = 0, lty = 2, col = 8)
  # Bayesian linear model
  points(1:n2, beta.blm[idx2,2], col = 4, pch = 19)
  arrows(1:n2, beta.blm[idx2,1], 1:n2, beta.blm[idx2,3], 
         angle = 90, code = 3, length = 0.02, col = 4)
  # Multiplicative Gamma process
  points(1:n2+0.2, beta.mgp[idx2,2], col = 2, pch = 19)
  arrows(1:n2+0.2, beta.mgp[idx2,1], 1:n2+0.2, beta.mgp[idx2,3], 
         angle = 90, code = 3, length = 0.02, col = 2)
  # Horseshoe prior
  points(1:n2+0.4, beta.hs[idx2,2], col = 3, pch = 19)
  arrows(1:n2+0.4, beta.hs[idx2,1], 1:n2+0.4, beta.hs[idx2,3], 
         angle = 90, code = 3, length = 0.02, col = 3)
  legend("bottomleft", col = c(4, 2, 3), lty = 1, lwd = 2,
         legend = c("BLM", "MGP", "HS"))
}

# Signal reconstruction error
rmse = function (x, y = 0) sqrt(mean((x - y)^2))

rmse.blm = tcrossprod(C, blm.fit$trace$beta) |> apply(2, function(x) rmse(x, fx + fy))
rmse.mgp = tcrossprod(C, mgp.fit$trace$beta) |> apply(2, function(x) rmse(x, fx + fy))
rmse.hs = tcrossprod(C, hs.fit$trace$beta) |> apply(2, function(x) rmse(x, fx + fy))
rmse.mgcv = rmse(mgcv.pred[,2], fx + fy)

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

