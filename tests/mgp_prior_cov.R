#' mgp_prior_cov.R
#' author: Cristian Castiglione
#' creation: 04/09/2023
#' last change: 04/09/2023

rm(list = ls())
graphics.off()

## Packages ----
library(ggplot2)
library(ggpubr)
library(dplyr)

source("spline_basis.R")
source("fourier_basis.R")
source("mgp_reg_utils.R")
source("blm_reg_fit.R")
source("mgp_reg_fit.R")

## Data simulation ----
n = 51; p = 20; a = 0; b = 2
x = sort(runif(n, min = a, max = b))
x = seq(from = a, to = b, length = n+2)[-c(1,n+2)]
f = 3 * cos(2 * pi * x) / (2 * pi * x + pi) + 1 + 0.5 * x
e = rnorm(n, mean = 0, sd = 0.1)
y = f + e

plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)

## Spline regression ----

# Compute the basis and penalty matrices
k = get.bspline.knots(x, p)
B = get.bspline.matrix(x, k, a, b)
P = get.bspline.penalty(k, a, b)
Z = get.ospline.matrix(x, B, P)$z
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
  plot(0, 0, type = "n", xlim = c(1, p+4), ylim = c(-15,+15),
       xlab = "j", ylab = expression(beta[j]), 
       main = expression(paste(
         "Posterior distribution of ", beta[j], " (j = 1, ..., p)")))
  abline(h = 0, lty = 2, col = 8)
  points(1:(p+4), beta.blm[,2], col = 4, pch = 19)
  arrows(1:(p+4), beta.blm[,1], 1:(p+4), beta.blm[,3], 
         angle = 90, code = 3, length = 0.03, col = 4)
  points(1:(p+4)+0.4, beta.mgp[,2], col = 2, pch = 19)
  arrows(1:(p+4)+0.4, beta.mgp[,1], 1:(p+4)+0.4, beta.mgp[,3], 
         angle = 90, code = 3, length = 0.03, col = 2)
}

# Let's have a look to the posterior correlation matrix of beta
corrplot::corrplot(cor(mgp.fit$trace$beta), method = "color")



## IMPLIED COVARIANCE FUNCTION ----

p = 30

k = get.bspline.knots(x, p)
B = get.bspline.matrix(x, k, a, b)
P = get.bspline.penalty(k, a, b)
Z = get.ospline.matrix(x, B, P)$z
X = cbind(1, x-1)
C = cbind(X, Z)

# Compute the average prior variances for u
phi = exp(-mgp.prior.sim(p+2, 100, init.mgp.prior(list()))$logphi) |> colMeans()
plot(phi, type = "o", pch = 19)

# Orthogonalize the spline basis via spectral decomposition
idx = (p+2):(p+4)
idz = 1:(p+2)

s = eigen(P)
Ux = s$vectors[, idx]
Uz = s$vectors[, idz]
dz = sqrt(s$values[idz])

# Compute the reconstructed and new covariance matrix
{
  A = tcrossprod(Uz %*% diag(dz), Uz %*% diag(dz))
  C = tcrossprod(Uz %*% diag(dz), Uz %*% diag(dz) %*% diag(phi[(p+2):1]))
  
  par(mfrow = c(1,2))
  # plot3D::image2D(P, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(P)))
  plot3D::image2D(A, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(A)))
  plot3D::image2D(C, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(C)))
  par(mfrow = c(1,1))
}

{
  A = tcrossprod(Uz %*% diag(1/dz), Uz %*% diag(1/dz))
  C = tcrossprod(Uz %*% diag(1/dz), Uz %*% diag(1/dz) %*% diag(phi[(p+2):1]))
  
  par(mfrow = c(1,2))
  # plot3D::image2D(P, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(P)))
  plot3D::image2D(A, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(A)))
  plot3D::image2D(C, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(C)))
  par(mfrow = c(1,1))
}

{
  Z = B %*% Uz %*% diag(1/dz)
  W = B %*% Uz %*% diag(1/dz) %*% diag(sqrt(phi[(p+2):1]))
  
  par(mfrow = c(1,2))
  matplot(x, Z, type = "o", pch = 19, lty = 1, ylim = c(-.23, +.23), ylab = "basis")
  matplot(x, W, type = "o", pch = 19, lty = 1, ylim = c(-.23, +.23), ylab = "basis")
  par(mfrow = c(1,1))
}


{
  BABt = B %*% A %*% t(B)
  BCBt = B %*% C %*% t(B)
  
  par(mfrow = c(1,2))
  plot3D::image2D(BABt, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(BABt)))
  plot3D::image2D(BCBt, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(BCBt)))
  par(mfrow = c(1,1))
}

{
  D = 1 / sqrt(diag(A))
  RA = diag(D) %*% A %*% diag(D)
  
  D = 1 / sqrt(diag(C))
  RC = diag(D) %*% C %*% diag(D)
  
  par(mfrow = c(1,2))
  plot3D::image2D(RA, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(RA)))
  plot3D::image2D(RC, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(RC)))
  par(mfrow = c(1,1))
}

{
  D = 1 / sqrt(diag(BABt))
  RBABt = diag(D) %*% BABt %*% diag(D)
  
  D = 1 / sqrt(diag(BCBt))
  RBCBt = diag(D) %*% BCBt %*% diag(D)
  
  par(mfrow = c(1,2))
  plot3D::image2D(RBABt, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(RBABt)))
  plot3D::image2D(RBCBt, col = hcl.colors(51, "RdBu", rev = TRUE), clim = c(-1,+1) * max(abs(RBCBt)))
  par(mfrow = c(1,1))
}

idx = 25
plot(x, RBABt[idx,], type = "l", col = 4)
lines(x, RBCBt[idx,], type = "l", col = 2)

idx = 12
plot(RA[idx,], type = "l", col = 4)
lines(RC[idx,], type = "l", col = 2)


