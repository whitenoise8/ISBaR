#' file: basis_fourier.R
#' author: Cristian Castiglione
#' creation: 10/08/2023
#' last change: 24/10/2023
#' reference:
#'    Ramsay, Silverman (2005) 
#'    Functional data analysis, Second edition
#'    Springer


#' @title Create the Fourier design matrix for scatter smoothing
#' 
#' @description
#' \code{get.fourier} compute the basis matrix for a Fourier basis expansion of dimension \code{n}
#' 
#' @param x description
#' @param n number of knots to include in the basis expansion
#' @param a lower bound of x
#' @param b upper bound of x
#' @param f frequency value
#' 
#' @return \code{get.fourier} returns a list containing the following elements:
#' \describe{
#'   \item{\code{d}}{vector of positive eigenvalues of the Fourier penalty matrix}
#'   \item{\code{X}}{basis matrix spanning the null-space of the Fourier penalty}
#'   \item{\code{Z}}{basis matrix spanning the non-linear part of the expansion}
#' }
#' 
#' @references 
#'   Ramsay, Silverman (2005) 
#'   Functional data analysis, Second edition
#'   Springer
#' 
#' @export
get.fourier = function (x, n, a = NULL, b = NULL, f = NULL) {
  if (is.null(a)) a = min(x) # lower bound
  if (is.null(b)) b = max(x) # upper bound
  if (is.null(f)) f = 1 # frequency
  m = floor((n-2) / 2)
  t = (2 * pi * f) * (x - a) / (b - a)
  d = rep(1:m, each = 2)^2
  X = cbind(1, t)
  B = matrix(NA, nrow = length(t), ncol = 2*m)
  for(h in 1:m) {
    B[, 2*h-1] = cos(h * t) / m^2
    B[, 2*h  ] = sin(h * t) / m^2
  }
  list(X = X, Z = Z, d = d)
}


TEST = FALSE
if (TEST) {
  
  # Simulate noisy observations from a smooth curve
  n = 500; p = 100; a = 0; b = 2
  x = sort(runif(n, min = a, max = b))
  f = 3 * cos(2 * pi * x) / (2 * pi * x + pi) + 1 + 0.5 * x
  e = rnorm(n, mean = 0, sd = 0.1)
  y = f + e
  
  plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)
  
  # Compute the O-spline design and penalty matrices
  fbasis = get.fourier(x, p)
  C = cbind(fbasis$X, fbasis$Z)
  D = diag(c(0, 0, rep(1, p-2)))
  
  # Compare the standard and orthogonalized B-spline basis functions
  matplot(x, Z, type = "o", lty = 1, pch = 20, cex = 0.8, 
          xlab = "", ylab = "", main = "Fourier basis")
  
  # Compute the sufficient statistics for Fourier regression
  lambda = 0.05
  ctc = crossprod(C, C)
  cty = crossprod(C, y)
  
  # Compute the penalized Fourier coefficients
  beta.hat = drop(solve(ctc + lambda * D, cty))
  
  # Compare the estimated curves
  plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)
  lines(x, drop(C %*% beta.hat), col = 4, lty = 2, lwd = 3)
  legend("bottomright", col = c(2, 4), lty = 1:2, lwd = 3,
         legend = c("TRUE", "Fourier"))
}

