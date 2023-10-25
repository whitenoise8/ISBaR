#' fourier_basis.R
#' author: Cristian Castiglione
#' creation: 10/08/2023
#' last change: 24/10/2023
#' reference:
#'    Ramsay, Silverman (2005) 
#'    Functional data analysis, Second edition
#'    Springer

# Create the basis matrix for third order B-spline regression 
get.fourier = function (x, n, a = NULL, b = NULL, f = NULL) {
  if (is.null(a)) a = min(x) # lower bound
  if (is.null(b)) b = max(x) # upper bound
  if (is.null(b)) f = 1 # frequency
  x = (0.5 * pi) * (x - a) / (b - a)
  # K = tcrossprod(x, 1:n)
  # B = matrix(NA, nrow = length(x), ncol = 2*n)
  # B[,seq(from = 1, to = 2*n, by = 2)] = cos(K)
  # B[,seq(from = 2, to = 2*n, by = 2)] = sin(K)
  # B = t(t(B) / rep((1:n)^2, each = 2))
  B = cos(tcrossprod(x, 1:n))
  B = sweep(B, 2, (1:n)^2, "/") # equivalent but faster than = t(t(B) / (1:n)^2)
  return (B)
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
  Z = get.fourier(x, p)
  C = cbind(1, x-1, Z)
  D = diag(c(0, 0, rep(1, p)))
  
  # Compare the standard and orthogonalized B-spline basis functions
  matplot(x, Z, type = "o", lty = 1, pch = 20, cex = 0.8, 
          xlab = "", ylab = "", main = "Fourier basis")
  
  # Compute the sufficient statistics for Fourier regression
  lambda = 0.0001
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
