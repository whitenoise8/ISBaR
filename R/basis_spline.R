# file: basis_spline.R
# author: Cristian Castiglione
# creation: 09/08/2023
# last change: 24/10/2023
# references: 
#   Wand, Ormerod (2008) 
#   On semiparametric regression with O'Sullivan penalized splines
#   Australian & New Zealand Journal of Statistics, 50(2): 179-198

#' @title Create a vector of basis knots
#' @keywords internal
get.bspline.knots = function (x, n, unif = FALSE) {
  knots = NULL
  if (unif) {
    # Build a vector of knots uniformly displaced over the interval
    knots = seq(from = min(x), to = max(x), length = n)
  } else {
    # Build an adaptive vector of knots concentrating
    # more information where we have more data
    knots = quantile(unique(x), seq(0, 1, length = n+2)[-c(1, n+2)])
  }
  return (knots)
}

#' @title Create the basis matrix for third order B-spline regression
#' @keywords internal
get.bspline.matrix = function (x, knots, a, b) {
  splines::bs(x, knots = knots, degree = 3, Boundary.knots = c(a, b), intercept = TRUE)
}

#' @title Create the penalty matrix for third order O'SUllivan spline regression
#' @keywords internal
get.bspline.penalty = function (knots, a, b) {
  
  # Build a vector of expanded knots for Simpson quadrature
  allknots = c(rep(a, 4), knots, rep(b, 4))
  K = length(knots); L = 3 * (K+8)
  xt = (rep(allknots, each = 3)[-c(1, L-1, L)] + rep(allknots, each = 3)[-c(1, 2, L)]) / 2
  wts = rep(diff(allknots), each = 3) * rep(c(1, 4, 1) / 6, K+7)
  
  # Compute the second derivative matrix for cubic B-spline
  # over the Simpson quadrature knots
  ddB = splines::spline.des(allknots, xt, derivs = rep(2, length(xt)), outer.ok = TRUE)$design
  
  # Form the weighted crossproduct to integrate out the time
  S = crossprod(ddB, wts * ddB)
  
  # Return the exact differential penalty matrix
  return (S)
}

#' @title Create othogonalized design and penalty matrices for O'Sullivan spline regression
#' @keywords internal
convert.to.ospline = function (x, B = NULL, P = NULL, check = TRUE) {
  
  # Number of basis
  n = ncol(P) - 4
  
  # Indices for separating the kernel and the null space of P
  idz = 1:(n+2)
  idx = (n+3):(n+4)
  
  # Extract the scaled eigenvectors of P
  eig = eigen(P)
  UZ = eig$vectors[,idz]
  UX = eig$vectors[,idx]
  LZ = sweep(UZ, 2, sqrt(eig$values[idz]), "/")
  # equivalent but faster than:
  # ... LZ = t(t(UZ) / sqrt(eig$values[idz])) or
  # ... LZ = UZ %*% diag(1/sqrt(eig$values[idz]))
  
  # Stability check
  if (check) {
    L = cbind(UX, LZ)
    S = t(crossprod(L,t(crossprod(L,P)))) # = t(L) %*% P %*% L
    if (sum(S^2) > 1.0001 * (n+2)) {
      warning("Numerical instabilities arising from spectral decomposition.")
    }
  }
  
  # Build the orthogonalized basis matrices and re-sort 
  # the columns of Z from low to high frequency basis
  X = cbind(rep(1,length(x)), x - mean(x)) # UX
  Z = (B %*% LZ)[,(n+2):1]
  d = sqrt(eig$values[idz])[(n+2):1]
  
  # Return the null basis X, the othogonalized basis Z 
  # and the non-null eigenvalues d
  list(d = d, X = X, Z = Z)
}

#' @title Create standard design and penalty matrices for B-spline regression
#' 
#' @description
#' \code{get.bspline} compute the basis and penalty matrices for a B-spline basis expansion of dimension \code{n}
#' 
#' @param x description
#' @param n number of knots to include in the basis expansion
#' @param a lower bound of x
#' @param b upper bound of x
#' @param unif if \code{TRUE} select the knots uniformly in \code{[a,b]} (default \code{FALSE})
#' 
#' @return description
#' \describe{
#'   \item{\code{k}}{vector of basis knots}
#'   \item{\code{B}}{B-spline basis matrix}
#'   \item{\code{P}}{B-spline penalty matrix}
#' }
#' 
#' @export
get.bspline = function (x, n = 20, a = 0, b = 1, unif = FALSE) {
  
  # Get the knots, the basis matrix and the penalty matrix
  k = get.bspline.knots(x, n, unif)
  B = get.bspline.matrix(x, k, a, b)
  P = get.bspline.penalty(k, a, b)
  
  # Return the knots, the basis matrix and the penalty matrix
  list(k = k, B = B, P = P)
}

#' @title Create the othogonalized design matrix for O'Sullivan spline regression
#' 
#' @description
#' \code{get.ospline} compute the basis matrices for an O'Sullivan spline basis expansion,
#' which consists of an orthogonalized version of B-spline representation where the basis
#' matrix is rotated and scaled using the eigenvectors and eigenvalues of the B-spline 
#' penalty matrix.
#' 
#' @param x description
#' @param n number of knots to include in the basis expansion
#' @param a lower bound of x
#' @param b upper bound of x
#' @param unif if \code{TRUE}, selects the knots uniformly in \code{[a,b]} (default \code{FALSE})
#' @param check if \code{TRUE}, perform a stability check on the orthogonalized basis
#' 
#' @return \code{get.ospline} returns a list containing the following elements:
#' \describe{
#'   \item{\code{d}}{vector of positive eigenvalues of the B-spline penalty matrix}
#'   \item{\code{X}}{basis matrix spanning the null-space of the B-spline penalty}
#'   \item{\code{Z}}{basis matrix spanning the non-linear part of the expansion}
#' }
#' 
#' @references 
#'  Wand, Ormerod (2008) 
#'  On semiparametric regression with O'Sullivan penalized splines
#'  Australian & New Zealand Journal of Statistics, 50(2): 179-198
#' 
#' @export
get.ospline = function (x, n = 20, a = 0, b = 1, unif = FALSE, check = TRUE) {
  
  # Get the knots, the basis matrix and the penalty matrix
  k = get.bspline.knots(x, n, unif)
  B = get.bspline.matrix(x, k, a, b)
  P = get.bspline.penalty(k, a, b)
  
  # Return the transformed basis mtrices
  convert.to.ospline(x, B, P, check)
}



TEST = FALSE
if (TEST) {
  
  # Simulate noisy observations from a smooth curve
  n = 500; p = 5; a = 0; b = 2
  x = sort(runif(n, min = a, max = b))
  f = 3 * cos(2 * pi * x) / (2 * pi * x + pi) + 1 + 0.5 * x
  e = rnorm(n, mean = 0, sd = 0.1)
  y = f + e
  
  plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)
  
  # Get the B- and O-spline matrices
  bsp = get.bspline(x, p, a, b)
  osp = get.ospline(x, p, a, b)
  
  B = bsp$B # B-spline basis matrix
  P = bsp$P # B-spline penalty matrix
  X = osp$X # Non-penalized O-spline basis 
  Z = osp$Z # Orthogonalized O-spline basis
  
  C = cbind(1, x-1, Z)
  D = diag(c(0, 0, rep(1, p+2)))
  
  # Compare the standard and orthogonalized B-spline basis functions
  par(mfrow = c(2,1), mar = rep(2,4))
  matplot(x, B, type = "o", lty = 1, pch = 20, cex = 0.8, xlab = "", ylab = "", main = "B-spline basis")
  matplot(x, Z, type = "o", lty = 1, pch = 20, cex = 0.8, xlab = "", ylab = "", main = "O-spline basis")
  par(mfrow = c(1,1), mar = rep(4,4))
  
  df = data.frame(
    time = rep(x, times = p+4),
    basis = as.factor(rep(1:(p+4), each = n)),
    bspline = c(B),
    ospline = c(cbind(.25 * X, Z))
  )

  library(ggplot2)  
  ggplt = ggplot(data = df, mapping = aes(x = time, color = basis)) + 
    theme_minimal() + theme(axis.title.x = element_blank())
  bplt = ggplt + geom_line(mapping = aes(y = bspline), lwd = .75) + ylab("B-spline basis")
  oplt = ggplt + geom_line(mapping = aes(y = ospline), lwd = .75) + ylab("O-spline basis")
  ggpubr::ggarrange(bplt, oplt, nrow = 2, ncol = 1, legend = FALSE)
  
  # Compute the sufficient statistics for B- and O-spline regression
  lambda = 0.01
  btb = crossprod(B, B); bty = crossprod(B, y)
  ctc = crossprod(C, C); cty = crossprod(C, y)
  
  # Compute the penalized spline coefficients
  beta.bsp = drop(solve(btb + lambda * P, bty))
  beta.osp = drop(solve(ctc + lambda * D, cty))
  
  # Check if the estimated curves are equal (as they should be, by construction)
  all.equal(drop(B %*% beta.bsp), drop(C %*% beta.osp))
  
  # Compare the estimated curves
  plot(x, y, col = 8); lines(x, f, col = 2, lwd = 2)
  lines(x, drop(B %*% beta.bsp), col = 3, lty = 2, lwd = 3)
  lines(x, drop(C %*% beta.osp), col = 4, lty = 3, lwd = 3)
  legend("bottomright", col = c(2, 3, 4, 7), lty = 1:4, lwd = 2,
         legend = c("TRUE", "B-spline", "O-spline"))
}
