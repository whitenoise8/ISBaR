# file: fast_rmnorm.R
# author: Cristian Castiglione
# creation: 16/10/2023
# last change: 25/10/2023

#' @title Draw a sample from a multivariate Gaussian distribution
#' 
#' @description
#' \code{rmnorm} draws a sample \eqn{x} from a multivariate Gaussian distribution 
#' \eqn{N_p(\mu, \Sigma)} having mean vector \eqn{\mu = A^{-1} b} and variance matrix 
#' \eqn{\Sigma = A^{-1}}, where \eqn{A} is a \eqn{p \times p} positive definite 
#' matrix and \eqn{b} is a \eqn{p \times 1} vector.
#' 
#' @param A a square positive definite numeric matrix
#' @param b a numeric vector
#' 
#' @details
#' \code{rmnorm} uses a fast sampling scheme based on the algorithm proposed by Rue (2001).
#' 
#' 
#' @return A vector containing the sampled Gaussian variables.
#' 
#' @references 
#'   Rue (2001)
#'   Fast sampling of Gaussian Markov random fields.
#'   Journal of the Royal Statistical Society, Series B, 63(2): 325-338
#' 
#' @keywords internal
rmnorm = function (A, b) {
  R = chol(A)
  z = rnorm(length(b))
  y = drop(forwardsolve(t(R), b, upper.tri = FALSE, transpose = FALSE))
  x = drop(backsolve(R, y + z, upper.tri = TRUE, transpose = FALSE))
  return (x)
}

rmnorm2 = function (A, b) {
  R = chol(A)
  z = rnorm(length(b))
  y = drop(solve(t(R), b))
  x = drop(solve(R, y + z))
  return (x)
}

#' @title Solve the linear system using the forward-backward Cholesky algorithm
#' 
#' @description
#' \code{lschol} implements the \code{LS-CHOL} algorithm to solve the linear system
#' \eqn{Ax = b} with respect to \eqn{x}, where \eqn{A} is a square positive definite
#' matrix and \eqn{b} is a vector or a matrix with compatible dimensions.
#' 
#' @param A a square positive definite numeric matrix
#' @param b a numeric vector/matrix
#' @param upper if \code{TRUE}, uses the upper Cholesky factor of \code{A}
#' 
#' @return A list containing:
#' \describe{
#'   \item{\code{sol}}{the solution of the linear system}
#'   \item{\code{chol}}{the upper/lower Cholesky factor of \code{A}}
#'   \item{\code{inv}}{the inverse Cholesky factor of \code{A}}
#' }
#' 
#' @keywords internal
lschol = function (A, b, upper = FALSE) {
  out = list()
  out = NULL
  if (upper) {
    R = chol(A)
    Q = backsolve(R, diag(length(b)), upper.tri = TRUE, transpose = TRUE)
    s = drop(crossprod(Q, Q %*% b))
    out = list(sol = s, chol = R, inv = t(Q))
  } else {
    L = t(chol(A))
    Q = forwardsolve(L, diag(length(b)), upper.tri = FALSE, transpose = FALSE)
    s = drop(crossprod(Q, Q %*% b))
    out = list(sol = s, chol = L, inv = Q)
  }
  return (out)
}


#' @title Solve the linear system using the forward-backward Cholesky algorithm
#' 
#' @description
#' \code{lschol2} is an alternative implementation of \code{lschol}.
#' 
#' @param A a square positive definite numeric matrix
#' @param b a numeric vector/matrix
#' @param upper if \code{TRUE}, uses the upper Cholesky factor of \code{A}
#' 
#' @return A list containing:
#' \describe{
#'   \item{\code{sol}}{the solution of the linear system}
#'   \item{\code{chol}}{the upper/lower Cholesky factor of \code{A}}
#'   \item{\code{inv}}{the inverse Cholesky factor of \code{A}}
#' }
#' 
#' @keywords internal
lschol2 = function (A, b, upper = FALSE) {
  out = list()
  out = NULL
  if (upper) {
    R = chol(A)
    Q = backsolve(R, diag(length(b)), upper.tri = TRUE, transpose = TRUE)
    s = backsolve(R, drop(Q %*% b), upper.tri = TRUE, transpose = FALSE)
    out = list(sol = s, chol = R, inv = t(Q))
  } else {
    L = t(chol(A))
    Q = forwardsolve(L, diag(length(b)), upper.tri = FALSE, transpose = FALSE)
    s = backsolve(L, drop(Q %*% b), upper.tri = FALSE, transpose = TRUE)
    out = list(sol = s, chol = L, inv = Q)
  }
  return (out)
}


#' @title Solve the linear system using the forward-backward Cholesky algorithm
#' 
#' @description
#' \code{lschol3} is equivalent to \code{lschol}, but instead of returning both the
#' solution of the linear system and the cholesky factor of \code{A}, it just returns
#' the solution, saving memory and execution time.
#' 
#' @param A a square positive definite numeric matrix
#' @param b a numeric vector/matrix
#' @param upper if \code{TRUE}, uses the upper Cholesky factor of \code{A}
#' 
#' @keywords internal
lschol3 = function (A, b, upper = TRUE) {
  s = NULL
  if (upper) {
    R = chol(A)
    r = forwardsolve(R, b, upper.tri = TRUE, transpose = TRUE)
    s = backsolve(R, r, upper.tri = TRUE, transpose = FALSE)
  } else {
    L = t(chol(A))
    r = forwardsolve(L, b, upper.tri = FALSE, transpose = FALSE)
    s = backsolve(L, r, upper.tri = FALSE, transpose = TRUE)
  }
  return (s)
}

TEST = FALSE
if (TEST) {
  
  n = 2000
  p = 250
  
  A = crossprod(matrix(rnorm(n * p), nrow = n, ncol = p))
  b = rnorm(p)
  sol  = solve(A, b)
  solR = lschol(A, b, TRUE)
  solL = lschol(A, b, FALSE)
  
  rmvnorm2 = function (A, b) {
    V = solve(A); m = V %*% b; drop(mvtnorm::rmvnorm(1, m, V))
  }
  
  print(all.equal(sol, solL$sol))
  print(all.equal(sol, solR$sol))
  print(all.equal(solL$sol, solR$sol))
  print(all.equal(solL$chol, t(solR$chol)))
  print(all.equal(solL$inv, t(solR$inv)))
  
  rbenchmark::benchmark(
    lschol(A, b, FALSE),
    lschol(A, b,  TRUE),
    lschol2(A, b, FALSE),
    lschol2(A, b,  TRUE),
    lschol3(A, b, FALSE),
    lschol3(A, b,  TRUE),
    columns = c("test", "replications", "elapsed", "relative"))
  
  rbenchmark::benchmark(
    rmnorm(A, b), rmnorm2(A, b), rmvnorm2(A, b),
    columns = c("test", "replications", "elapsed", "relative"))
}
