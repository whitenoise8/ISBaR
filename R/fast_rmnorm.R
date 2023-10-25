#' fast_rmnorm.R
#' author: Cristian Castiglione
#' creation: 16/10/2023
#' last change: 16/10/2023

# Sample from a Gaussian distribution N(m, V) with m = inv(A)*b, V = inv(A)
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

# Solve the linear system Ax = b using the forward-backward algorithm
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
