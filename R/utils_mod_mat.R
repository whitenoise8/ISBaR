#' file: utils_mod_mat.R
#' author: Cristian Castiglione
#' creation: 05/01/2024
#' last change: 05/01/2024


#' Get the column-dimension of each block
#' @keywords internal
get.block.dim = function (X, Z) {
  c(ncol(X), unlist(lapply(Z, ncol)))
}

#' Get the column-indices of each block
#' @keywords internal
get.block.idx = function (p) {
  n = length(p)
  start = cumsum(c(1, p[-n]))
  end = cumsum(p)
  idx = apply(cbind(start, end), 1, function(x) x[1]:x[2])
  names(idx) = c("X", paste("Z", 1:(n-1), sep = ""))
  return (idx)
}

#' Get the column-indices of each block concatenated by row
#' @keywords internal
get.flat.idx = function (idx, q) {
  keep = c()
  for (h in 1:length(q)) {
    keep = c(keep, idx[[h]][1:q[h]])
  }
  return (keep)
}

#' Build the model matrix and the sufficient statistics
#' @keywords internal
get.matrix = function (X, Z) {
  C = cbind(X, do.call(cbind, Z))
  return (C)
}

#' Get the vector of penalty parameters
#' @keywords internal
get.penalty = function (tau, phi, idx) {
  lambda = rep(0, length = length(phi))
  for (h in 2:length(idx)) {
    ih = idx[[h]]
    lambda[ih] = tau[h] * phi[ih]
  }
  return (lambda)
}

#' Print simulation status
#' @keywords internal
print.status = function (iter, maxiter, keep, logpost, report, burn) {
  if (iter %% report == 0) {
    cat(gettextf(" iter: %5d /%5d  ", iter, maxiter),
        gettextf(" ncomp: %2d  ", length(keep)),
        gettextf(" logp: %.4f \n", logpost))
  }
  if (iter == burn) {
    cat(c(rep("-", 53), "\n"), sep = "")
  }
}
