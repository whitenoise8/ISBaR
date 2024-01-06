
#' file: basis_wavelet.R
#' author: Cristian Castiglione
#' creation: 01/12/2023
#' last change: 01/12/2023

#' @title Create a wavelet design matrix of a certain order (default = 6)
#' 
#' @description
#' \code{get.wavelet} description
#' 
#' @param x vector of evaluation points
#' @param rng x range
#' @param n number of levels
#' @param f filter number
#' @param res resolution 
#' 
#' @return \code{W} description
#' 
#' @export
get.wavelet <- function(
    x, rng = range(x), n = 6, f = 5, res = 16384
) {
  # Load required package:
  require(wavethresh)
  
  # If only intercept
  if (n == 1) {
    W <- matrix(1,length(x),1)
  }
  
  if (n > 1) {
    # Check that x within support limits:
    if (any(x < rng[1]) | any(x > rng[2]))
      stop("All abscissae should be within rng values.")
    
    # Ensure that the number of levels is `allowable'.
    if (!any(n == c(1:10)))
      stop("Number of levels should be between 2 and 10.")
    
    # Ensure the res value is `allowable'.
    if (!any(res == c(2^(10:20))))
      stop("Resolution value should be a power of 2, with the power between 10 and 20.")
    
    # Transform x to the unit interval and obtain variables
    # required for linear interpolation:
    xUnit <- (x - rng[1]) / (rng[2] - rng[1])
    xUres <- xUnit * res
    fXuRes <- floor(xUres)
    
    # Set filter and wavelet family  
    family <- "DaubExPhase"
    K <- 2^n - 1
    
    # Create a dummy wavelet transform object
    wdObj <- wd(rep(0,res), filter.number = f, family = "DaubExPhase")
    
    Z <- matrix(0, length(x), K)
    for (k in 1:K) {
      # Create wobj so that it contains the Kth basis
      # function of the Z matrix with `res' regularly 
      # spaced points:
      putCobj <- putC(wdObj, level = 0, v = 0)
      putCobj$D <- putCobj$D * 0
      putCobj$D[res - k] <- 1
      
      # Obtain kth column of Z via linear interpolation
      # of the wr(putCobj) grid values:
      wtVec <- xUres - fXuRes
      wvVec <- wr(putCobj)
      wvVec <- c(wvVec, rep(wvVec[length(wvVec)], 2))
      Z[,k] <- sqrt(res) * ((1 - wtVec) * wvVec[fXuRes + 1] + wtVec*wvVec[fXuRes + 2])
    }
    
    # Create column indices to impose "left-to-right" ordering
    # within the same level:
    newColInds <- 1
    for (ell in 1:(n - 1))
      newColInds <- c(newColInds, (2^(ell + 1)-1):(2^(ell)))
    
    Z <- Z[, newColInds]
    W <- cbind(rep(1, length(x)), Z)
    
  }
  
  # Build the selection matrix
  S = matrix(c(1, rep(0, 2^n-1)), nrow = 2^n, ncol = 1)
  if (n > 1) {
    S = cbind(S, c(0, rep(1, 3), rep(0, 2^n - 4)))
    if (n > 2) {
      for (h in 2:(n-1)) {
        S = cbind(S, c(rep(0, 2^h), rep(1, 2^h), rep(0, 2^n - 2^(h+1))))
      }
    }
  }
  
  
  list(W = W, S = S)
}



