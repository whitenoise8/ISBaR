% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_rmnorm.R
\name{lschol3}
\alias{lschol3}
\title{Solve the linear system using the forward-backward Cholesky algorithm}
\usage{
lschol3(A, b, upper = TRUE)
}
\arguments{
\item{A}{a square positive definite numeric matrix}

\item{b}{a numeric vector/matrix}

\item{upper}{if \code{TRUE}, uses the upper Cholesky factor of \code{A}}
}
\description{
\code{lschol3} is equivalent to \code{lschol}, but instead of returning both the
solution of the linear system and the cholesky factor of \code{A}, it just returns
the solution, saving memory and execution time.
}
\keyword{internal}
