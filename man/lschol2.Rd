% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_rmnorm.R
\name{lschol2}
\alias{lschol2}
\title{Solve the linear system using the forward-backward Cholesky algorithm}
\usage{
lschol2(A, b, upper = FALSE)
}
\arguments{
\item{A}{a square positive definite numeric matrix}

\item{b}{a numeric vector/matrix}

\item{upper}{if \code{TRUE}, uses the upper Cholesky factor of \code{A}}
}
\value{
A list containing:
\describe{
  \item{\code{sol}}{the solution of the linear system}
  \item{\code{chol}}{the upper/lower Cholesky factor of \code{A}}
  \item{\code{inv}}{the inverse Cholesky factor of \code{A}}
}
}
\description{
\code{lschol2} is an alternative implementation of \code{lschol}.
}
\keyword{internal}
