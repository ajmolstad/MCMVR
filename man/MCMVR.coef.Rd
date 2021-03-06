\name{MCMVR.coef}
\alias{MCMVR.coef}
\title{MCMVR.coef}
\description{A function to extract the estimated regression coefficient matrix and intercept from an fitted model from \code{MCMVR.cv}.}
\usage{
MCMVR.coef(fit, lambda = NULL, tau = NULL)
}

\arguments{
\item{fit}{A fitted model from \code{MCMVR.cv}.}
\item{lambda}{The value of \eqn{\lambda} used to estimate \eqn{\beta}. If \code{NULL} and \code{fit} was obtained using \code{nfolds} non-\code{NULL}, will use tuning parameters which minimized prediction error in cross-validation.  }
\item{tau}{The value of \eqn{\tau} used to estimate \eqn{\beta}. If \code{NULL} and \code{fit} was obtained using \code{nfolds} non-\code{NULL}, will use tuning parameters which minimized prediction error in cross-validation. }
}

\value{

\item{\code{beta0}:}{The estimated \eqn{q}-dimensional intercept vector.  }
\item{\code{beta}: }{The estimated \eqn{p \times q} matrix of regression coefficients. }
}

