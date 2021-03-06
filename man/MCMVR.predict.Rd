\name{MCMVR.predict}
\alias{MCMVR.predict}
\title{MCMVR.predict}
\description{A function to make predictions based on a fitted model from \code{MCMVR.cv}.}
\usage{
MCMVR.predict(Xnew, fit, lambda = NULL, tau = NULL)
}

\arguments{
\item{Xnew}{A new \eqn{n_{\rm new} \times p} matrix of predictors to make predictions.}
\item{fit}{A fitted model from \code{MCMVR.cv}.}
\item{lambda}{The value of \eqn{\lambda} used to estimate \eqn{\beta}. If \code{NULL} and \code{fit} was obtained using \code{nfolds} non-\code{NULL}, will use tuning parameters which minimized prediction error in cross-validation.  }
\item{tau}{The value of \eqn{\tau} used to estimate \eqn{\beta}. If \code{NULL} and \code{fit} was obtained using \code{nfolds} non-\code{NULL}, will use tuning parameters which minimized prediction error in cross-validation. }
.}

\value{
\item{\code{preds}:}{An \eqn{n_{\rm new} \times q} matrix of predicted values.}

\item{\code{beta0}:}{The estimated \eqn{q}-dimensional intercept vector.  }

\item{\code{beta}: }{The estimated \eqn{p \times q} matrix of regression coefficients. }
}

