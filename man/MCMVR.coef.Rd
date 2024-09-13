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


\examples{
set.seed(1)
p <- 50
q <- 10
n <- 100

beta <- matrix(rnorm(p*q)*sample(c(0,1), p*q, prob = c(.9, .1), replace=TRUE), nrow=p, ncol=q)

Z <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- tcrossprod(Z, t(beta)) + matrix(rnorm(n*q, sd=1), nrow=n)
X <- Z + matrix(rnorm(n*p, sd=sqrt(0.5)), nrow=n)

Znew <- matrix(rnorm(n*p), nrow=n, ncol=p)
Xnew <- Znew + matrix(rnorm(n*p, sd=sqrt(0.5)), nrow=n)
Ynew <- tcrossprod(Znew, t(beta)) + matrix(rnorm(n*q, sd=1), nrow=n) 

# ---------------------------------------------------  
# Perform 5-fold CV for a grid of tuning parameters 
# ---------------------------------------------------
tau.vec <- 10^seq(3, 0, by=-.5)
fit <- MCMVR.cv(X = X, Y = Y, tau.vec = tau.vec, nlambda = 50, nfolds = 5, 
  delta = .01, tol = 1e-8, quiet= TRUE, inner.quiet= TRUE, penalty="L1")

# ------------------------------------------------------
# Extract coefficients and make predictions
# ------------------------------------------------------
betas <- MCMVR.coef(fit)$beta
preds <- MCMVR.predict(Xnew, fit)
}
