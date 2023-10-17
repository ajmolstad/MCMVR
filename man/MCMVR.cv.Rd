\name{MCMVR.cv}
\alias{MCMVR.cv}
\title{MCMVR.cv}
\description{A function to perform cross-validation and model fitting. }\usage{
MCMVR.cv(X, Y, tau.vec, nlambda, lambda.vec = NULL, nfolds = NULL, 
  delta = .01, tol = 1e-8, quiet = TRUE, inner.quiet= TRUE, penalty="L1")
}

\arguments{
\item{X}{An \eqn{n \times p} matrix of predictors.}
\item{Y}{An \eqn{n \times q} matrix of response variables.}
\item{tau.vec}{A vector of candidate tuning parameters \eqn{\tau}. }
\item{nlambda}{The number of candidate tuning parameters \eqn{\lambda} to be used for model fitting and cross-validation. }
\item{lambda.vec}{A user specified set of candidate tuning parameters: should not be utilized for general purpose model fitting. }
\item{nfolds}{When \code{NULL}, does not perform cross-validation but fits the model for all combinations of \eqn{\tau} and \eqn{\lambda}}.
\item{delta}{Positive scalar less than one which indicates the ratio of smallest to large candidate \eqn{\lambda}. Default is .01. }.
\item{tol}{Convergence tolerance for inner accelerate proximal gradient descent algorithm.}
\item{quiet}{If \code{TRUE}, prints progress. }.
\item{inner.quiet}{If \code{TRUE}, prints objective function values for inner functions. Should be rarely used unless to help determine an appropriate \code{tol}. }
\item{penalty}{A character, either \code{"L1"} or \code{"NN"} for the \eqn{\ell_1}-penalized version of our method or the nuclear norm penalty version, respectively.}}

\value{
A fitted model and cross-validation diagnostics. To extract regression coefficients for certain values of \eqn{\tau} and \eqn{\lambda}, use \code{MCMVR.coef}. To perform prediction, use \code{MCMVR.pred}. 

\item{\code{beta}}{A matrix of dimension \eqn{pq}\code{length(tau.vec)}\eqn{ \times }\code{nlambda}. To extract coefficients for a specific pair of tuning parameters, use \code{MCMVR.coef}.}

\item{sparsity.mat}{A matrix of dimension \code{nlambda} \eqn{ \times } \code{length(tau.vec)} which provides the number of nonzero elements in \eqn{\beta} for each tuning parameter pair or the nuclear norm of the regression coefficient matrix. }

\item{\code{err.pred}}{An array of dimension \code{nlambda} \eqn{ \times } \code{length(tau.vec)} \eqn{ \times } \code{nfold} giving the mean cross-validation error for each fold. }

\item{\code{err.wpred}}{An array of dimension \code{nlambda} \eqn{ \times } \code{length(tau.vec)} \eqn{ \times } \code{nfold} giving the mean cross-validation weighted error for each fold (i.e., responses transformed to same scale). }

\item{\code{Y.offset}}{A \eqn{q}-dimensional vector giving the sample mean for the response in the training data.}

\item{\code{Z.offset}}{A \eqn{p}-dimensional vector giving the sample mean for the predictor in the training data.}

\item{\code{lambda.vec}}{A \code{nlambda}-dimensional vector giving the values of \eqn{\lambda} at which \eqn{\beta} were computed.}

\item{\code{tau.vec} }{The same as the input vector \code{tau.vec}.}

\item{\code{tau.min}}{The value of the tuning parameter \eqn{\tau} which minimizes prediction error in cross-validation. }

\item{\code{lam.min}}{The value of the tuning parameter \eqn{\lambda} which minimizes prediction error in cross-validation.}
}


