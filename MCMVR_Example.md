MCMVR Example
================
Aaron J. Molstad (<amolstad@ufl.edu>)
7/18/2019

In this document, we provide a short tutorial on how to use the R package. We download this package from GitHub.

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("ajmolstad/MCMVR")
library(MCMVR)
```

First, we generate data from the \`\`errors-in-variables'' data generating model described in Section 3.1 of the article.

``` r
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
```

We fit the model using the cross-validation function. There are a number of key arguments: the first is `tau.vec`, where a user must specify a vector of candidate tuning parameters *τ* over which to fit the model:
$$ \\arg \\min\_{\\beta} {\\rm tr}\\left\\{n^{-1}(Y - X\\beta) (\\beta'\\beta + \\tau I\_q)^{-1}(Y - X\\beta)' \\right\\} + \\frac{\\lambda}{\\tau}{\\rm Pen}(\\beta).$$
 Another important argument is `nfolds`. If set to NULL, cross-validation is not performed, the model is fit to the complete data without cross-validation. Finally, a user must also decide which penalty to use. The current options are `penalty="L1"` and `penalty="NN"` which set ${\\rm Pen}(\\beta) = \\sum\_{j,k}|\\beta\_{j,k}|$ and ${\\rm Pen}(\\beta) = \\sum\_{j=1}^{\\min (p,q)} \\varphi\_j(\\beta)$, respectively, where *φ*<sub>*j*</sub>(*β*) is the *j*th largest singular value of *β*. Note that one only needs to input the number of candidate *λ*, `nlambda`, and the ratio of max to min lambda `delta`.

``` r
# ---------------------------------------------------
# Perform 5-fold CV for a grid of tuning parameters 
# ---------------------------------------------------
tau.vec <- 10^seq(3, 0, by=-.5)
fit <- MCMVR.cv(X = X, Y = Y, tau.vec = tau.vec, nlambda = 50, nfolds = 5, 
  delta = .01, tol = 1e-8, quiet= TRUE, inner.quiet= TRUE, penalty="L1")
str(fit)
```

    ## List of 10
    ##  $ beta        : num [1:3500, 1:50] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ sparsity.mat: num [1:50, 1:7] 0 1 2 2 2 2 5 6 8 8 ...
    ##  $ err.pred    : num [1:50, 1:7, 1:5] 6.66 6.58 6.51 6.45 6.4 ...
    ##  $ err.wpred   : num [1:50, 1:7, 1:5] 0.969 0.964 0.959 0.955 0.951 ...
    ##  $ Y.offset    : num [1:10] 0.2821 -0.1201 0.0626 -0.1167 -0.2071 ...
    ##  $ X.offset    : num [1:50] -0.03591 0.00857 -0.18562 -0.15828 -0.01975 ...
    ##  $ lambda.vec  : num [1:50] 6.12 5.57 5.07 4.62 4.2 ...
    ##  $ tau.vec     : num [1:7] 1000 316.2 100 31.6 10 ...
    ##  $ tau.min     : num 3.16
    ##  $ lam.min     : num 0.44
    ##  - attr(*, "class")= chr "EIVMR"

``` r
# ---------------------------------------------------
# visualize CV error
# ---------------------------------------------------
library(lattice)
levelplot(apply(fit$err.pred, c(1,2), mean), col.regions=grey(100:0/100), xlab=expression(lambda), ylab=expression(tau), aspect="fill")
```

![](MCMVR_Example_files/figure-markdown_github/fit%20model-1.png)

``` r
# ---------------------------------------------------
# get coefs from model which minimized CV error
# ---------------------------------------------------
betas <- MCMVR.coef(fit)$beta

par(mfrow=c(1,2))
image(betas!=0, col=grey(100:0/100), main=paste("Nonzero entries of estimate"))
image(beta!=0, col=grey(100:0/100), main=paste("Nonzero entries of truth"))
```

![](MCMVR_Example_files/figure-markdown_github/fit%20model-2.png)

``` r
# -------------------------------------------------
# make predictions for Xnew
# ------------------------------------------------
preds <- MCMVR.predict(Xnew, fit)

par(mfrow=c(1,2))
plot(preds$preds[,1], Ynew[,1], pch=20, main="Response 1", ylab="Our prediction", xlab="True response")
abline(0,1, lty=2)
plot(preds$preds[,2], Ynew[,2], pch=20, main="Response 2", ylab="Our prediction", xlab="True response")
abline(0,1, lty=2)
```

![](MCMVR_Example_files/figure-markdown_github/fit%20model-3.png)

``` r
# ---------------------------------------------------
# get coefs from model with seperate tuning parameter
# ---------------------------------------------------
betas <- MCMVR.coef(fit, tau = fit$tau.vec[1], lambda=fit$lambda.vec[12])$beta

par(mfrow=c(1,2))
image(betas!=0, col=grey(100:0/100), main=paste("Nonzero entries of estimate"))
image(beta!=0, col=grey(100:0/100), main=paste("Nonzero entries of truth"))
```

![](MCMVR_Example_files/figure-markdown_github/fit%20model-4.png)

We now also fit the model for the nuclear norm penalized version. Note that it is sometimes useful to relax the convergence tolerance for this version. We also turn off the `quiet` option. Note that `inner.quiet` should only be used for determining an appropriate value of `tol`.

``` r
# ---------------------------------------------------
# Perform 5-fold CV for a grid of tuning parameters 
# ---------------------------------------------------
tau.vec <- 10^seq(3, 0, by=-.5)
fit <- MCMVR.cv(X = X, Y = Y, tau.vec = tau.vec, nlambda = 20, nfolds = 5, 
  delta = .01, tol = 1e-10, quiet = FALSE, inner.quiet= TRUE, penalty="NN")
```

    ## 1 1 ; Nuclear norm =  0 
    ## 1 2 ; Nuclear norm =  0.4969919 
    ## 1 3 ; Nuclear norm =  1.491432 
    ## 1 4 ; Nuclear norm =  2.895799 
    ## 1 5 ; Nuclear norm =  4.394054 
    ## 1 6 ; Nuclear norm =  6.037959 
    ## 1 7 ; Nuclear norm =  7.633972 
    ## 1 8 ; Nuclear norm =  9.085348 
    ## 1 9 ; Nuclear norm =  10.42884 
    ## 1 10 ; Nuclear norm =  11.61476 
    ## 1 11 ; Nuclear norm =  12.65818 
    ## 1 12 ; Nuclear norm =  13.56972 
    ## 1 13 ; Nuclear norm =  14.36509 
    ## 1 14 ; Nuclear norm =  15.04664 
    ## 1 15 ; Nuclear norm =  15.62766 
    ## 1 16 ; Nuclear norm =  16.11633 
    ## 1 17 ; Nuclear norm =  16.5221 
    ## 1 18 ; Nuclear norm =  16.85975 
    ## 1 19 ; Nuclear norm =  17.13038 
    ## 1 20 ; Nuclear norm =  17.35113 
    ## 2 1 ; Nuclear norm =  0 
    ## 2 2 ; Nuclear norm =  0.5005064 
    ## 2 3 ; Nuclear norm =  1.498667 
    ## 2 4 ; Nuclear norm =  2.905169 
    ## 2 5 ; Nuclear norm =  4.403946 
    ## 2 6 ; Nuclear norm =  6.048851 
    ## 2 7 ; Nuclear norm =  7.645307 
    ## 2 8 ; Nuclear norm =  9.097865 
    ## 2 9 ; Nuclear norm =  10.44415 
    ## 2 10 ; Nuclear norm =  11.6329 
    ## 2 11 ; Nuclear norm =  12.67971 
    ## 2 12 ; Nuclear norm =  13.59708 
    ## 2 13 ; Nuclear norm =  14.39732 
    ## 2 14 ; Nuclear norm =  15.08625 
    ## 2 15 ; Nuclear norm =  15.67322 
    ## 2 16 ; Nuclear norm =  16.16929 
    ## 2 17 ; Nuclear norm =  16.58112 
    ## 2 18 ; Nuclear norm =  16.92016 
    ## 2 19 ; Nuclear norm =  17.19877 
    ## 2 20 ; Nuclear norm =  17.42263 
    ## 3 1 ; Nuclear norm =  0 
    ## 3 2 ; Nuclear norm =  0.5115146 
    ## 3 3 ; Nuclear norm =  1.520727 
    ## 3 4 ; Nuclear norm =  2.933977 
    ## 3 5 ; Nuclear norm =  4.434609 
    ## 3 6 ; Nuclear norm =  6.080995 
    ## 3 7 ; Nuclear norm =  7.679731 
    ## 3 8 ; Nuclear norm =  9.135377 
    ## 3 9 ; Nuclear norm =  10.48816 
    ## 3 10 ; Nuclear norm =  11.68609 
    ## 3 11 ; Nuclear norm =  12.74612 
    ## 3 12 ; Nuclear norm =  13.68037 
    ## 3 13 ; Nuclear norm =  14.498 
    ## 3 14 ; Nuclear norm =  15.20645 
    ## 3 15 ; Nuclear norm =  15.81302 
    ## 3 16 ; Nuclear norm =  16.3281 
    ## 3 17 ; Nuclear norm =  16.75909 
    ## 3 18 ; Nuclear norm =  17.11421 
    ## 3 19 ; Nuclear norm =  17.40677 
    ## 3 20 ; Nuclear norm =  17.64404 
    ## 4 1 ; Nuclear norm =  0 
    ## 4 2 ; Nuclear norm =  0.5472109 
    ## 4 3 ; Nuclear norm =  1.588968 
    ## 4 4 ; Nuclear norm =  3.02009 
    ## 4 5 ; Nuclear norm =  4.523484 
    ## 4 6 ; Nuclear norm =  6.172042 
    ## 4 7 ; Nuclear norm =  7.773269 
    ## 4 8 ; Nuclear norm =  9.235778 
    ## 4 9 ; Nuclear norm =  10.6041 
    ## 4 10 ; Nuclear norm =  11.82657 
    ## 4 11 ; Nuclear norm =  12.92126 
    ## 4 12 ; Nuclear norm =  13.90005 
    ## 4 13 ; Nuclear norm =  14.77008 
    ## 4 14 ; Nuclear norm =  15.53716 
    ## 4 15 ; Nuclear norm =  16.20631 
    ## 4 16 ; Nuclear norm =  16.7826 
    ## 4 17 ; Nuclear norm =  17.2738 
    ## 4 18 ; Nuclear norm =  17.68696 
    ## 4 19 ; Nuclear norm =  18.03076 
    ## 4 20 ; Nuclear norm =  18.31303 
    ## 5 1 ; Nuclear norm =  0 
    ## 5 2 ; Nuclear norm =  0.6627698 
    ## 5 3 ; Nuclear norm =  1.788089 
    ## 5 4 ; Nuclear norm =  3.248576 
    ## 5 5 ; Nuclear norm =  4.740305 
    ## 5 6 ; Nuclear norm =  6.378854 
    ## 5 7 ; Nuclear norm =  7.967467 
    ## 5 8 ; Nuclear norm =  9.424426 
    ## 5 9 ; Nuclear norm =  10.80619 
    ## 5 10 ; Nuclear norm =  12.06395 
    ## 5 11 ; Nuclear norm =  13.21911 
    ## 5 12 ; Nuclear norm =  14.2844 
    ## 5 13 ; Nuclear norm =  15.26726 
    ## 5 14 ; Nuclear norm =  16.17295 
    ## 5 15 ; Nuclear norm =  17.00502 
    ## 5 16 ; Nuclear norm =  17.76579 
    ## 5 17 ; Nuclear norm =  18.45903 
    ## 5 18 ; Nuclear norm =  19.08842 
    ## 5 19 ; Nuclear norm =  19.65529 
    ## 5 20 ; Nuclear norm =  20.16657 
    ## 6 1 ; Nuclear norm =  0 
    ## 6 2 ; Nuclear norm =  1.038415 
    ## 6 3 ; Nuclear norm =  2.249844 
    ## 6 4 ; Nuclear norm =  3.677558 
    ## 6 5 ; Nuclear norm =  5.069798 
    ## 6 6 ; Nuclear norm =  6.647228 
    ## 6 7 ; Nuclear norm =  8.119999 
    ## 6 8 ; Nuclear norm =  9.474282 
    ## 6 9 ; Nuclear norm =  10.78757 
    ## 6 10 ; Nuclear norm =  12.01057 
    ## 6 11 ; Nuclear norm =  13.17215 
    ## 6 12 ; Nuclear norm =  14.29402 
    ## 6 13 ; Nuclear norm =  15.3968 
    ## 6 14 ; Nuclear norm =  16.50087 
    ## 6 15 ; Nuclear norm =  17.61678 
    ## 6 16 ; Nuclear norm =  18.73428 
    ## 6 17 ; Nuclear norm =  19.83772 
    ## 6 18 ; Nuclear norm =  20.9153 
    ## 6 19 ; Nuclear norm =  21.96227 
    ## 6 20 ; Nuclear norm =  22.97429 
    ## 7 1 ; Nuclear norm =  0 
    ## 7 2 ; Nuclear norm =  1.401798 
    ## 7 3 ; Nuclear norm =  2.744489 
    ## 7 4 ; Nuclear norm =  3.88227 
    ## 7 5 ; Nuclear norm =  5.222365 
    ## 7 6 ; Nuclear norm =  6.613793 
    ## 7 7 ; Nuclear norm =  7.758001 
    ## 7 8 ; Nuclear norm =  8.889974 
    ## 7 9 ; Nuclear norm =  10.02008 
    ## 7 10 ; Nuclear norm =  11.08257 
    ## 7 11 ; Nuclear norm =  12.118 
    ## 7 12 ; Nuclear norm =  13.14494 
    ## 7 13 ; Nuclear norm =  14.17462 
    ## 7 14 ; Nuclear norm =  15.21416 
    ## 7 15 ; Nuclear norm =  16.26758 
    ## 7 16 ; Nuclear norm =  17.33527 
    ## 7 17 ; Nuclear norm =  18.41499 
    ## 7 18 ; Nuclear norm =  19.50372 
    ## 7 19 ; Nuclear norm =  20.59899 
    ## 7 20 ; Nuclear norm =  21.69823 
    ## Through CV fold 1 
    ## Through CV fold 2 
    ## Through CV fold 3 
    ## Through CV fold 4 
    ## Through CV fold 5

``` r
str(fit)
```

    ## List of 10
    ##  $ beta        : num [1:3500, 1:20] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ sparsity.mat: num [1:20, 1:7] 0 500 500 500 500 500 500 500 500 500 ...
    ##  $ err.pred    : num [1:20, 1:7, 1:5] 6.61 6.4 6.07 5.58 5.11 ...
    ##  $ err.wpred   : num [1:20, 1:7, 1:5] 0.979 0.954 0.906 0.836 0.775 ...
    ##  $ Y.offset    : num [1:10] 0.2821 -0.1201 0.0626 -0.1167 -0.2071 ...
    ##  $ X.offset    : num [1:50] -0.03591 0.00857 -0.18562 -0.15828 -0.01975 ...
    ##  $ lambda.vec  : num [1:20] 9.99 7.84 6.15 4.83 3.79 ...
    ##  $ tau.vec     : num [1:7] 1000 316.2 100 31.6 10 ...
    ##  $ tau.min     : num 1000
    ##  $ lam.min     : num 1.44
    ##  - attr(*, "class")= chr "EIVMR"

``` r
levelplot(apply(fit$err.pred, c(1,2), mean), col.regions=grey(100:0/100), xlab=expression(lambda), ylab=expression(tau), aspect="fill")
```

![](MCMVR_Example_files/figure-markdown_github/generate%20RR%20data-1.png)
