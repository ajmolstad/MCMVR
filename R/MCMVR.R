# --------------------------------------------------------
# Please contact amolstad@fredhutch.org with issues  
# --------------------------------------------------------
tr.op1 <- function(x){
  sum(x^2)
}

tr.op2 <- function(x,y){
  sum(y*x)
}

svd2 <- function(x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE)
{
  #print("LINPACK:"); print(LINPACK)  ## added so you can see it's changed
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite or missing values in 'x'")
  dx <- dim(x)
  n <- dx[1L]
  p <- dx[2L]
  if (!n || !p)
    stop("a dimension is zero")
  La.res <- La.svd(x, nu, nv)   ## your problem line
  res <- list(d = La.res$d)
  if (nu)
    res$u <- La.res$u
  if (nv) {
    if (is.complex(x))
      res$v <- Conj(t(La.res$vt))
    else res$v <- t(La.res$vt)
  }
  res
}


# ----------------------------------------------------
# inner iteration function for nuclear norm penalty 
# ---------------------------------------------------
Acc.Prox.NN <- function(x, y, xtx, xty, D, beta.old, tau, lambda, weight, inner.quiet = TRUE, epsilon.beta = 1e-8, max.iter.inner = 2e4){
  
  # ---------------------------------
  # Nuclear norm specific functions 
  # ---------------------------------
  prox <- function(input, tau){
    keep <- svd2(input)
    return(list("out" = tcrossprod(keep$u*(tcrossprod(rep(1, dim(keep$u)[1]),pmax(keep$d - tau, 0))),keep$v), "d" = pmax(keep$d - tau, 0)))
  }

  pen <- function(input, weight){
    return(sum(svd2(input)$d))
  }


  # -- set backtrack parameter and compute initial obj.func val
  Dinner <-  rep(1, dim(y)[2])%*%t(diag(D))
  gam <- 10
  rho <- 1e-8
  old <- MCnegloglikCpp(beta = beta.old, x = x, y = y, tau = tau, d = Dinner)
  loglik <- rep(0, max.iter.inner)
  loglik[1] <- old + lambda*pen(beta.old, weight)
  # -- orig is obj function at initial value 
  orig <- abs(loglik[1])
  
  # -- preliminaries for iteration 
  iterating <- TRUE
  beta.iter <- 1 
  alphakm1 <- 0
  alphak <- 1
  betak <- beta.old
  betakm1 <- beta.old
  tildebetak <- beta.old 
  
  while(iterating){
    
    # ------------------------------------------
    # extrapolation step -- compute gradient and obj.func value 
    # ------------------------------------------
    gammak <- betak + (alphakm1/alphak)*(tildebetak - betak) + ((alphakm1 - 1)/(alphak))*(betak - betakm1) 
    
    # ------------------------------------------
    # get BB step sizes 
    # ------------------------------------------
    if(beta.iter > 1){
      
      # ----------------
      sk <- tildebetak - gammakm1; rk <- MCgradbetaCpp(beta = tildebetak, x = x, y = y, xtx = xtx, xty = xty, tau = tau, d = Dinner) - temp.gamma
      tx <- pmax(abs(tr.op1(sk)/tr.op2(sk, rk)), abs(tr.op2(sk, rk)/tr.op1(rk)))
      sk <- thetak - betakm1; rk <- MCgradbetaCpp(beta = thetak, x = x, y = y, xtx = xtx, xty = xty, tau = tau, d = Dinner) - temp.beta
      ty <- pmax(abs(tr.op1(sk)/tr.op2(sk, rk)), abs(tr.op2(sk, rk)/tr.op1(rk)))
      
      if(is.na(ty) | is.na(tx) | ty == 0 | tx == 0 | ty == Inf | tx == Inf | ty == -Inf | tx == -Inf){
        tx <- 1; ty <- 1
      }
      
      } else {
      
       tx <- 1; ty <- 1
      
      }
    
    
    # ------------------------------------------
    # line search for tildebetakp1
    # ------------------------------------------
    temp.gamma <- MCgradbetaCpp(beta = gammak, x = x, y = y, xtx = xtx, xty = xty, tau = tau, d = Dinner)
    old.F <- MCnegloglikCpp(beta = gammak, x = x, y = y, tau = tau, d = Dinner) + lambda*sum(svd2(gammak)$d)
    updating <- TRUE
    
    # --- backtracking line search 
    while(updating){
      
      # -- new iterate proximal step
      upTemp <- prox(gammak - ty*temp.gamma, ty*lambda)
      tildebetakp1 <- upTemp$out
      up.F.tildebeta <- MCnegloglikCpp(beta = tildebetakp1, x = x, y = y, tau = tau, d = Dinner) + lambda*sum(upTemp$d)#pen(tildebetakp1, weight, gamma = alpha)
      
      if(up.F.tildebeta < old.F - rho*sum((tildebetakp1 - gammak)^2)){
        updating <- FALSE
      } else {
        ty <- ty/gam
      }
      
      if(ty < 1e-8){
        updating <- FALSE
      }
    }
    
    # ------------------------------------------
    # line search for thetakp1
    # ------------------------------------------
    temp.beta <- MCgradbetaCpp(beta = betak, x = x, y=y, xtx=xtx, xty = xty, tau = tau, d = Dinner)
    old.F <- MCnegloglikCpp(beta = betak, x = x, y = y, tau = tau, d = Dinner) + lambda*pen(betak, weight)
    updating <- TRUE
    
    # --- backtracking line search 
    while(updating){
      
      # -- new iterate proximal step
      upTemp <- prox(betak - tx*temp.beta, tx*lambda)
      thetakp1 <- upTemp$out
      up.F.theta <- MCnegloglikCpp(beta = thetakp1, x = x, y = y, tau = tau, d = Dinner) + lambda*sum(upTemp$d)#pen(thetakp1, weight, gamma = alpha)
      inner <- (up.F.theta < old.F - rho*sum((thetakp1 - betak)^2))

      if(inner){
        updating <- FALSE
      } else {
        tx <- tx/gam
      }
      
      if(tx < 1e-8){
        updating <- FALSE
      }
    }
    
    # --------------------------------------
    # update step size
    # --------------------------------------
    alphakp1 <- (sqrt(4*alphak^2 + 1) + 1)/2
    
    
    # --------------------------------------
    # update beta
    # --------------------------------------
    if(up.F.tildebeta <= up.F.theta){
        betakp1 <- tildebetakp1
       loglik[beta.iter+1] <-up.F.tildebeta
    } else {
      betakp1 <- thetakp1
      loglik[beta.iter+1] <- up.F.theta
    }
    
    # --------------------------------------
    # check convergence
    # --------------------------------------
    #loglik[beta.iter+1] <- MCnegloglikCpp(betakp1, x, y, tau, Dinner) + lambda*pen(betakp1, weight, gamma = alpha)

    if(beta.iter > 3){
      if(abs(loglik[beta.iter-1] - loglik[beta.iter+1])/orig < epsilon.beta){
          iterating <- FALSE
      }
    }
    
    if(!inner.quiet){
      cat(loglik[beta.iter+1], "\n")
    }
    
    
    # --------------------------------------
    # update variables 
    # --------------------------------------
    alphakm1 <- alphak
    alphak <- alphakp1
    tildebetak <- tildebetakp1
    thetak <- thetakp1
    betakm1 <- betak
    betak <- betakp1
    gammakm1 <- gammak
    beta.iter <- beta.iter + 1
    
    if(beta.iter > max.iter.inner){
      iterating = FALSE
    }
    
  }
  
  temp <- MCgradbetaCpp(beta = betak, x = x, y = y, xtx = xtx, xty = xty, tau = tau, d = Dinner)
  #check1 <- sum((abs(tau*temp) <= tau*alpha*weight*lambda)[which(betak == 0)])/sum(betak == 0)
  #check2 <- max(abs((temp + alpha*weight*lambda*sign(betak) - (1-alpha)*lambda*betak)[which(betak != 0)]))
  
  return(list("beta" = betak, "objfunc" = loglik[1:(beta.iter - 2)], "check1" = sum(svd2(betak)$d > 1e-8)))#, "check2" = check2))
}

# ----------------------------------------------------
# inner iteration function for lasso penalty 
# ---------------------------------------------------
Acc.Prox.L1 <- function(x, y, xtx, xty, D, beta.old, tau, lambda, weight, inner.quiet = TRUE, epsilon.beta = 1e-8, max.iter.inner = 2e4){
  
  # ---------------------------------
  # Nuclear norm specific functions 
  # ---------------------------------
  prox <- function(input, tau){
    return(list("out" = pmax(abs(input) - tau, 0)*sign(input)))
  }

  pen <- function(input, weight){
    return(sum(abs(input)))
  }


  # -- set backtrack parameter and compute initial obj.func val
  Dinner <-  rep(1, dim(y)[2])%*%t(diag(D))
  gam <- 10
  rho <- 1e-8
  old <- MCnegloglikCpp(beta = beta.old, x = x, y = y, tau = tau, d = Dinner)
  loglik <- rep(0, max.iter.inner)
  loglik[1] <- old + lambda*pen(beta.old, weight)
  # -- orig is obj function at initial value 
  orig <- abs(loglik[1])
  
  # -- preliminaries for iteration 
  iterating <- TRUE
  beta.iter <- 1 
  alphakm1 <- 0
  alphak <- 1
  betak <- beta.old
  betakm1 <- beta.old
  tildebetak <- beta.old 
  
  while(iterating){
    
    # ------------------------------------------
    # extrapolation step -- compute gradient and obj.func value 
    # ------------------------------------------
    gammak <- betak + (alphakm1/alphak)*(tildebetak - betak) + ((alphakm1 - 1)/(alphak))*(betak - betakm1) 
    
    # ------------------------------------------
    # get BB step sizes 
    # ------------------------------------------
    if(beta.iter > 1){
      
      # ----------------
      sk <- tildebetak - gammakm1; rk <- MCgradbetaCpp(beta= tildebetak, x = x, y=y, xtx= xtx, xty = xty, tau = tau, d = Dinner) - temp.gamma
      tx <- pmax(abs(tr.op1(sk)/tr.op2(sk, rk)), abs(tr.op2(sk, rk)/tr.op1(rk)))
      sk <- thetak - betakm1; rk <- MCgradbetaCpp(beta= thetak, x = x, y = y, xtx = xtx, xty = xty, tau = tau, d = Dinner) - temp.beta
      ty <- pmax(abs(tr.op1(sk)/tr.op2(sk, rk)), abs(tr.op2(sk, rk)/tr.op1(rk)))
      
      if(is.na(ty) | is.na(tx) | ty == 0 | tx == 0 | ty == Inf | tx == Inf | ty == -Inf | tx == -Inf){
        tx <- 1; ty <- 1
      }
      
      } else {
      
       tx <- 1; ty <- 1
      
      }
    
    
    # ------------------------------------------
    # line search for tildebetakp1
    # ------------------------------------------
    temp.gamma <- MCgradbetaCpp(beta = gammak, x = x, y = y, xtx = xtx, xty = xty, tau = tau, d = Dinner)
    old.F <- MCnegloglikCpp(beta = gammak, x = x, y = y, tau = tau, d = Dinner) + lambda*pen(gammak, weight)
    updating <- TRUE
    
    # --- backtracking line search 
    while(updating){
      
      # -- new iterate proximal step
      upTemp <- prox(gammak - ty*temp.gamma, ty*lambda)
      tildebetakp1 <- upTemp$out
      up.F.tildebeta <- MCnegloglikCpp(beta = tildebetakp1, x = x, y = y, tau = tau, d = Dinner) + lambda*pen(tildebetakp1, weight)#pen(tildebetakp1, weight, gamma = alpha)
      
      if(up.F.tildebeta < old.F - rho*sum((tildebetakp1 - gammak)^2)){
        updating <- FALSE
      } else {
        ty <- ty/gam
      }
      
      if(ty < 1e-10){
        updating <- FALSE
      }
    }
    
    
    # ------------------------------------------
    # line search for thetakp1
    # ------------------------------------------
    temp.beta <- MCgradbetaCpp(beta = betak, x = x, y = y, xtx = xtx, xty = xty, tau = tau, d = Dinner)
    old.F <- MCnegloglikCpp(beta = betak, x = x, y = y, tau = tau, d = Dinner) + lambda*pen(betak, weight)
    updating <- TRUE
    
    # --- backtracking line search 
    while(updating){
      
      # -- new iterate proximal step
      upTemp <- prox(betak - tx*temp.beta, tx*lambda)
      thetakp1 <- upTemp$out
      up.F.theta <- MCnegloglikCpp(beta = thetakp1, x = x, y = y, tau = tau, d = Dinner) + lambda*pen(thetakp1, weight)#pen(thetakp1, weight, gamma = alpha)
      inner <- (up.F.theta < old.F - rho*sum((thetakp1 - betak)^2))

      if(inner){
        updating <- FALSE
      } else {
        tx <- tx/gam
      }
      
      if(tx < 1e-10){
        updating <- FALSE
      }
    }
    
    # --------------------------------------
    # update step size
    # --------------------------------------
    alphakp1 <- (sqrt(4*alphak^2 + 1) + 1)/2
    
    
    # --------------------------------------
    # update beta
    # --------------------------------------
    if(up.F.tildebeta <= up.F.theta){
        betakp1 <- tildebetakp1
       loglik[beta.iter+1] <-up.F.tildebeta
    } else {
      betakp1 <- thetakp1
      loglik[beta.iter+1] <- up.F.theta
    }
    
    # --------------------------------------
    # check convergence
    # --------------------------------------
    #loglik[beta.iter+1] <- MCnegloglikCpp(betakp1, x, y, tau, Dinner) + lambda*pen(betakp1, weight, gamma = alpha)

    if(beta.iter > 3){
      if(abs(loglik[beta.iter-1] - loglik[beta.iter+1])/orig < epsilon.beta){
          iterating <- FALSE
      }
    }
    
    if(!inner.quiet){
      cat(loglik[beta.iter+1], "\n")
    }
    
    
    # --------------------------------------
    # update variables 
    # --------------------------------------
    alphakm1 <- alphak
    alphak <- alphakp1
    tildebetak <- tildebetakp1
    thetak <- thetakp1
    betakm1 <- betak
    betak <- betakp1
    gammakm1 <- gammak
    beta.iter <- beta.iter + 1
    
    if(beta.iter > max.iter.inner){
      iterating = FALSE
    }
    
  }
  
  temp <- MCgradbetaCpp(beta = betak, x = x, y = y, xtx = xtx, xty= xty, tau = tau, d = Dinner)
  check1 <- NULL#sum((abs(tau*temp) <= tau*alpha*weight*lambda)[which(betak == 0)])/sum(betak == 0)
  check2 <- NULL#max(abs((temp + alpha*weight*lambda*sign(betak) - (1-alpha)*lambda*betak)[which(betak != 0)]))
  
  return(list("beta" = betak, "objfunc" = loglik[1:(beta.iter - 2)], "check1" = check1, "check2" = check2))
}


# -- Cross-validation for tau and lambda grid (initialize at lasso solution)
MCMVR.cv <- function(X, Y, tau.vec, nlambda, lambda.vec = NULL, nfolds = NULL, 
  delta = .01, tol = 1e-8, quiet= TRUE, inner.quiet= TRUE, penalty="L1"){
    
    # ---- assign algorithm based on penalty
    if(penalty=="L1"){
      Acc.Prox <- Acc.Prox.L1
    } 
    if(penalty=="NN"){
      Acc.Prox <- Acc.Prox.NN
    }

    if(penalty!="L1" & penalty!="NN"){
      stop("penalty argument must be \"NN\" or \"L1\" ")
    }

    # ---- preliminaries
    n <- dim(X)[1]
    p <- dim(X)[2]
    q <- dim(Y)[2]
    
    # ----- standardize if necessary 
    x <- (X - rep(1, n)%*%t(apply(X, 2, mean)))
    xtx <- crossprod(x)
    y <- (Y - rep(1, n)%*%t(apply(Y, 2, mean)))
    xty <- crossprod(x, y)

    # --- D matrix: can be adjusted for generalized version -- would require some new args 
    D <- diag(1, p)
    
    # ---- select tuning parameters for lambda 
    if(is.null(lambda.vec)){
      lambda.vec <- rep(0, length=nlambda)
      if(penalty=="L1"){
        lambda.max <- 2*(dim(x)[1]^(-1))*max(abs(xty)) + 1e-6
      }
      if(penalty=="NN"){
        lambda.max <- 2*(dim(x)[1]^(-1))*max(svd2(xty)$d) + 1e-6
      }
      lambda.min <- delta*lambda.max
      for(kk in 1:nlambda){
        lambda.vec[kk] <- lambda.max^((nlambda-kk)/(nlambda-1))*lambda.min^((kk-1)/(nlambda-1))
      }
    }
    
    # --- for each lambda, initialize largest tau value at the lasso solution ------
    beta.full <- matrix(0, nrow = p*q*length(tau.vec), ncol = length(lambda.vec))
    sparsity.mat <- matrix(0, length(lambda.vec), length(tau.vec))
    beta.old <- matrix(0, nrow=p, ncol=q)
    beta.temp <- beta.old

    for(jj in 1:length(tau.vec)){
      beta.old <- matrix(0, nrow=p, ncol=q)
      for(kk in 1:length(lambda.vec)){

            temp <- Acc.Prox(x, y, xtx, xty, D,beta.old, tau.vec[jj], lambda.vec[kk]/tau.vec[jj], 
                                weight = weight, inner.quiet = inner.quiet, 
                                epsilon.beta = tol*tau.vec[jj], max.iter.inner = 1e5)
            beta.old <- temp$beta
            beta.full[(p*q*(jj-1)+1):(p*q*jj),kk] <- c(beta.old)
            sparsity.mat[kk,jj] <- sum(beta.old!=0)

          if(penalty=="NN"){
            if(!quiet){
              cat("tau =", tau.vec[jj],": lambda =",  lambda.vec[kk], ": nuclear norm =", sum(svd2(beta.old)$d) , "\n")
            }
          }
            
          if(penalty=="L1"){
            if(!quiet){
              cat("tau =", tau.vec[jj],": lambda =",  lambda.vec[kk], ": nonzero entries =", sum(beta.old!=0),"\n")
            }
          }
      }
    }
    
    if(!is.null(nfolds)){
      
      # --- perform cross validation
      fold <- sample(rep(1:nfolds, length=n))
      cv.index <- split(1:n, fold)
      errs_wpred <- array(Inf, dim=c(length(lambda.vec),length(tau.vec), nfolds))
      errs_pred <- array(Inf, dim=c(length(lambda.vec),length(tau.vec), nfolds))
      
      for(k in 1:nfolds){
        
        # --- center X and Y --------------------
        ntrain <- dim(X[-cv.index[[k]],])[1]
        x.inner <- (X[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(X[-cv.index[[k]],], 2, mean)))#/(rep(1, ntrain)%*%t(apply(X[-cv.index[[k]], ], 2, sd)))
        xtx.inner <- crossprod(x.inner)
        y.inner <- (Y[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(Y[-cv.index[[k]],], 2, mean)))#/(rep(1, ntrain)%*%t(apply(Y[-cv.index[[k]], ], 2, sd)))
        xty.inner <- crossprod(x.inner,  y.inner)
        
        # ------------------------------------------------
        n.test <- length(cv.index[[k]])
        x.test <- (X[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(X[-cv.index[[k]],], 2, mean)))#/(rep(1, n.test)%*%t(apply(X[-cv.index[[k]], ], 2, sd)))
        y.test <- (Y[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(Y[-cv.index[[k]],], 2, mean)))#/(rep(1, n.test)%*%t(apply(Y[-cv.index[[k]], ], 2, sd)))
        weight <- matrix(1, nrow=p, ncol=q)
        # if(weighted){
        #   weight <- rep(1, p)%*%t(1/apply(y.test, 2, sd))
        #   #weight <- weight/max(weight)
        # }
        # ------------------------------------------------
        for(jj in 1:length(tau.vec)){
          beta.old <- matrix(0, nrow=p, ncol=q)
          for(kk in 1:length(lambda.vec)){
            temp <- Acc.Prox(x.inner, y.inner, xtx.inner, xty.inner, D, 
              beta.old, tau.vec[jj], lambda.vec[kk]/tau.vec[jj], weight = weight, 
              inner.quiet = inner.quiet, epsilon.beta = tol*tau.vec[jj], max.iter.inner = 1e5)
             beta.temp <- temp$beta
             beta.mat <- beta.temp
             residuals <- (y.test - x.test%*%temp$beta)
             errs_pred[kk, jj, k] <- mean(residuals^2)
             errs_wpred[kk, jj, k] <- mean((residuals%*%diag(apply(y.test, 2, sd)^(-1)))^2)
             beta.old <- temp$beta
          }
        }
        if(!quiet){
          cat("Through CV fold", k, "\n")
        }
      }
      
    } else {
      errs_pred <- NULL
      errs_wpred <- NULL
    }

    
    if(!is.null(errs_pred)){
      inds <- which(apply(errs_pred, c(1,2), mean) == min(apply(errs_pred, c(1,2), mean)), arr.ind = TRUE)
      tau.min <- tau.vec[inds[1,2]]
      lam.min <- lambda.vec[inds[1,1]]
    } else {
      tau.min <- NULL
      lam.min <- NULL
    }
    
    fit <-  list("beta" = beta.full, 
                 "sparsity.mat" = sparsity.mat,
                 "err.pred" = errs_pred, 
                 "err.wpred" = errs_wpred,
                 "Y.offset" = apply(Y, 2, mean), 
                 "X.offset" = apply(X, 2, mean),
                 "lambda.vec" = lambda.vec, 
                 "tau.vec" = tau.vec,
                 "tau.min" = tau.min,
                 "lam.min" = lam.min)
    
    class(fit) <- "MCMVR"
    return(fit)
}



MCMVR.coef <- function(fit, lambda = NULL, tau = NULL){
  
  if(is.null(lambda) && is.null(tau)){
    tau <- fit$tau.min
    lambda <- fit$lam.min
    if(is.null(fit$tau.min)){
      stop('No tuning parameters selected by CV')
    }
  } 
  
  
  if(class(fit)!="MCMVR"){
    stop('fit needs to be of class MCMVR (obtained using MCMVR.cv)')
  }
  
  tau.ind <- which(fit$tau.vec == tau)
  lam.ind <- which(fit$lambda.vec == lambda)
  p <- length(fit$X.offset)
  q <- length(fit$Y.offset)
  
  beta.vec <- fit$beta[(p*q*(tau.ind-1)+1):(p*q*tau.ind),lam.ind]
  beta.mat <- matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)

  # --- get intercept
  B0 <- fit$Y.offset - crossprod(beta.mat, fit$X.offset)
  return(list("beta0" = B0, "beta" = beta.mat))
}



MCMVR.predict <- function(Xnew, fit, lambda = NULL, tau = NULL){
  
  if(is.null(lambda) && is.null(tau)){
    tau <- fit$tau.min
    lambda <- fit$lam.min
    if(is.null(fit$tau.min)){
      stop('No tuning parameters selected by CV')
    }
  } 
  
  
  if(class(fit)!="MCMVR"){
    stop('fit needs to be of class MCMVR (obtained using MCMVR.cv)')
  }
  
  tau.ind <- which(fit$tau.vec == tau)
  lam.ind <- which(fit$lambda.vec == lambda)
  p <- length(fit$X.offset)
  q <- length(fit$Y.offset)
  
  beta.vec <- fit$beta[(p*q*(tau.ind-1)+1):(p*q*tau.ind),lam.ind]
  beta.mat <- matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)

  # --- get intercept 
  B0 <- fit$Y.offset - crossprod(beta.mat, fit$X.offset)

  # ---- prediction 
  if(dim(Xnew)[1] > 1){
    preds <- rep(1, dim(Xnew)[1])%*%t(B0) + Xnew%*%beta.mat
  } else {
    preds <- B0 + crossprod(beta.mat, Xnew)
  }

  return(list("preds" = preds, "beta0" = B0, "beta" = beta.mat))

}



