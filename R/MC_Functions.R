# -----------------------------------------------------
# Complete set of functions to fit MC estimator
# Please contact amolstad@fredhutch.org with issues  
# --------------------------------------------------------

# - Proximal operator (can be modified by user)
	prox_L1  <- function(input,tau){return(pmax(abs(input) - tau, 0)*sign(input))}

# - inner iteration function 
	Acc.Prox.L1 <- function(z, y, ztz, zty, beta.old, tau, lambda, L = NULL, weight, inner.quiet = TRUE, epsilon.beta = 1e-8, max.iter.inner = 2e4){
	    
	    # -- set step size if NULL
		L0 <- 1/tau
		L <- 1/tau

		# -- set backtrack parameter and compute initial obj.func val
	    gam <- 2
	    old <- MCnegloglikCpp(beta.old, z, y, tau)
	    oldpen <- old + lambda*sum(abs(weight*beta.old))

	    # -- orig is obj function at initial value 
	    orig <- abs(old)
	    temp <- MCgradbetaCpp(beta.old, z, y, ztz, zty, tau)

	    # -- preliminaries for iteration 
	    iterating <- TRUE
	    beta.iter <- 1 
	    alphakm1 <- 1
	    alphak <- 1
	    betak <- beta.old
	    betakm1 <- beta.old
	    loglik <- rep(0, max.iter.inner)

	    # -------------------------------------------------------------
	    while(iterating){
	        
	        # -- extrapolation step -- compute gradient and obj.func value 
	        gammak <- betak + ((alphakm1 - 1)/(alphak))*(betak - betakm1)
	        temp <- MCgradbetaCpp(gammak, z, y, ztz, zty, tau)
	        old.F <- MCnegloglikCpp(gammak, z, y, tau)
	        updating <- TRUE
		        
	        # --- backtracking line search 
	        while(updating){

	        	# -- new iterate proximal step
	            beta.new <- prox_L1(gammak - L^(-1)*temp, weight*L^(-1)*lambda)
	            up.F <- MCnegloglikCpp(beta.new, z, y, tau)
	            
	            if(up.F < old.F - sum(diag(crossprod(temp, gammak - beta.new))) + (L/2)*sum((gammak - beta.new)^2)){
	                updating <- FALSE
	            } else {
	                L <- L*gam
	            }

	            if(L > 1e10){
	                updating <- FALSE
	            }
	        }

	        # --- check if new iterate decreasing obj.func value 
	        negloglik.temp <- MCnegloglikCpp(betak, z, y, tau)
	        if(up.F + lambda*sum(abs(weight*beta.new)) <   negloglik.temp + lambda*sum(abs(weight*betak))){
	            betakp1 <- beta.new
	            negloglik.temp <- up.F
	        } else {
	            betakp1 <- betak
	        }
	        
	        # --- adjust extrapolation parameter, iteration counter
	        alphakp1 <- (1 + (1 + 4*alphak^2)^.5)/2
	        beta.iter <- beta.iter + 1
	        
	        if(beta.iter > max.iter.inner){
	            iterating <- FALSE
	        }
	        
	        # -------------------------
	        betakm1 <- betak
	        betak <- betakp1
	        alphakm1 <- alphak
	        alphak <- alphakp1
	       	L <- L0

	        loglik[beta.iter-1] <- negloglik.temp + lambda*sum(abs(weight*betak))

	        if(beta.iter%%5000 == 0){
	            alphakm1 <- 1
	            alphak <- 1
	        }
	        if(!inner.quiet){
		        if(beta.iter%%500 == 0){
		            cat("Iteration:", beta.iter, "loglik = ", loglik[beta.iter-1] , "\n")
		        }
		    }

	        if(beta.iter > 3){
	            if(abs(loglik[beta.iter-3] - loglik[beta.iter-1])/orig < epsilon.beta)
	                iterating <- FALSE
	        }
	    }
	    
	    return(list("beta" = betak, "objfunc" = loglik[1:(beta.iter - 2)]))
	}




	# - Cross-validation for tau and lambda grid (initialize at lasso solution)
	MCMVR.cv <- function(Z, Y, tau.vec, nlambda, weight = NULL, standardize = FALSE, nfolds = NULL, delta = .01, tol = 1e-8, quiet=TRUE, inner.quiet=TRUE){
	  
	  	# ---- preliminaries
		n <- dim(Z)[1]
		p <- dim(Z)[2]
		q <- dim(Y)[2]
		if(is.null(weight)){
	    	weight <- matrix(1, nrow=p, ncol=q)
	    }

		# ----- standardize if necessary 
		if(!standardize){
			z <- Z - rep(1, n)%*%t(apply(Z, 2, mean))
			ztz <- crossprod(z)
			y <- Y - rep(1, n)%*%t(apply(Y, 2, mean))
			zty <- crossprod(z, y)
		} else {
			z <- (Z - rep(1, n)%*%t(apply(Z, 2, mean)))/(rep(1, n)%*%t(apply(Z, 2, sd)))
			ztz <- crossprod(z)
			y <- (Y - rep(1, n)%*%t(apply(Y, 2, mean)))/(rep(1, n)%*%t(apply(Y, 2, sd)))
			zty <- crossprod(z, y)
		}

		# ---- select tuning parameters for lambda 
		lambda.vec <- rep(0, length=nlambda)
		lambda.max <-1.01*dim(z)[1]^(-1)*max(abs(crossprod(z,y))) + 1e-6
		lambda.min <- delta*lambda.max
		for(kk in 1:nlambda){
			lambda.vec[kk] <- lambda.max^((nlambda-kk)/(nlambda-1))*lambda.min^((kk-1)/(nlambda-1))
		}

		# --- for each lambda, initialize largest tau value at the lasso solution 
		beta.full <- Matrix(0, nrow = p*q*length(tau.vec), ncol = length(lambda.vec), sparse=TRUE)
		sparsity.mat <- matrix(0, length(lambda.vec), length(tau.vec))

		for(jj in 1:length(tau.vec)){
			beta.old <- matrix(0, nrow=p, ncol=q)

			for(kk in 1:length(lambda.vec)){
		   		temp <- Acc.Prox.L1(z, y, ztz, zty, beta.old, tau.vec[jj], lambda.vec[kk]/tau.vec[jj], 
		   			L = NULL, weight = weight, inner.quiet = inner.quiet, epsilon.beta = tol, max.iter.inner = 1e4)
		   		beta.old <- temp$beta
				beta.full[(p*q*(jj-1)+1):(p*q*jj),kk] <- c(beta.old)
		     	sparsity.mat[kk,jj] <- sum(beta.old!=0)
		     	if(!quiet){
		     		cat(jj, kk, ": non-zero = ", sum(beta.old!=0), "\n")
		     	}
		  	}
		}
		
		if(!is.null(nfolds)){

			# --- perform cross validation
			fold <- sample(rep(1:nfolds, length=n))
			cv.index <- split(1:n, fold)
			errs_wpred <- array(0, dim=c(length(lambda.vec),length(tau.vec), nfolds))
			errs_pred <- array(0, dim=c(length(lambda.vec),length(tau.vec), nfolds))

			for(k in 1:nfolds){

				if(!standardize){
					
					# --- center Z and Y
					ntrain <- dim(Z[-cv.index[[k]],])[1]
					z.inner <- Z[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(Z[-cv.index[[k]],], 2, mean))
					ztz.inner <- crossprod(z.inner)
					y.inner <- Y[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(Y[-cv.index[[k]],], 2, mean))
					var.y <- apply(y.inner, 2, var)
					zty.inner <- crossprod(z.inner, y.inner)

					n.test <- length(cv.index[[k]])
					z.test <- Z[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(Z[-cv.index[[k]],], 2, mean))
					y.test <- Y[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(Y[-cv.index[[k]],], 2, mean))
				
				} else {	

					# --- center Z and Y --------------------
					ntrain <- dim(Z[-cv.index[[k]],])[1]
					z.inner <- (Z[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(Z[-cv.index[[k]],], 2, mean)))/(rep(1, ntrain)%*%t(apply(Z[-cv.index[[k]], ], 2, sd)))
					ztz.inner <- crossprod(z.inner)
					y.inner <- (Y[-cv.index[[k]], ] - rep(1, ntrain)%*%t(apply(Y[-cv.index[[k]],], 2, mean)))/(rep(1, ntrain)%*%t(apply(Y[-cv.index[[k]], ], 2, sd)))
					var.y <- apply(y.inner, 2, var)
					zty.inner <- crossprod(z.inner, y.inner)

					# ------------------------------------------------
					n.test <- length(cv.index[[k]])
					z.test <- (Z[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(Z[-cv.index[[k]],], 2, mean)))/(rep(1, n.test)%*%t(apply(Z[-cv.index[[k]], ], 2, sd)))
					y.test <- (Y[cv.index[[k]], ] - rep(1, n.test)%*%t(apply(Y[-cv.index[[k]],], 2, mean)))/(rep(1, n.test)%*%t(apply(Y[-cv.index[[k]], ], 2, sd)))
				
				}

				# ------------------------------------------------
				for(jj in 1:length(tau.vec)){
					beta.old <- matrix(0, nrow=p, ncol=q)
					for(kk in 1:length(lambda.vec)){
			   			temp <- Acc.Prox.L1(z.inner, y.inner, ztz.inner, zty.inner, beta.old, tau.vec[jj], lambda.vec[kk]/tau.vec[jj], L = NULL, weight = weight, inner.quiet = inner.quiet, epsilon.beta = tol, max.iter.inner = 1e4)
			   			errs_pred[kk, jj, k] <- sum((y.test - z.test%*%temp$beta)^2)/(n.test*q)
			   			beta.old <- temp$beta
			  		}
			 	}
			 	if(!quiet){
			 		cat("Through CV fold", k, "\n")
			 	}
			}

		} else {
			errs_pred <- NULL
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
		     		"Y.offset" = apply(Y, 2, mean), 
		     		"Z.offset" = apply(Z, 2, mean),
		     		"Y.sd" = apply(Y, 2, sd), 
		     		"Z.sd" = apply(Z, 2, sd),
					"lambda.vec" = lambda.vec, 
					"tau.vec" = tau.vec,
					"tau.min" = tau.min,
					"lam.min" = lam.min,
					"standardize" = standardize)

		class(fit) <- "MC_MVR"
		return(fit)
	}
	


	MCMVR.predict <- function(Znew, fit, lambda = NULL, tau = NULL){

		if(is.null(lambda) && is.null(tau)){
			tau <- fit$tau.min
			lambda <- fit$lam.min
			if(is.null(fit$tau.min)){
				stop('No tuning parameters selected by CV')
			}

		} 
		

		tau.ind <- which(fit$tau.vec == tau)
		lam.ind <- which(fit$lambda.vec == lambda)
		p <- length(fit$Z.offset)
		q <- length(fit$Y.offset)
		if(!fit$standardize){
			beta.vec <- fit$beta[(p*q*(tau.ind-1)+1):(p*q*tau.ind),lam.ind]
			beta.mat <- matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)
		} else {
			beta.vec <- fit$beta[(p*q*(tau.ind-1)+1):(p*q*tau.ind),lam.ind]
			beta.mat <- tcrossprod(fit$Z.sd, rep(1, q))*matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)*tcrossprod(rep(1, p), fit$Y.sd)
		}
		# --- get intercept 	
		B0 <- fit$Y.offset - crossprod(beta.mat, fit$Z.offset)
		if(dim(Znew)[1] > 1){
			preds <- tcrossprod(rep(1, dim(Znew)[1]), B0) + tcrossprod(Znew, t(beta.mat))
		} else {
			preds <- B0 + crossprod(beta.mat, Znew)
		}

		return(list("pred" = preds, "beta0" = B0, "beta" = beta.mat))

	}

	MCMVR.coef <- function(fit, lambda = NULL, tau = NULL){

		if(is.null(lambda) && is.null(tau)){
			tau <- fit$tau.min
			lambda <- fit$lam.min
			if(is.null(fit$tau.min)){
				stop('No tuning parameters selected by CV')
			}
		} 

		tau.ind <- which(fit$tau.vec == tau)
		lam.ind <- which(fit$lambda.vec == lambda)
		p <- length(fit$Z.offset)
		q <- length(fit$Y.offset)
		if(!fit$standardize){
			beta.vec <- fit$beta[(p*q*(tau.ind-1)+1):(p*q*tau.ind),lam.ind]
			beta.mat <- matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)
		} else {
			beta.vec <- fit$beta[(p*q*(tau.ind-1)+1):(p*q*tau.ind),lam.ind]
			beta.mat <- tcrossprod(fit$Z.sd, rep(1, q))*matrix(beta.vec, byrow=FALSE, nrow=p, ncol=q)*tcrossprod(rep(1, p), fit$Y.sd)
		}
		# --- get intercept 	
		B0 <- fit$Y.offset - crossprod(beta.mat, fit$Z.offset)
		return(list("beta0" = B0, "beta" = beta.mat))
	}
