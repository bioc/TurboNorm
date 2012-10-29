################################################################################
##Fast P-spline smoothing ... 
##
##
##
################################################################################
#interface to Paul's C implementation for fast summation
binsum <- function(x, y = 0 * x + 1, n = max(as.integer(x)))
  {
    ## Sum values of y in bins given by x
    
    ## Create working arrays
    a <- as.integer(x)
    b <- rep(0.0, n + 1)
    m <- as.integer(length(x))
    ## Call C function
    som <- .C("binsum", a, y, b, m, CLASSES = c("integer", "numeric", "numeric", "integer"))[[3]]
    return(som[-1])
  }

#the work horse based on P.H.C. Eilers' implementation
turbotrend <- function(x, y, w = rep(1, length(y)), n = 100, lambda=10^seq(-10, 10, length=1000), iter = 0, method=c("original", "demmler"))
{
  my.call <- match.call()	
  method <- match.arg(method)	

  maxLambda <- 10^12
	
  if(length(x) != length(y))
    stop("x and y are of different length!")
	
  if(length(x) != length(w))
    stop("weights must be of same length as data!")	

  ##order x and y by x
  ##x <- x[order(x)]
  ##y <- y[order(x)]
  ##xmin <- x[1]
  ##xmax <- x[length(x)]
			
  xmin <- min(x)
  xmax <- max(x)
  dx <- 1.0001 * (xmax - xmin) / n
  xt <- xmin + ((1:n) - 0.5) * dx
  xi <- floor((x - xmin) / dx + 1)

  ##Collect counts and sums (handle empty bins correctly)
  ##tel <- tapply(0 * xi + 1, xi, sum) #BTB
  ##tel <- tapply(w, xi, sum)  #BTWB
  ##p <- as.integer(names(tel))
  ##u <- rep(0, n)
  ##u[p] <- tel
  ##using Paul's C implementation
  u <- binsum(xi, as.numeric(w))       #as binners.c expects "numeric"
	
  ##som <- tapply(y, xi, sum) #BTy
  ##som <- tapply(w*y, xi, sum) #BTWy
  ##v <- rep(0, n)
  ##
  ##v[p] <- som
  ##using Paul's C implementation
  v <- binsum(xi, y*w)

  ## Build penalty and solve system
  D <- diff(diag(n), diff = 2)
	
  ## Perform generalized cross-validation in order to find the optimal lambda if not given
  if(length(lambda) > 1){				

    maxLambda <- lambda[length(lambda)]
		
    ##old slower search for optimal lambda
    ##gcv <- switch(method, 
    ##                 original = sapply(lambda, function(x) originalBasis(x, xi, y, u, v, D)$gcv),    
    ##                 demmler = sapply(lambda, function(x) DemmlerReinschBasis(x, xi, y, u, v, D)$gcv))			

    ##lambda <- lambda[which.min(gcv)]

    ##is this faster?
    GCV.org <- function(x, xi, y, u, v, D) originalBasis(x, xi, y, u, v, D)$gcv
    GCV.DR <- function(x, xi, y, u, v, D) DemmlerReinschBasis(x, xi, y, u, v, D)$gcv
    lambda <- switch(method, 
                     original = optimize(GCV.org, interval=range(lambda), xi, y, u, v, D)$minimum,
                     demmler = optimize(GCV.DR, interval=range(lambda), xi, y, u, v, D)$minimum)
  }
	
  if(lambda <= 0 & lambda > maxLambda)
    warning("lambda <= 0 or > maxLambda!")

  solved <- switch(method, 
                   original = originalBasis(lambda, xi, y, u, v, D),
                   demmler = DemmlerReinschBasis(lambda, xi, y, u, v, D))              

  ##catch error
  if(inherits(solved, "try-error"))
    return(solved)
  
  ## Robustifying iterations 
  if(iter > 0) {
    P <- lambda * t(D) %*% D
    yhat <- solved$coeff[xi]
    for (it in 1:iter) {
      r <- y - yhat
      s <- median(abs(r))
      t <- r / (5 * s + 1e-4)
      wr <- (1 - t ^ 2) ^ 2
      wr[abs(t) > 1] <- 0
      u <- binsum(xi, wr)
      v <- binsum(xi, wr * y)
      W <- diag(u)
      z <- solve(W + P, v)
      yhat <- z[xi]
    }
  } 
  else {
    yhat <- solved$coeff[xi]
  }
  
  ## Return list with results
  object <- list(x = x, y = yhat, w=w, n=n, xtrend = xt, ytrend = solved$coeff, lambda=solved$lambda, iter=iter, gcv=solved$gcv, 
                 edf=solved$edf, method=method, call=my.call) 
  class(object) <- "turbotrend"
  object
}

################################################################################
##print method in a similar way as print.smoothPspline from pspline-package
##
##
##
################################################################################

print.turbotrend <- function(x, ...) {
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
      dput(cl)
   }
   cat("\nEffective degrees of freedom:", format(x$edf), "\n")
   cat("Number of bins:", format(x$n), "\n")
   cat("Penalty value:", format(x$lambda),  "\n") 
   cat("Number of robustifying iterations:", format(x$iter), "\n")              
   cat("GCV :", format(x$gcv), "\n")
   invisible(x)
}


################################################################################
##solve 
##
##
################################################################################
originalBasis <- function(lambda, xi, y, u, v, D)
{
  ##implementation without calculating the 'smoother' matrix explicitly because that involves nxn matrices i.s.o mxm
  tmp <- diag(u) + lambda * t(D) %*% D
  z <- try(solve(tmp, v))
        
  ##catch error
  if(inherits(z, "try-error"))
    return(z)

  ##cyclic permutation of matrices as it doesn't change the trace but makes the computation more efficient
  edf <- sum(diag(solve(tmp) %*% diag(u))) 
  ##gcv <- crossprod(y - z[xi])/(length(u) - edf)^2
	
  gcv <- sum((y - z[xi])^2)/(length(u) - edf)^2
 
  ##z0 <- solve(diag(u), v) #penalty = 0
  ##var0 <- sum((y - z0[xi])^2)/(length(u))^2 #residuals variances	
  ##aic <- sum((y - z[xi])^2)/var0 + 2*length(u)*log(sqrt(var0)) - 2*length(y)*log(2*pi)

  list(lambda=lambda, edf=edf, gcv=gcv, coeff=z)
}

################################################################################
##solve 
##
##
################################################################################
DemmlerReinschBasis <- function(lambda, xi, y, u, v, D)
{
	X <- chol(diag(u)) ## X^(-t)X^(-1) == diag(u)
	X <- diag(1/diag(X)) #inverse

	E <- eigen(X%*% t(D)%*%D%*%t(X), symmetric=TRUE)	

	T <- E$vectors ##T should be orthogonal and all.equal(T%*%diag(c)%*%t(T), X%*% t(D)%*%D%*%t(X))
	c <- E$values	
	
	z <- t(X)%*%T%*% solve(diag(1 + lambda*c), t(T)%*%X%*%v)
	edf <- sum(1/(1 + lambda*c))
	#gcv <- crossprod(y - z[xi])/(length(u) - edf)^2
	gcv <- sum((y - z[xi])^2)/(length(u) - edf)^2

	list(lambda=lambda, edf=edf, gcv=gcv, coeff=z)
}


  
     
