################################################################################
##the work horse based on P.H.C. Eilers' implementation
##
##
##
################################################################################

turbotrend <- function(x, y, w = rep(1, length(y)), n = 100, lambda=10^seq(-10, 10, length=1000), method=c("original", "demmler"))
{
	my.call <- match.call()	
	method <- match.arg(method)	

	maxLambda <- 10^12
	
	if(length(x) != length(y))
		stop("x and y are of different length!")
	
	if(length(x) != length(w))
		stop("weights must be of same length as data!")	

	xmin <- min(x)
  xmax <- max(x)
  dx <- 1.0001 * (xmax - xmin) / n
  xt <- xmin + ((1:n) - 0.5) * dx
  xi <- floor((x - xmin) / dx + 1)

	#this should optimize the tapply
	xi <- as.integer(xi)

  # Collect counts and sums (handle empty bins correctly)
  #tel <- tapply(0 * xi + 1, xi, sum) #BTB
	tel <- tapply(w, xi, sum)  #BTWB
  p <- as.integer(names(tel))
  u <- rep(0, n)
  u[p] <- tel

  #som <- tapply(y, xi, sum) #BTB
	som <- tapply(w*y, xi, sum) #BTWy
  v <- rep(0, n)
  v[p] <- som

  # Build penalty and solve system
  D <- diff(diag(n), diff = 2)
	
	# Perform generalized cross-validation in order to find the optimal lambda if not given
  if(length(lambda) > 1){				

		maxLambda <- lambda[length(lambda)]
		
		#old slower search for optimale lambda
		#gcv <- switch(method, 
    #                 original = sapply(lambda, function(x) originalBasis(x, xi, y, u, v, D)$gcv),    
		#                 demmler = sapply(lambda, function(x) DemmlerReinschBasis(x, xi, y, u, v, D)$gcv))			

	  #lambda <- lambda[which.min(gcv)]

		#is this faster?
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

  # Return list with results
  object <- list(x = x, y = solved$coeff[xi], w=w, n=n, xtrend = xt, ytrend = solved$coeff, lambda=solved$lambda, gcv=solved$gcv, 
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
	#implementation without calculating the 'smoother' matrix explicitly because that involves nxn matrices i.s.o mxm
	tmp <- diag(u) + lambda * t(D) %*% D
	z <- solve(tmp, v)	 
	#cyclic permutation of matrices as it doesn't change the trace but makes the computation more efficient
	edf <- sum(diag(solve(tmp) %*% diag(u))) 
	#gcv <- crossprod(y - z[xi])/(length(u) - edf)^2
	gcv <- sum((y - z[xi])^2)/(length(u) - edf)^2

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

