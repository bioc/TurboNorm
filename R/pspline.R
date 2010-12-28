################################################################################
##wrapper function for two colour microarray normalization
##input either MarrayRaw object or limma RGList
##
##
################################################################################

pspline <- function(object, background = c("none", "substract"), weights = NULL, nintervals = 100, subset = NULL, showArrays=0, verbose=FALSE, line.col=2, line.lty=1, line.lwd=2, ...)
{
	#parameter checking
	background <- match.arg(background)
	
	#determine class of object
	RawData <- switch(class(object),         
	  				  		marrayRaw = object,
	   					RGList = as(object, "marrayRaw"),
	   					stop("Unknown type: ", class(object)))

	#set background too zero
	if(background == "none")
		maRb(RawData) <- maGb(RawData) <- matrix(0, nrow=nrow(RawData), ncol=ncol(RawData))

	#extract M- and A-values and thus optionally background corrected
	A <- maA(RawData)	
	M <- maM(RawData)
	
	#determine weight
	if(is.null(weights))
		weights <- rep(1, nrow(A))

	#test nintervals
	if(nintervals < 10 | nintervals > nrow(RawData))
		stop("The nintervals need to be between 10 and the number of datapoints!")

	#determine subset
	if(!is.null(subset))
		weights[!subset] <- 0

	#test weights
	if(length(weights) != nrow(A))
		stop("Weights are of wrong length!")

	#construct marrayNorm-object for holding the results
	NormData <- as(RawData, "marrayNorm")

	nArrays <- length(showArrays)
	if(min(showArrays) < 0 | max(showArrays) > ncol(A))
		stop("Wrong number of arrays selected for plotting!")

	if(nArrays > 0 & showArrays[1] != 0){
		layout(matrix(1:nArrays, nrow=n2mfrow(nArrays)[1], ncol=n2mfrow(nArrays)[2]))
	}

	#apply P-spline normalization
	for(j in 1:ncol(A))
	{
		tt <- turbotrend(A[,j], M[,j], w = weights, n = nintervals) 

		if(verbose)
			print(tt)	

		maM(NormData)[,j] <- M[,j] - tt$y

		if(is.element(j, showArrays))
		{
			plot(A[,j], M[,j], xlab="A", ylab="M", main=paste("Raw data array:", j), ...)			
			lines(tt$x[order(tt$x)], tt$y[order(tt$x)], col=line.col, lty=line.lty, lwd=line.lwd)
		}
	}

	#conversion too normalized object 
	#maybe set weights to zero in MAList otherwise weights will be used in lmFit (see limma userguide page 29)
	conversion <- list(RGList="MAList", marrayRaw="marrayNorm")
	as(NormData, conversion[[class(object)]])
}

