\name{methylation}
\alias{methylation}
\docType{data}

\title{CpG island DNA methylation array data}

\description{CpG island DNA methylation array data of a neuro-ectodermal cell line that was treated with a demethylating agent}

\usage{
	data(methylation)
}

\format{
  	"RGList" as defined in the package limma containing data from CpG island DNA methylation array data of a neuro-ectodermal cell line that was 	treated with a demethylating agent. The element "weights" of the "RGList" contains the subset of invariant fragments, those without methylation-sensitive restriction sites, as a logical vector. 
}

\details{
   	The data is extracted from a larger experiment described in van Iterson et al. Because the data is from a high-dense tiling array a random subset of the data was chosen for convenience in making the vignette.
}

\references{
van Iterson et al. A novel and fast normalization method for high-density arrays. SAGMB 11(4) 2012.
\emph{in preparation}. 
}

\examples{
	data(methylation)
}

\keyword{datasets}
