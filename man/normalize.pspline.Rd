\name{normalize.pspline}
\alias{normalize.pspline}
\alias{normalize.AffyBatch.pspline}
\title{Functions for single-colour microarray data normalization using the P-splines}
\description{Modified version of normalize.loess and normalize.AffyBatch.pspline from the affy package uses the P-spline smoother in stead of the loess algorithm}

\usage{
normalize.pspline(mat, epsilon = 10^-2, maxit = 1, log.it = TRUE,
    verbose = TRUE, weights = rep(1, nrow(mat)), ...)
normalize.AffyBatch.pspline(abatch,
    type=c("together","pmonly","mmonly","separate"), ...)
}

\arguments{
  \item{mat}{a matrix with columns containing the values of the chips to normalize.}
  \item{abatch}{an \code{\link[affy:AffyBatch-class]{AffyBatch}} object.}
  \item{epsilon}{a tolerance value (supposed to be a small value - used as a stopping criterion).}
  \item{maxit}{maximum number of iterations.}
  \item{log.it}{logical. If \code{TRUE} it takes the log2 of \code{mat}}
  \item{verbose}{logical. If \code{TRUE} displays current pair of chip being worked on.}
  \item{weights}{For weighted normalization. The default is NULL, so there are no weights used.}
	\item{type}{A string specifying how the normalization should be applied. See details for more.}
  \item{\dots}{Graphical parameters can be supplied.}
}

\details{
This function is a modified version of the function \code{\link[affy]{normalize.loess}} from the affy package. In stead of the loess algorithm the function uses the P-spline algorithm. 
The type argument should be one of \code{"separate","pmonly","mmonly","together"} which indicates whether to normalize only one probe type(PM,MM) or both together or separately.
}

\value{Normalized \code{\link[affy:AffyBatch-class]{AffyBatch}} }

\references{ 
Laurent Gautier, Leslie Cope, Benjamin M. Bolstad and Rafael A. Irizarry (2004). affy -analysis of Affymetrix GeneChip data at the probe level. Bioinformatics, Vol. 20, no. 3, 307-315.

Paul .H.C. Eilers and Brain D. Marx (1996). Flexible smoothing with B-splines and Penalties. Statistical Science, Vol 11, No. 2, 89-121.
}

\author{Maarten van Iterson and Chantal van Leeuwen}

\seealso{\code{\link[affy]{normalize.loess}}}

\examples{

library(affydata)
data(Dilution)
PM <- log2(pm(Dilution[,c(1,3)]))
M <- PM[,1]-PM[,2]
A <- 0.5*(PM[,1]+PM[,2]) 

nPM <- log2(normalize.pspline(pm(Dilution[,c(1,3)])))
nM <- nPM[,1]-nPM[,2]
nA <- 0.5*(nPM[,1]+nPM[,2]) 

par(mfcol=c(2,1))
plot(M~A)
plot(nM~nA)

norm <- normalize.AffyBatch.pspline(Dilution, type="pmonly")

}

\keyword{smooth}

