\name{pspline}
\alias{pspline}
\title{Function for two-colour microarray data normalization using the P-splines}
\description{Wrapper function for two colour microarray data normalization using the P-spline smoother suitable for a RGList- or MarrayRaw-objects.}

\usage{pspline(object, background = c("none", "substract"), weights = NULL, nintervals = 100, subset=NULL, showArrays = 0, verbose=FALSE, line.col=2,  line.lty=1, line.lwd=2, ...)}

\arguments{
  \item{object}{either a RGList or an MarrayRaw-object.}
  \item{background}{for background substraction use 'substract'. Default is no backgroud substraction.}
  \item{weights}{vector of weights that will be used a for a weighted normalization. The default \code{NULL} assumse equal weight \code{1} for all data points.}
  \item{nintervals}{number of bins in which the data will be divided. The default is 100 bins.}
  \item{showArrays}{either a integer( > 0)  or a vector of integers indicating the arrays for which a MA-plot will be produced.}
  \item{subset}{subset of the data on which the normalization will be based. A special case of weighted normalization.}
	\item{verbose}{if \code{TRUE} gives additional information on the fit}.
	\item{line.col, line.lty, line.lwd}{line colour, type and width that will be used in the plots, defaults are col=2, lty=1 and lwd=2.}
	\item{\dots}{additional graphical arguments for plotting.}
}

\details{if necessary?}

\value{  
  The value that will be returned is either a MAList or MarrayNorm-object dependening on the input type.
}

\references{ 
Paul .H.C. Eilers and Brain D. Marx (1996). Flexible smoothing with B-splines and Penalties. Statistical Science, Vol 11, No. 2, 89-121.
}

\author{ Chantal van Leeuwen and Maarten van Iterson}

\seealso{\code{\link[limma]{normalizeWithinArrays}, \link[marray]{maNormMain}}}

\examples{

library(marray)
data(swirl)

x <- pspline(swirl, showArrays=2, pch=20, col="grey")
x <- pspline(swirl, showArrays=2:4, line.col="green")

}

\keyword{smooth}

