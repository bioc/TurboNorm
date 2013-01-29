\name{panel.pspline}
\alias{panel.pspline}
\title{Panel function for adding a P-spline smoothed curves to a lattice graphics panel}
\description{The function panel.pspline is similar to panel.loess but show the P-spline smoothed curve.}

\usage{panel.pspline(x, y, weights = rep(1, length(y)), nintervals = 100, type, horizontal = FALSE, col.line=1, lty=1, lwd=1, ...)}

\arguments{
  \item{x, y}{vectors giving the coordinates of the points in the scatter plot}
  \item{weights}{vector of weights of with same length as the data for a weighted smoothing. Default all weights are 1.}
  \item{nintervals}{an integer indicating the number of intervals equal to 1 + number of knots. Currently the intervals must be langer than 10.}
  \item{type}{see \code{\link[lattice]{panel.loess}}}
  \item{horizontal}{see \code{\link[lattice]{panel.loess}}}
  \item{col.line, lty, lwd}{line colour, type and width that will be used in the plots, defaults are col=1, lty=1 and lwd=1.}
  \item{\dots}{see \code{\link[lattice]{panel.loess}}}
}

\details{?panel.loess}

\references{
  Deepayan Sarkar (2009). lattice: Lattice Graphics. R package version
  0.17-26. http://CRAN.R-project.org/package=lattice
  
  van Iterson M, Duijkers FA, Meijerink JP, Admiraal P, van Ommen GJ, Boer JM, van Noesel MM, Menezes RX (2012). A novel and fast normalization method for
  high-density arrays. SAGMB, 11(4).
  
  Paul .H.C. Eilers and Brain D. Marx (1996). Flexible smoothing with B-splines and Penalties. Statistical Science, Vol 11, No. 2, 89-121.
}

\author{Maarten van Iterson and Chantal van Leeuwen}

\seealso{\code{\link[lattice]{panel.loess}}}

\examples{

library(marray)
library(lattice)
data(swirl)
data <- data.frame(M=as.vector(maM(swirl)), A=as.vector(maA(swirl)), Sample=rep(paste("Array", 1:4), each=nrow(swirl)))

xyplot(M~A|Sample, data=data,
                        panel = function(x, y) {
            panel.grid(h=-1, v= 2)
            panel.xyplot(x, y)
            panel.loess(x, y, span=0.25, col="black")
                                                panel.pspline(x, y, col="red", lwd=2)})

}

\keyword{hplot}
