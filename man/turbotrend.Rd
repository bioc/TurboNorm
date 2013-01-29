\name{turbotrend}
\alias{turbotrend}
\title{turbotrend: a fast scatterplot smoother}

\description{A fast scatterplot smoother based on B-splines with second order difference penalty}

\usage{turbotrend(x, y, w = rep(1, length(y)), n = 100, lambda=10^seq(-10, 10, length=1000), iter=0, method=c("original", "demmler"))}

\arguments{
  \item{x,y}{vectors giving the coordinates of the points in the scatter plot.}
  \item{w}{vector of weights of with same length as the data for a weighted smoothing. Default all weights are 1.}
  \item{n}{an integer indicating the number of intervals equal to 1 + number of knots. Currently the intervals must be langer than 10.}
  \item{lambda}{Optionally a user-defined penalty parameter can be provided, if not generalized cross-validation is used to find the optimal penalty parameter.}
  \item{iter}{Number of robustifying iterations similar as lowess.}
  \item{method}{method for solving the system of linear equations either using the data in the original space or transformed to the Demmler-Reinsch basis.}
}

\details{some details about implementation}

\value{
  An object of type \code{pspline} is returned as a list with the following items:
  \item{x}{original data vector x}
  \item{y}{fitted y-values with same length as vector x}
  \item{w}{vector of weights}
  \item{n}{number of bins}
  \item{ytrend}{binnend fitted y-values}
  \item{xtrend}{binned x-values}
  \item{lambda}{if scalar penalty parameter used else if vector of two lower and upper bound of the grid}
  \item{iter}{number of robustifying iterations}
  \item{gcv}{generalized cross-validation}
  \item{edf}{effective degrees of freedom (trace of the smoother matrix)}
  \item{call}{function call which produced this output}
}

\references{
  van Iterson M, Duijkers FA, Meijerink JP, Admiraal P, van Ommen GJ, Boer
  JM, van Noesel MM, Menezes RX (2012). A novel and fast normalization method for
  high-density arrays. SAGMB, 11(4).

  Paul .H.C. Eilers and Brain D. Marx (1996). Flexible smoothing with B-splines and Penalties. Statistical Science, Vol 11, No. 2, 89-121.
}

\author{Maarten van Iterson, Chantal van Leeuwen}


\seealso{\code{\link{loess}},\code{\link{lowess}}, \code{\link{smooth}}, \code{\link{smooth.spline}} and  \code{\link[pspline]{smooth.Pspline}}}

\examples{

library(marray)
data(swirl)

x <- maA(swirl)[,1]
y <- maM(swirl)[,1]
xord <- x[order(x)]
yord <- y[order(x)]

plot(xord, yord, main = "data(swirl)  &  smoothing splines + lowess")
lines(turbotrend(xord, yord), col = "red", lwd=2)
lines(smooth.spline(xord, yord), col = "green", lwd=2)
lines(lowess(xord, yord), col = "purple", lwd=2)
legend("topleft", c("piecewise constant P-splines", "Cubic B-splines", "lowess"), text.col=c("red","green","purple"), bty="n")

}

\keyword{smooth}
\keyword{regression}
