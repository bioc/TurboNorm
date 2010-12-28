\name{TurboNorm-package}
\alias{TurboNorm-package}
\docType{package}

\title{A fast scatterplot smoother with applications for microarray normalization}

\description{A fast scatterplot smoother based on B-splines with second order difference penalty. Functions for microarray normalization of single-colour data i.e. Affymetrix/Illumina and two-colour data supplied as marray MarrayRaw-objects or limma RGList-objects are available.}

\details{
\tabular{ll}{
Package: \tab TurboNorm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-09-15\cr
License: \tab LGPL\cr
LazyLoad: \tab yes\cr
}

This package contains a implementation of piecewise constant P-splines of Eilers and Marx (1996) that can be used for normalization of either single- or two-colour data. For two-colour data objects of type \code{RGList} from the \code{limma} package and \code{MarrayRaw} from the package \code{marray} can be normalized using the function \code{pspline()}. For single colour microarray data wrapper functions are written based on the \code{affy} package functions \code{normalize.loess()} and \code{normalize.AffyBatch.loess()} namely \code{normalize.pspline()} and \code{normalize.AffyBatch.pspline()}. Also a \code{panel.pspline()} is available for adding the smoothed curve to \code{lattice} graphics panels.

}
\author{ 
Chantal van Leeuwen and Maarten van Iterson 
Maintainer: Maarten van Iterson<M.van_Iterson.HG@lumc.nl>}

\references{
Paul .H.C. Eilers and Brain D. Marx (1996). Flexible smoothing with B-splines and Penalties. Statistical Science, Vol 11, No. 2, 89-121.
}

\note{The package pspline (S original by Jim Ramsey, R port by Brian Ripley) implements the B-spline/Natural Cubic Spline smoother}

\seealso{\code{\link{turbotrend}, \link{pspline}, \link{normalize.pspline}, \link{normalize.AffyBatch.pspline}, \link{panel.pspline}}}

\keyword{package}
