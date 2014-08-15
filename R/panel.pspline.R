################################################################################
##Mimic affy's pairwise loess normalization function
##
##
##
################################################################################

panel.pspline <- function (x, y, weights = rep(1, length(y)), nintervals = 100, type, horizontal = FALSE, col.line=1, lty=1, lwd=1, ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 1)
    return()
  plot.line <- trellis.par.get("plot.line")

  if (horizontal) {
    smooth <- turbotrend(y, x, w = weights, n = nintervals)
    panel.lines(x = smooth$x[order(smooth$x)], y = smooth$y[order(smooth$x)], col.line=col.line, lty=lty, lwd=lwd, ...)
  }
  else {
    smooth <- turbotrend(x, y, w = weights, n = nintervals)
    panel.lines(x = smooth$x[order(smooth$x)], y = smooth$y[order(smooth$x)], col.line=col.line, lty=lty, lwd=lwd, ...)
  }
}
