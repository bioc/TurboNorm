################################################################################
##taken from the affy package: batch loess normalization function
##replaced the loess with the P-spline smoother
##
##
################################################################################

normalize.AffyBatch.pspline <- function (abatch, type = c("together", "pmonly", "mmonly", "separate"), ...)
{
  type <- match.arg(type)

  apar <- list(...)
  if(length(apar) != 0)
    {
      if(any(weights %in% names(apar)))
        weights <- apar$weights
      else
        weights <- rep(1, nrow(exprs(abatch)))      
    }
   else
        weights <- rep(1, nrow(exprs(abatch)))

  if (type == "separate") {
    Index <- unlist(indexProbes(abatch, "pm"))
    intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], weights = weights[Index])
    Index <- unlist(indexProbes(abatch, "mm"))
    intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], weights = weights[Index])
  }
  else if (type == "together") {
    Index <- unlist(indexProbes(abatch, "both"))
    intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], weights = weights[Index])
  }
  else if (type == "pmonly") {
    Index <- unlist(indexProbes(abatch, "pm"))
    intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], weights = weights[Index])
  }
  else if (type == "mmonly") {
    Index <- unlist(indexProbes(abatch, "mm"))
    intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], weights = weights[Index])
  }
  return(abatch)
}
