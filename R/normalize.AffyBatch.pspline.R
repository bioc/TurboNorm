################################################################################
##taken from the affy package: batch loess normalization function
##replaced the loess with the P-spline smoother
##
##
################################################################################

normalize.AffyBatch.pspline <- function (abatch, type = c("together", "pmonly", "mmonly", "separate"), ...) 
{
    type <- match.arg(type)
    if (type == "separate") {
        Index <- unlist(indexProbes(abatch, "pm"))
        intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], ...)
        Index <- unlist(indexProbes(abatch, "mm"))
        intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], ...)
    }
    else if (type == "together") {
        Index <- unlist(indexProbes(abatch, "both"))
        intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], ...)
    }
    else if (type == "pmonly") {
        Index <- unlist(indexProbes(abatch, "pm"))
        intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], ...)
    }
    else if (type == "mmonly") {
        Index <- unlist(indexProbes(abatch, "mm"))
        intensity(abatch)[Index, ] <- normalize.pspline(intensity(abatch)[Index,], ...)
    }
    return(abatch)
}
