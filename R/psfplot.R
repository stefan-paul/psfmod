##' Plotting function
##'
##' @param out output of \code{\link{psfmod}}
##' @param legend logigal, determining whether a legend should be drawn
##' @param ... futher arguments passed to matplot
##' @return plot
##' @export

psfplot <- function(out, legend  = FALSE, ...){
  mycol <- rainbow(6)
  matplot(out[, 1], out[,c(-1,-8,-9,-10)], type = "l", col = mycol, lty = 1,
        xlab = "Time (years)", ylab = "Values", ...)
  if(legend){
  legend("topright", lty = 1, col = mycol, legend = c("Biomass A", "Biomass B", "N", "Litter A", "Litter B", "Light"))
}
    }


#' @export
plot.psfmod <- function(out, ...){
  psfplot(out, ...)
}
