#' Fits a mixture model on an array of images
#'
#' This function computes a k component gamma mixture model marginally for each channel and returns a mixfit object to perform image normalization using the XX function.
#' @param image Required character vector of image locations or raster object.
#' @param mask Required character vector of image location or raster object of the tissue mask.
#' @param k number of mixture components. Can be thought of as the number of tissue classes distinguishable by this channel. Probably better to have too many than too few. With too too many convergence issues may occur.
#' @param lambda k length vector of mixing proportions.
#' @param alpha k length vector of shape parameters
#' @param beta k length vector of beta parameters
#' @param mc.cores Number of cores for parallel things. With large images and a large number of cores this may create memory issues. With many channels, it may be better to distributed computing over a cluster.
#' @param subsamp Proportion of pixels within mask to use to estimate mixture model. gammamixEM can fail with a large number of pixels.
#' @param maxpix Maximum number of pixels to use to estimate mixture model.
#' @param ... Other arguments passed to gammamixEM
#' @keywords gamma mixture model
#' @return Returns a mixfit object with the following values:
#' @importFrom raster raster stack
#' @importFrom mixtools gammamixEM
#' @importFrom parallel mclapply
#' @export
fitmix = function(image, mask, lambda=NULL, k=3, alpha=c(.333, 3.33, 33.3), beta=c(c(.222, 2.22, 22.2)), mc.cores=getOption("mc.cores", 2L), subsamp=1, maxpix=NULL, ...){
  if(all(is.character(image))) image = raster::stack(image)
  # will error if multiple masks are passed
  if(all(is.character(mask))) mask = raster::raster(mask)

  mixfit = parallel::mclapply(as.list(image), function(img){
    x = img[ mask!=0]
    nx = round(length(x) * subsamp)
    if(nx>maxpix) nx = maxpix
    x = sample(x, size=nx, replace=FALSE)
    # this doesn't work real well
    mixfit = mixtools::gammamixEM(x, lambda=lambda, alpha=alpha, beta=beta, k=k)
    }, mc.cores = mc.cores)

}
