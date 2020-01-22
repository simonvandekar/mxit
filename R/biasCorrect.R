#' Performs inhomogeneity correction (IHC) on an array of images.
#'
#' Adjusts for slow frequency fluctuations in image intensities, i.e. field bias in MRI.
#' This function computes a k component gamma mixture model marginally for each channel and returns a mixfit object to perform image normalization using the XX function.
#' @param image Required character vector of image locations or raster object.
#' @param mask Required character vector of image location or raster object of the tissue mask.
#' @param splineParam Controls smoothness of bias field. Argument passed to ANTsRCore::n4BiasFieldCorrection.
#' @param mc.cores Number of cores to use for parallel things.
#' @param ... Other arguments passed to ANTsRCore::n4BiasFieldCorrection.
#' @keywords bias correction
#' @return Returns the bias corrected images as a raster stack object.
#' @importFrom raster as.raster raster stack
#' @importFrom ANTsRCore as.antsImage n4BiasFieldCorrection
#' @importFrom parallel mclapply
#' @export
biasCorrect = function(image, mask, splineParam=2000, mc.cores=getOption("mc.cores", 2L), ...){
  if(all(is.character(image))) image = raster::stack(image)
  # will error if multiple masks are passed
  if(!missing(mask) & all(is.character(mask))) mask = raster::raster(mask)
  # out = mclapply(as.list(image), function(img){
  #   n4BiasFieldCorrection(img=as.antsImage(as.matrix(img)), 
  #                         mask=as.antsImage(as.matrix(mask)), 
  #                         splineParam = splineParam, ...)
  # }, mc.cores=mc.cores)
  
  out = lapply(as.list(image), function(img){
    n4BiasFieldCorrection(img=as.antsImage(as.matrix(img)), 
                          mask=as.antsImage(as.matrix(mask)), 
                          splineParam = splineParam)
  })
  lapply(out, as.raster)
  stack(out)
}
