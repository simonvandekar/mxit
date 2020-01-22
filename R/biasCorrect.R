#' Performs inhomogeneity correction (IHC) on an array of images.
#'
#' Adjusts for slow frequency fluctuations in image intensities, i.e. field bias in MRI.
#' This function computes a k component gamma mixture model marginally for each channel and returns a mixfit object to perform image normalization using the XX function.
#' @param image Required list of image locations or raster objects.
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
biasCorrect = function(image, 
                       mask, 
                       splineParam=2000, 
                       mc.cores=getOption("mc.cores", 2L), 
                       simple=FALSE,
                       ...){
  if(simple){
    out = n4BiasFieldCorrection(img=image,
                                mask=mask,
                                splineParam = splineParam)
    return(raster(t(as.matrix(out))))
  }
  
  ## if not file location strings, setup ANTS images
  if(!all(is.character(image[[1]]))){
    image = lapply(image, function(img){
      as.antsImage(as.matrix(img))
    })
  }
  
  if(!missing(mask) & #will error if multiple masks are passed
     !all(is.character(mask))){
    mask = as.antsImage(as.matrix(mask))
  }
  
  ## apply bias correction
  out = mclapply(image, function(img){
    n4BiasFieldCorrection(img=img,
                          mask=mask,
                          splineParam = splineParam)
  }, mc.cores=mc.cores)
  
  #memory cleaning - prevents abort?
  out1 <- out
  rm(out)
  
  ## change format to return as list of rasters in correct orientation
  out1 = lapply(out1, function(i){
    raster(t(as.matrix(i)))
  })
  out1
}
