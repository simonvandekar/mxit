#' Calculates density of image intensity.
#'
#' @param image Required character vector of image locations or raster object.
#' @param mask Required character vector of image location or raster object of the tissue mask.
#' @param bias_corrected logical to mark whether the density has been bias corrected
#' @param minpix Minimum number of pixels in `mask` to calculate density in `image`
#' @param fact NOT SURE -- check w Simon
#' @param logIt logical to return log density if TRUE
#' @keywords density of image intensity
#' @return Returns the density values of image intensity for the `image`
#' @importFrom raster raster
#' @importFrom spatstat cellStats
#' @importFrom stringr str_extract
#' @export
imageDensity = function(img,
                        maskimg,
                        bias_corrected=FALSE,
                        minpix=100,
                        fact=NULL,
                        logIt=FALSE){
  img_name = str_extract(img,"[ \\w-]+?(?=\\.)")
  if(is.na(maskimg)){ ## quit if mask image is null
    NA
  } else {
    mask = raster(maskimg)
    if(cellStats(mask, stat='sum')<minpix){ ## quit if mask image doesn't have enough data
      NA
    } else {
      img = raster(as.character(img))
      if(!is.null(fact)){ ## not entirely sure what this does
        mask = aggregate(mask, fact=fact)
        img = aggregate(img, fact=fact)
      }
      ## --- BREAK POINT --- (e.g. high computation time)
      vals = img[ mask!=0] ## subset image values only where mask exists
      if(logIt){
        out = density(log(vals/mean(vals))) #calculate log density
      } else{
        out = density(vals/mean(vals)) ## calculate reg density
      }
      attr(out,"name") = img_name
      attr(out,"bias_correct") = bias_corrected
      out ## add attribute for file name?
    }
  }
}
