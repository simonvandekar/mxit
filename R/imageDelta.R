#' Determines differences between two images and returns raster of difference.
#'
#' @param loc1 location of first image
#' @param loc2 location of second image
#' @keywords difference between two raster images
#' @return Returns the raster image of difference between `loc1` and `loc2`
#' @importFrom raster raster
#' @importFrom stringr str_extract
#' @export
imageDelta <- function(loc1,loc2){
  img1 = raster(as.matrix(raster(loc1)))
  img2 = raster(as.matrix(raster(loc2)))
  
  img.delta = img2 - img1
  attr(img.delta, "name") = str_extract(loc1,"[ \\w-]+?(?=\\.)")
  img.delta
}
