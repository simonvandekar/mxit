#' Performs n4 inhomogeneity correction on an image.
#'
#' Adjusts for slow frequency fluctuations in image intensities, i.e. field bias in MRI.
#' This function computes a k component gamma mixture model marginally for each channel and returns a mixfit object to perform image normalization using the XX function.
#' @param image Required character of image location or raster object.
#' @param mask Optional character of image location or raster object of the tissue mask.
#' @param splineParam Controls smoothness of bias field. Argument passed to ANTsRCore::n4BiasFieldCorrection.
#' @param mc.cores Number of cores to use for parallel things.
#' @param ... Other arguments passed to ANTsRCore::n4BiasFieldCorrection.
#' @keywords bias correction
#' @return Returns the bias corrected images as a raster stack object.
#' @importFrom raster raster t
#' @importFrom ANTsRCore as.antsImage n4BiasFieldCorrection
#' @export
n4 = function(image, mask, splineParam=2000, mc.cores=getOption("mc.cores", 2L), ...){
  if(is.character(image)){
    out=raster(image)
    if(missing(mask)){
      out[,] = as.matrix(n4BiasFieldCorrection(img=image, splineParam = splineParam, ...))
    } else {
      out[,] = as.matrix(n4BiasFieldCorrection(img=image, mask=mask, splineParam = splineParam, ...))
    }
  } else if(is.raster(image)){
    out = image
    if(missing(mask)){
      out[,] = as.matrix(n4BiasFieldCorrection(img=as.antsImage(as.matrix(image)), splineParam = splineParam, ...))
    } else {
      out[,] = as.matrix(n4BiasFieldCorrection(img=as.antsImage(as.matrix(image)), mask=mask, splineParam = splineParam, ...))
    }
  }
  return(raster::t(out))
}





#' Performs inhomogeneity correction on an image using a penalized GAM.
#'
#' Adjusts for slow frequency fluctuations in image intensities, i.e. field bias in MRI.
#' This function computes a k component gamma mixture model marginally for each channel and returns a mixfit object to perform image normalization using the XX function.
#' @param image Optional character of image location or raster object.
#' @param mask Optional character of image location or raster object of the tissue mask.
#' @param df degrees of freedom for the basis. This is the total degrees of freedom for the entire image, not each dimension's degrees of freedom.
#' @param ... named arguments passed to bam.
#' @keywords bias correction
#' @return Returns the bias corrected image.
#' @importFrom raster raster t
#' @importFrom mgcv bam
#' @export
sbc.bam = function(image, mask, df=16, seq.by=3, ...){
    if(is.character(image)){ image = raster(image) }
    dims = dim(image)
    image = log10(image + 1)
    if(missing(mask)){
      # sample data
      x = as.data.frame(lapply(expand.grid(row=seq(1, dims[1], by = seq.by), col=seq(1, dims[2], by = seq.by)), as.numeric))
    } else {
      if(is.character(mask)){ mask = raster(mask) }
      x = as.data.frame(which(mask>0, arr.ind = TRUE))
      image[ which(mask<=0) ] = 0
    }
    # fit the model
    x$y = image[ as.matrix(x)]
    freml = bam(y ~ s(row, col, k=df), data=x, ...)


    f = image
    if(missing(mask)){
      fits = c(predict(freml, newdata=as.data.frame(lapply(expand.grid(row=seq(1, dims[1]), col=seq(1, dims[2])), as.numeric))) -  freml$coefficients[1])
      f[,] = fits
      f = raster::t(f)
    } else {
      fits = fitted(freml) - freml$coefficients[1]
      f[ as.matrix(x[,c('row', 'col')]) ] = fits
    }
    image = image - f
    result = list(uncorrupted=(10^(image)-1), field=10^f)
}

#' Performs inhomogeneity correction on an image using a linear model unpenalized cubic splines.
#'
#' Adjusts for slow frequency fluctuations in image intensities, i.e. field bias in MRI.
#' This function computes a k component gamma mixture model marginally for each channel and returns a mixfit object to perform image normalization using the XX function.
#' @param image Optional character of image location or raster object.
#' @param mask Optional character of image location or raster object of the tissue mask.
#' @param df degrees of freedom for the basis. This is each dimension's degrees of freedom.
#' @keywords bias correction
#' @return Returns the bias corrected image.
#' @importFrom raster raster t
#' @importFrom splines ns
#' @export
sbc.lm = function(image, mask, df=4, seq.by=3, ...){
    if(is.character(image)){ image = raster(image) }
    dims = dim(image)
    image = log10(image + 1)
    if(missing(mask)){
      # sample data
      x = as.data.frame(lapply(expand.grid(row=seq(1, dims[1], by = seq.by), col=seq(1, dims[2], by = seq.by)), as.numeric))
    } else {

      if(is.character(mask)){ mask = raster(mask) }
      x = as.data.frame(which(mask>0, arr.ind = TRUE))
      image[ which(mask<=0) ] = 0
    }
    # fit the model
    x$y = image[ as.matrix(x)]
    ksr = round(seq(1, dims[1], length.out=df+1)[c(2:(df))])
    ksc = round(seq(1, dims[2], length.out=df+1)[c(2:(df))])
    model = lm(y ~ ns(row, knots=ksr)*ns(col, knots=ksc), data=x)


    f = image
    if(missing(mask)){
      fits = c(predict(model, newdata=as.data.frame(lapply(expand.grid(row=seq(1, dims[1]), col=seq(1, dims[2])), as.numeric))) -  model$coefficients[1])
      f[,] = fits
      f = raster::t(f)
    } else {
      fits = fitted(model) - model$coefficients[1]
      f[ as.matrix(x[,c('row', 'col')]) ] = fits
    }
    image = image - f
    result = list(uncorrupted=(10^(image)-1), field=10^f)
}

