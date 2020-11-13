#' Computes a locally weighted correlation
#'
#' This function computes a locally weight correlation between imaging channels
#' @param dist A list of distances from the center cell.
#' @param data Variables to compute the correlation between.
#' @param fwhm The full-width at half maximum for the weights.
#' @keywords weighted regression
#' @return Returns a list with the off-diagonal weighted correlation between the channels
#' @importFrom stats cor
#' @export
lwc = function(dist, data, fwhm){
  sd = fwhm/2/sqrt(2*log(2))
  weights = dnorm(c(0, dist), sd=sd )
  weights = weights/sum(weights)
  weightedData = as.matrix(sweep(data, 2, colSums(data*weights), FUN='-')) * sqrt(weights)
  res = cov2cor(crossprod(weightedData))
  nms = outer(rownames(res), colnames(res), paste, sep='_WITH_')[upper.tri(res, diag=FALSE)]
  res = res[upper.tri(res, diag=FALSE)]
  names(res) = nms
  res
}

#' MSE of locally weighted correlation versus cormat
#'
#' This function computes a locally weight correlation between imaging channels
#' @param dist A list of distances from the center cell.
#' @param data Variables to compute the correlation between.
#' @param fwhm The full-width at half maximum for the weights.
#' @param cormat The reference correlation matrix.
#' @keywords weighted regression
#' @return Returns a list with the off-diagonal weighted correlation between the channels
#' @importFrom stats cor
#' @export
corMSEc = function(dist, data, fwhm, cormat){
  sd = fwhm/2/sqrt(2*log(2))
  weights = dnorm(c(0, dist), sd=sd )
  weights = weights/sum(weights)
  weightedData = as.matrix(sweep(data, 2, colSums(data*weights), FUN='-')) * sqrt(weights)
  res = cov2cor(crossprod(weightedData))
  res = mean((res-cormat)[lower.tri(res, diag=FALSE)]^2)
  names(res) = 'corMSE'
  res
}

#' Computes a locally weighted PCA
#'
#' This function computes a locally weight correlation between imaging channels
#' @param dist A list of distances from the center cell.
#' @param data Variables to compute the correlation between.
#' @param fwhm The full-width at half maximum for the weights.
#' @param ncomp Number of components to return
#' @keywords weighted regression
#' @return Returns a list with the off-diagonal weighted correlation between the channels
#' @importFrom stats cor
#' @export
pcac = function(dist, data, fwhm, ncomp=5){
  sd = fwhm/2/sqrt(2*log(2))
  weights = dnorm(c(0, dist), sd=sd )
  weights = weights/sum(weights)
  weightedData = as.matrix(sweep(data, 2, colSums(data*weights), FUN='-')) * sqrt(weights)
  ssqrt = sqrt(colMeans(weightedData^2))
  res = svd(scale(weightedData, center=FALSE, scale=ifelse(ssqrt==0, 1, ssqrt)), nu = 0, nv=0)$d^2
  # proportion of variance for each component
  res = res[1:ncomp]/sum(res)
  names(res) = paste('pc', 1:ncomp, sep='')
  res
}

#' Computes Coupling matrix
#'
#' This function computes the coupling metric which summarizes intracellular correlations within local tissue 
#' @param cellID A vector of cellIDs.
#' @param coords A vector of coordinates of cell centroids.
#' @param neighbors A list of vectors of cellIDs which are neighbors for each element of cellID.
#' @param dist A list of vectors of distances between the elements in neighbors and the cell identified by cellID.
#' @param fwhm The full-width at half maximum for the weights.
#' @param data Variables to compute the correlation between.
#' @param k Number of nearest neighbors for knn algorithm.
#' @param subset Subset of CellIDs to extract for viewing local relationships.
#' @param coup A coupling function that takes at least distance, data, and a FWHM arguments and returns an object
#' @param ... Arguments passed to coup function
#' @keywords weighted regression
#' @return Returns a list with the off-diagonal weighted correlation between
#' @importFrom raster raster stack
#' @importFrom FNN get.knn
#' @importFrom parallel mclapply
#' @export
coupling = function(cellID, coords=NULL, neighbors=NULL, dist=NULL, fwhm, data, k=10, subset=NULL, coup=lwc, ...){
  if(is.null(neighbors) | is.null(dist)){
    if(is.null(coords)){
      stop('must supply coords or neighbors and dist.')
    }
    nd = FNN::get.knn(coords, k=k)
    neighbors = apply(nd$nn.index, 1, function(neighs) cellID[neighs])
    # convert to list
    neighbors = lapply(seq_len(ncol(neighbors)), function(i) neighbors[,i])
    dist = lapply(seq_len(nrow(nd$nn.dist)), function(i) nd$nn.dist[i,])
    rm(nd)
  }
  if(!is.null(subset)){
    subsetLogical=TRUE
    subset = which(cellID %in% subset)
  } else {
    subset = 1:length(cellID)
    subsetLogical=FALSE
  }
  res = parallel::mclapply(subset, function(ind) coup(dist[[ind]], data=data[c(ind, which(cellID %in% neighbors[[ind]]) ),], fwhm=fwhm, ...) )
  res = data.frame(cellID=cellID[subset], do.call(rbind, res))
  if(subsetLogical){
    res = list(cors=res, data=lapply(subset, function(ind){out = data[c(ind, which(cellID %in% neighbors[[ind]]) ),]
                                                        out$dist = c(0, dist[[ind]])
                                                        out})
    )
  } else {
    res = list(cors = res, data = NULL)
  }
}


plotdata = function(result, cellIDs=NULL){
  corvars = unique(do.call(c, strsplit(names(result$cors)[-1], split = '_WITH_')))
  cormat = matrix(NA, nrow=length(corvars), ncol=length(vars))
  # heatmap
  sapply(cellIDs, function(cellID){
    cormat[upper.tri(cormat, diag = FALSE)] = as.numeric(result$cors[ result$cors$cellID==cellID, -1])
    diag(cormat) = 1
    cormat[lower.tri(cormat, diag = FALSE)] = as.numeric(result$cors[ result$cors$cellID==cellID, -1])
    rownames(cormat) = colnames(cormat) = corvars # gsub('Median_Cell_', '', corvars)
    col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                   "#4393C3", "#2166AC", "#053061")))(200)
    corrplot::corrplot(cormat, type = "upper", tl.col = "black", method='number', diag=FALSE, col=col2, number.cex=0.8)
    plot(result$data[[which(result$cors$cellID==cellID) ]][, -which(names(result$data[[1]]) %in% 'dist')], col=gray(1-result$data[[which(result$cors$cellID==cellID)]]$dist/max(result$data[[which(result$cors$cellID==cellID)]]$dist)) )
  })
}
