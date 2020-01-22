#' Utilizes biasCorrect.R to experiment and work with IHC.
#'
#'
#'
#'

## ---- DATA

## access to all images
# path.all = "~/Box/multiplex/MxIF\ Data/H051_Ileum_H_1/AFRemoved"
# imgs.all = list.files(path.all,pattern="tif$",full.names = TRUE,recursive = TRUE)

## local access to test data
path.imgs = "~/Documents/Current/Research/mxit/test-dat"
imgs.test = list.files(path.imgs,pattern="tif$",full.names = TRUE,recursive = TRUE)

## local access to test epi masks
path.masks = "~/Documents/Current/Research/mxit/test-masks"
masks.test = list.files(path.masks,pattern="png$",full.names = TRUE,recursive = TRUE)


## ----- PACKAGES
require(parallel)
require(ANTsR)
require(raster)
require(ggplot2)
require(stringr)
require(viridis)
require(spatstat)

## ----- IMPORT FUNCTIONS
setwd("~/Documents/Current/Research/mxit/R")

source("biasCorrect.R")
source("imageDensity.R")
source("plotDensity.R")

## ----- ALL TEST IMAGES

imgs = list() ## list to hold all results

## loop through imgs.test
for(i in 1:length(imgs.test)){
  ## bias correct all images
  bc_img = biasCorrect(imgs.test[i],
                       masks.test[i])
  
  ## calculate image density on original images
  orig_dens = imageDensity(imgs.test[i],
                           masks.test[i])
  
  ## calculate image density on bias-corrected images
  bc_dens = imageDensity(bc_img,
                           masks.test[i], #adjusted mask?
                           bias_corrected = TRUE)
  
  imgs = list(imgs,
              orig_dens,
              bc_dens)
}

## plot density comparison
plotDensity(imgs)

### --- WORKING TEST EXAMPLE (NO BIAS CORRECTION)

i1 = imgs.test[1]; m1 = masks.test[1]
i2 = imgs.test[2]; m2 = masks.test[2]

d1 = imageDensity(i1,m1); d2 = imageDensity(i2,m2)
attr(d2,"name") <- attr(d1,"name")
attr(d2,"bias_correct") <- TRUE

imgs = list(d1,d2) ## list to hold all results

## plot density comparison
plotDensity(imgs)

### --- bc working ex - not working ATM

biasCorrect(imgs.test[1],
            masks.test[1])
## ABORTS R

i1 = imgs.test[1]
m1 = masks.test[1]
splineParam=2000

image = raster::stack(i1)
mask = raster::raster(m1)
out = n4BiasFieldCorrection(img=as.antsImage(as.matrix(image)), 
                            mask=as.antsImage(as.matrix(mask)), 
                            splineParam = splineParam)
## ERROR
## Requested region is (at least partially) outside the largest possible region.