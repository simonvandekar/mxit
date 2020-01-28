#' Utilizes biasCorrect.R to experiment and work with IHC.

## ---- DATA

## use other images (original)

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
source("imageDelta.R")

image = imgs.test[i]
mask = masks.test[i]

dens = imageDensity(image,mask)
p = plotDensity(dens)
p

## ----- ALL TEST IMAGES

imgs = list() ## list to hold all results
img_deltas = list()

## loop through imgs.test
for(i in 1:length(imgs.test)){
  ## bias correct all images
  
  # ## calculate image density on original images
  # orig_dens = imageDensity(imgs.test[i],
  #                          masks.test[i])
  # 
  # ## calculate image density on bias-corrected images
  # bc_dens = imageDensity(bc_img,
  #                          masks.test[i],
  #                          bias_corrected = TRUE)
  
  # imgs = list(imgs,
  #             orig_dens,
  #             bc_dens)
  
  #find deltas
  img_bc = biasCorrect(imgs.test[i],
                       mask,
                       simple=TRUE,
                       splineParam = 2000)
  
  img_og = raster(as.matrix(raster(imgs.test[i])))
  name_var = stringr::str_extract(imgs.test[i],"[ \\w-]+?(?=\\.)")
  
  img_del = img_og - img_bc
  writeRaster(img_bc, paste0("../test-results/",name_var,"_bias_correct.tif"))
  writeRaster(img_del, paste0("../test-results/",name_var,"_delta.tif"))
  img_deltas = list(img_deltas,
                    img_del)
}

## plot density comparison
plotDensity(imgs)

### --- WORKING TEST EXAMPLE (NO BIAS CORRECTION)

i1 = imgs.test[1]
m1 = masks.test[1]
i2 = biasCorrect(i1,m1)

d1 = imageDensity(i1,m1); d2 = imageDensity(i2,m1)
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


out = n4BiasFieldCorrection(img=imgs.test[1], 
                            mask=masks.test[1], 
                            splineParam = splineParam)
out = as.matrix(out)
out = out/max(out)
out = t(out)
tiff::writeTIFF(out,"b1p.tif")

out = as.matrix(out)
out = t(out)
out = round(out)
out1 = raster::raster(out)
writeRaster(out1, "raster3.tif")
rm(out); rm(out1)


## delta image
r_bc = raster(as.matrix(raster("raster3.tif")))
r_og = raster(as.matrix(raster(imgs.test[1])))
r_del = r_og - r_bc

## MIXTURE MODELS
## log transform -> gaussianmmixture model
## untransformed -> gamma mixture model
## ERROR
## Requested region is (at least partially) outside the largest possible region.

##testing biasCorrect

# a1 = biasCorrect(imgs.test[1],
#                  masks.test[1])
# a2 = biasCorrect(list(raster(imgs.test[1])),
#                  raster(masks.test[1]))
# # b1 = biasCorrect(imgs.test[1:3],
# #                  masks.test[1])
# # b2 = biasCorrect(lapply(imgs.test[1:3], function(img) raster(img)),
# #                  masks.test[1])