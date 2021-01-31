# R_code_ecosystem_functions.r

# R code to view biomass over the world and calculate changes in ecosystem functions
# energy
# chemical cycling
# proxies

install.packages("rasterdiv")                                                          # packages to assess Rao variations in different ecosystems
library(rasterdiv)
install.packages("rasterVis")                                                          # raster visualisation
library (rasterVis)                                                                    ## Loading required package:  raster 
library(ratser)                                                                        ## Loading required package:  sp
                                                                                       ## Loading required package:  lattice
                                                                                       ## Loading required package:  latticeExtra

data(copNDVI)
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253, 255, NA), right=TRUE)                        #removing water pixels using cbind argument for 253, 255
levelplot(copNDVI)                                                                     #plot the same image for factors=10 and factors=100

copNDVI10 <- aggregate(copNDVI, fact=10)
levelplot(copNDVI10)

copNDVI100 <- aggregate(copNDVI, fact=100)                                             # see: High resolution satellite imagery for tropical biodiversity studies: The devil is in the detail
levelplot(copNDVI100)

library(ggplot2)

myPalette <- colorRampPalette(c('white','green','dark green'))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))                 # to plot NDVI score on the map, very impressive map

ggR(copNDVI, geom_raster = TRUE)+
scale_fill_gradientn(name = "NDVI", colours = myPalette(100))+
labs(x="Longitude",y="Latitude", fill="")+
#   theme(legend.position = "bottom") +
  NULL
# +
# ggtitle("NDVI")

                                                                                        ## deforestation project
setwd("D:/lab/")
defor1 <- brick("defor1_.jpg")                                                          #to import images and link them with bands
defor2 <- brick("defor2_.jpg")                                                          # band1: NIR, defor1_.1# band2: red, defor1_.2# band3: green

plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")                                           #par function plot to compare
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")


par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

                                                                                        #classify the DVI for both images in different class

dvi1 <- defor1$defor1_.1 - defor1$defor1_.2
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2                                             # to plot NDVI score on the map, very impressive map
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100)
par(mfrow=c(2,2))
plot(dvi1, col=cl)
plot(dvi2, col=cl)

difdvi <- dvi1 - dvi2                                                                   ## Warning in dvi1 - dvi2:  Raster objects have different extents.  Result for their intersectionis returned
cld <- colorRampPalette(c('blue','white','red'))(100) 

plot(difdvi, col=cld)                                                                   # to see the loss of ecosystem services
hist(difdvi)                                                                            # high loss in primary production ecosystem services: hist=histogram

sessionInfo ()                                                                          # R session information, matrix, packages
