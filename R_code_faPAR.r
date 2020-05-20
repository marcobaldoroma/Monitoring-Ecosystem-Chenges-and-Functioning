# how to look energy flow and CO2 cycle from satellite: faPAR fraction of absorbed photosynthetically active radiation

library(raster)                                             # install.packages("raster")
library(rasterVis)                                          # for levelplot
library(rasterdiv)

                                                            #copDVI= copernicus DIV = World wide analysis DVI different variation index
setwd("C:/lab/")
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA))          #reclissifaty coppenNDVI   remuving data from 255-255 and putting No Data= NA
levelplot(copNDVI)                                          #levelplot cran
                                                            #different vegetation index... you can make an avarage of the pixel variable value, very high NDVI are in primary equator forest and broad-laef and conifer northern forest

faPAR10 <- raster("faPAR10.tif")                            #give name 10 per pixel aggregation 10=1 ,   raster = data   faPAR10= name of image
levelplot(faPAR10)

pdf("copNDVI.pdf")                                          #make a pdf for DVI and faPAR
levelplot(copNDVI)
dev.off()

pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()

######################################### second part 

setwd("C:/lab/")

load("faPAR.RData")
                                                             # the original faPAR from copernicus is 2gb
                                                             # let's see how much space needs for 8-bits image

ls()                                                         #lookin for faPAR10
faPAR10

writeRaster(copNDVI, "copNDVI.tif")                          # write the data copNDVI in .tif, 5.3MB

faPAR <- stretch(faPAR10,minv=0,maxv=250)                    # faPAR10 in in bits from 0 to 0.93. Change from 0 to 255

writeRaster (faPAR, "faPAR.tif")
faPAR

levelplot(faPAR)                                             # EX. make a levelplot for this set



#### third part
                                                             ## regression model between faPAR and NDVI
                                                             ## make before a general example with erosion and heavy metal.

erosion <- c(12,14,16,24,26, 40,55,67)

                                                             # heavy metal in ppm
hm <- c(30, 100, 150, 200, 260, 340, 460, 600)

plot(erosion, hm, col="red", pch=19, xlab = "Erosion", ylab= "Heavy Metal")   
                                                             #point character is just the symble pch=19 
#linear model function
model1 <- lm(hm ~ erosion)
summary(model1)
                                                             # the line describe by the ab ( slope and intercept)
abline(model1)

                                                             # make the linear model for the relationship between biomass and ecosystem function # we want the same analysis for faPAR and NDVI
library (raster)
library(rasterdiv)

setwd("C:/lab/")

faPAR10 <- raster("faPAR10.tif")
faPAR10                                                      # 56M of points/cells
plot(faPAR10)
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)

install.packages ("sf")
library (sf)                                                 # to call st_* functions 
                                                             # random point we will cal our new function
random.points <- function(raster,n)                          # algorithm to select the point on our raster image randomly
{
lin <- rasterToContour(is.na(raster))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}

                                                             # then select 1000 points randomly from faPAR10 rasters
pts <- random.points(faPAR10, 1000)
plot(faPAR10$pts)


copNDVIp <- extract(copNDVI, pts)
faPAR10p <- extract(faPAR10, pts)
                                                             # photosynthesis vs biomass
model2 <- lm(faPAR10p ~ copNDVIp)
summary(model2)                                              # having a look of the correlation of the two variables

plot(copNDVIp, faPAR10p, col="green", xlab="biomass", ylab="photosynthesis")
abline(model2, col="red")

levelplot(copNDVI)
levelplot(faPAR10)

