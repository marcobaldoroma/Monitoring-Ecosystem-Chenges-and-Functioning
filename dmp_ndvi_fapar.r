library(GGally)
library(rgdal)
library(gdalUtils) #to transform the file format
library(raster)
library(RStoolbox)
library(ncdf4)
#install.packages("sf")
library(sf) 

setwd ("D:/toy/")

dmp <- raster ("rasterDMP.tif")
ndvi <- raster ("rasterNDVI.tif")
fapar <- raster ("rasterFAPAR.tif")

cl <- colorRampPalette(c('white', 'pink', 'red', 'dark green'))

par(mfrow=c(3,1))
plot(dmp, col=cl, main= "Dry Matter Productivity 300")   # plot(d, col=cl, main="Densities of covid-19") 
plot(ndvi, col=cl, main="NDVI300")
plot(fapar, col=cl, main="FAPAR300")

#plot(dmp, col=cl, main= "Dry Matter Productivity 300")
#plot(ndvi,col=cl, main="NDVI300")
#plot(fapar, col=cl, main="FAPAR300")
# crop my AOI
ext <- c(75.10, 75.60, 12.05, 12.55)

dmp.up <- crop(dmp,ext)
ndvi.up <- crop(ndvi,ext)
fapar.up <- crop(fapar,ext)

par(mfrow=c(3,1))

plot(dmp.up, col=cl, main= "Dry Matter Productivity 300")
plot(ndvi.up,col=cl, main="NDVI300")
plot(fapar.up, col=cl, main="FAPAR300")

random.points <- function(raster,n)                          
{
lin <- rasterToContour(is.na(raster))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') 
pts <- spsample(pol[1,], n, type = 'random')
}

# then select 1000 points randomly from faPAR10 rasters

pts <- random.points(dmp.up, 1000)
plot(dmp.up$pts)


############# new part
DMPp <- extract(dmp.up, pts)
NDVIp <- extract(ndvi.up, pts) # extract 1000 random points from copNDVI
faPARp <- extract(fapar.up,pts)
#build the linear model between copNDVIp and faPAR10 (copNDVIp because the calculation is faster with less values)
# the line is calculated by reducing the distance between the points (x;y) in the graph
# phothosythesis vs biomass
model2 <- lm(dmp.up ~ ndvi.up) # R = 0.4 because in conifer forest biomass is high and phothosynthesis but not that high. p = 2 ^-16 two variables are related each other
plot(faPAR10p, NDVIp)
abline(model2, col="red")

