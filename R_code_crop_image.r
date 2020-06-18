




library(raster)
library(ncdf4)   #the library to read our pixels

snow <- raster( "c_gls_SCE500_202005180000_CEURO_MODIS_V1.0.1.nc")

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)        #colorRamp with the meanest colors

plot(snow, )
# making a crop of our images i.g. italy 

ext <- c(0, 20, 35, 50)   # this is the extent in which we want to zoom on.
   # we have to say 1 the image and second the extent 

 zoom(snow, ext=ext)
 
but we can also cutting the image that is calls crop the function 

# in this case the extent is direct related at the previus
# this is very useful to use the area of interest

snowitaly <- crop(snow, ext)
plot(snowitaly, col=cl)

# you can drow one specific area in other way. drowing the extent
# to give a geometry for example a rectangle
# () this command to say we are working in the image as a vector 

zoom(snow, ext=drawExtent())  # drop and zoom can be done from exent and with draw function


setwd("~/lab")
library(raster)
library(ncdf4) 

snow <- raster("c_gls_SCE500_202005180000_CEURO_MODIS_V1.0.1.nc")

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow, col=cl)

ext <- c(0, 20, 35, 50)
zoom(snow, ext=ext)

# crop and create a new image
snowitaly <- crop(snow, ext) 
plot(snowitaly, col=cl) 

# zoom(snow, ext=drawExtent())
# crop(snow, drawExtent())


Interpolation of the data thanks spatstat library and species distribution modelling to understand the position and structure of the forest. 
Another nice analysis is about diameter and height of the forest.


## Interpolation: spatstat library
# library(dbmss)
library(spatstat)

inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T)


attach(inp)
plot(X,Y)
inppp <- ppp(x=X,y=Y,c(716000,718000),c(4859000,4861000))

marks(inppp) <- Canopy.cov
canopy <- Smooth(inppp)
plot(canopy)

marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp)

par(mfrow=c(1,2))
plot(s)
points(inppp)
plot(lichs)
points(inppp)

#########

# Dati psammofile Giacomo
inp.psam <- read.table("Dataset_psammofile_giacomo.csv", sep=";", head=T)

attach(inp.psam)

summary(inp.psam)

plot(E,N)
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))

marks(inp.psam.ppp) <- C_org
C <- Smooth(inp.psam.ppp)
plot(C)
points(inp.psam.ppp)
