# how to look energy flow and CO2 cycle from satellite: faPAR fraction of absorbed photosynthetically active radiation

library(raster)                                           # install.packages("raster")
library(rasterVis)                                        # for levelplot
library(rasterdiv)

                                                          #copDVI= copernicus DIV = World wide analysis DVI different variation index
setwd("C:/lab/")
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA))        #reclissifaty coppenNDVI   remuving data from 255-255 and putting No Data= NA
levelplot(copNDVI)                                        #levelplot cran
                                                          #different vegetation index... you can make an avarage of the pixel variable value, very high NDVI are in primary equator forest and broad-laef and conifer northern forest

faPAR10 <- raster("faPAR10.tif")                          #give name 10 per pixel aggregation 10=1 ,   raster = data   faPAR10= name of image
levelplot(faPAR10)

pdf("copNDVI.pdf")                                        #make a pdf for DVI and faPAR
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
library(raster)
library(rasterdiv)
library(rasterVis)

ls() #lookin for faPAR10
faPAR10

writeRaster(copNDVI, "copNDVI.tif")                         # write the data copNDVI in .tif, 5.3MB

faPAR <- stretch(faPAR10,minv=0,maxv=250)                   # faPAR10 in in bits from 0 to 0.93. Change from 0 to 255

writeRaster(faPAR, "faPAR.tif")
faPAR

                                                            # EX. make a levelplot for this set


 
