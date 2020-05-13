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

setwd("C./lab/")

load("faPAR.RData")
                                                             # the original faPAR from copernicus is 2gb
                                                             # let's see how much space is for 8bits
library(raster)
library(rasterdiv)
library(rasterVis)

writeRaster(copNDVI, "copNDVI.tif")
levelplot (writeRaster)                                      #5.3MB
                                                             # EX.levelplot this set
