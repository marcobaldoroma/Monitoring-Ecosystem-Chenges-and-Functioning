# R_code_snow.r

# setwd("~/lab/") #linux
# setwd("/Users/utente/lab") #mac
setwd("C:/lab/") 

install.packages("ncdf4")                                             # new library to read a different format data files
library(ncdf4)
library(raster)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.NC")     #giving the name for visualization process, we get 
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)        # we make a colorramppalette for snow changes blue to white

# Exercise: plot snow cover with the cl palette
plot(snowmay,col=cl)  

                                                                      ##### import snow data # for several layers all together in one time 
                                                                      # he aggregates and clumped several pixel to make the images less heavy
                                                                      # you can create a new folder i.e. Snow
                                                                      # Slow manner to import set
                                                                      # new working directory setwd("C:/lab/snow")
# setwd("~/lab/snow") #linux
# setwd("/Users/utente/lab/snow") #mac

setwd("C:/lab/snow") # windows

snow2000 <- raster("snow2000r.tif")                                   # for example snow in 2000 and with raster you can select the layer of interest, and you can import all together
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))                                                     # we want plot all images together 
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

                                                                      # import the set in faster way
                                                                      # lapply you can repet the same fuction for the whole set
                                                                      # you can list files of the all folder. i.e. snow_2000, snow_2005, ..... snow_2020

rlist <- list.files(pattern="snow")                                   # lapply function repet the raster fuction to import the interest layer of the all set snow
rlist

import <- lapply(rlist, raster)   
                                                                      # we import the layer for all the single files with the function raster fuction repet for the whole set thanks at the lapply fuction
                                                                      # this way of work is call stack or raster stack
                                                                      # stack of all the dataset that we import 
snow.multitemp <- stack(import)                                       # we give the proper name at the vector (snow multitemp, because we are analysing snow cover in a temporal scale
 
plot(snow.multitemp, col=cl)                                          # it is important to make several analysis all together at the same time with much powerful information range and in faster way.

#########################################################################################################################################
# make a prediction of snow cover, making a graph with time(x) and snow cover %(y)
# with a regression functions we can fit our set to predict snow cover for the next 5 years using the trend of the last 20 years
# Big Data structure of complex data, we want see how make the importation of data.
# lets look at the function "prediction"

## let make a prediction of snow cover with a simple line code

source("prediction.r")                                                # with source function you can select a file like a R script and it will work by themself for all my analysis.
                                 
############################################################ what there are inside the prediction function that we created?
#require(raster)
#require(rgdal)

# define the extent
#ext <- c(-180, 180, -90, 90)
#extension <- crop(snow.multitemp, ext)
    
# make a time variable (to be used in regression)
# time <- 1:nlayers(snow.multitemp)

# run the regression
#fun <- function(x) {if (is.na(x[1])){ NA } else {lm(x ~ time)$coefficients[2] }} 
#predicted.snow.2025 <- calc(extension, fun) # time consuming: make a pause!
#predicted.snow.2025.norm <- predicted.snow.2025*255/53.90828

##################### day second

setwd("C:/lab/snow/")                                                  # set working directory with snow folder

library(raster)                                                        # to make our work with lapply function

# Exercise: import all of the snow cover images all together (snow cover elaborated layer)

rlist <- list.files(pattern="snow")
rlist

import <- lapply(rlist, raster)                                        # https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/lapply
                                                                       # we import the layer for all the single files with the raster fuction repet for the whole set thanks at the lapply fuction
                                                                       # this way of work is call stack or raster stack
                                                                       # stack of all the dataset that we import 
snow.multitemp <- stack(import)                                        # we give the proper name at the vector (snow multitemp, because we are analysing snow cover in a temporal scale
                                                                       # stack do a job similar at par, but faster
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)         # we make a colorramppalette for snow changes blue to white
plot(snow.multitemp, col=cl)                                           # it is important to make several analysis all together at the same time with much powerful information range and in faster way.


load("R_code_snow.r.RData")      # load preview work/script

prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)

# the idea is to make a regretion model to understand the future snow cover in continuity with the trend until nowaday


# export the output
# you made the calculation and you want to send the output to a collegue, you can use writeRaster func.

writeRaster(prediction, "final.tif")  # is good func for transform images in tif, It create new data 1. we import data (5images), 2. we stack all that ims and added prediction, 3. with writeRaster
# writeRaster is the oppost of Raster function!! that import data into R. at the opp. writeRaster you can export the image from R to your PC folder.
# final stack (to see if the extent of the image is different)
# https://gdal.org/ to have a look in which format will work my function


final.stack <- stack(snow.multitemp, prediction)   # we are going to stack the all multitemp+ prediction of snow cover
plot(final.stack, col=cl)
# export the R graph for a thesis in pdf format

pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()

png("my_final_exciting_graph.png")
plot(final.stack, col=cl)
dev.off()
