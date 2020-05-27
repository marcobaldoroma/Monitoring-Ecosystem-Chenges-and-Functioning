### Enviroment biodiversity variation measurement
window <- matrix(1, nrow = 3, ncol = 3)
window
std_sntr <- focal(snt_r$prova5_.4, w=window, fun=sd)

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) # 
plot(std_sntr, col=cl)

par(mfrow=c(1,2))
plot(std_sntr, col=cl)
# plot(std_snt8bit, col=cl)

std_sntr1 <- focal(snt_r$prova5_.1, w=window, fun=sd)

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) # 
plot(std_sntr1, col=cl)

### PCA related sd
library(RStoolbox)
sntrpca <- rasterPCA(snt_r)

summary(sntrpca$model) 

clp <- colorRampPalette(c('dark grey','grey','light gray'))(100) # 
plot(sntrpca$map,col=clp)

plotRGB(sntrpca$map,1,2,3,stretch="lin")

std_sntrpca <- focal(sntrpca$map$PC1, w=window, fun=sd)

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) # 
plot(std_sntrpca, col=cl)

##############



# R code Ecosystem Biodiversity Variable

library (raster)

setwd("C:/lab/")

# ratser function can import a single layer band but could be useful
# brick can import all the layer bands
snt <- brick("sentinel.tif")

snt

# 2^8 # bits
[1] 256
> 2^12
[1] 4096
> 2^13
[1] 8192

# Boa Vista Brazil

# plot this images
plot(snt)

# B1 blue
# B2 green
#B3 red
#B4 NIR

# R G B
plotRGB(snt,3,2,1, stretch="lin")
plotRGB(snt,4,3,2, stretch="lin") #putting the nir on the top layer

# standard deviation is possible to calculate only in one layer at time.
# we can do it with the multivariate analysis, in the way to have only 1 dimention instead than one PCA
#have a look of the relationship of the bands
pairs(snt)

# PCA analysis

# install.packages("RStoolbox")
library(raster)
library(RStoolbox) # this is for PCA

sntpca <- rasterPCA(snt)  # to select our rater PCA
sntpca

# we want link the model at our images

summary(sntpca$model) 
#Cumulative Proportion    0.7015076  means that we are using 70% of the original information

plot(sntpca$map)

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")  # in this way we can have different coloration of the area based on the bands 

# set the moving window, we say that there are 5x5 pixels at time, mesuring the standard deviation of window

window <- matrix(1, nrow= 5, ncol = 5)



# function focal is using, see the description, mathematic analysis that we can do with focal func.

sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)


cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) # 
plot(sd_snt, col=cl)

# we want plot the two images

par(mfrow= c(1,2))
plotRGB(snt,4,3,2, stretch="lin", main="original image") 
plot(sd_snt, col=cl, main="diversity")


#############  Day 2  R_code_EBVs.r second part on cladonia_stellaris_calaita

setwd("C:/lab/")

library (raster) # we need for use brick and raster function

library(RStoolbox) # we use for the PCA of the image

# we use brick function to import all the layers of our image
clad <- brick ("cladonia_stellaris_calaita.JPG")

# first blue, second green, third red
plotRGB(clad, 1,2,3, stretch = "lin")

# we have to decide our moving window, to calculate the SD and report the SD on our pixel. 1 is to set the value of the moving window pixel
window <- matrix(1, nrow= 3, ncol= 3)
window

# focal function have the capacity to make the calcolation of SD
cladpca <- rasterPCA(clad)
cladpca  # we can see the output(informations) for cladpca

summary(cladpca$model)
# 98% of the images informations is showing from my model

# inside the clad pca the is the map of pca component, so link the pca with map with PC1 # calcolate focal
sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd)

# we want aggregate the image for make the sd analysis faster and make the focal of our aggregation

PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)

# make the powerful col for our function and use par to plot a frame 
# plot the calculation both sd and sd agg

par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl)

# we want see the original image and secon the variability of the image. to see the variability of each individual, so the inner variation with sd calculation

par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) #
plotRGB(clad, 1,2,3, stretch = "lin")
plot(sd_clad, col=cl)





