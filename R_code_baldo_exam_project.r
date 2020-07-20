# Marco Baldo Exam Project

#install.packages("GGally") 
#install.packages("raster")    
#install.packages("RStoolbox")  
#install.packages("rasterdiv")  
#install.packages("rasterVis")
#install.packages("ncdf4")

library(GGally)
library(raster)
library(RStoolbox)
library(rasterdiv)
library(rasterVis)
library(ncdf4) 

# The product is displayed in a regular latitude/longitude grid with the ellipsoïd WGS 1984 (Terrestrial radius=6378km). The resolution of the grid is 1°/336. The reference is the centre of the pixel. Subtract (for longitude) or add (for latitude) half the angular resolution to get the pixel top left coordinates.
# Images have been downloaded from VITO Earth Observation website. Mission Proba-V ESA, sensor Vegetation, resolution= moderate, grid size=1/3KM, datum=WGS84, coordinates in UTM, type of files NC,tiff.
# First of all, I'm interested to have a world view of Dry Matter Productivity (Kg/ha/day).
# Definition: Dry matter Productivity (DMP) represents the overall growth rate or dry biomass increase of the vegetation and is directly related to ecosystem Net Primary Productivity (NPP), however with units customized for agro-statistical purposes (kg/ha/day). Similarly the Gross Dry Matter Productivity (GDMP) is equivalent to Gross Primary Productivity (GPP).The main difference between DMP and GDMP lies in the inclusion of the autotrophic respiration.

setwd("D:/exam/DMP/")                                                                     # set wd for DMP300 files

dmp2014 <- raster("c_gls_DMP300-RT5_QL_201406100000_GLOBE_PROBAV_V1.0.1.TIFF")            # import all the temporal images with raster function
dmp2015 <- raster("c_gls_DMP300-RT5_QL_201506100000_GLOBE_PROBAV_V1.0.1.TIFF")            # the day is always Jun 10th, I didn't use at the begging the Stack function because I needed the files separated for the presentation
dmp2016 <- raster("c_gls_DMP300-RT5_QL_201606100000_GLOBE_PROBAV_V1.0.1.TIFF")  
dmp2017 <- raster("c_gls_DMP300-RT5_QL_201706100000_GLOBE_PROBAV_V1.0.1.TIFF")  
dmp2018 <- raster("c_gls_DMP300-RT5_QL_201806100000_GLOBE_PROBAV_V1.0.1.TIFF")  
dmp2019 <- raster("c_gls_DMP300-RT5_QL_201906100000_GLOBE_PROBAV_V1.0.1.TIFF")  

par(mfrow=c(2,3)) 
cl <- colorRampPalette(c("yellow","light green","dark green", "black"))(100)              # we want plot all images together for this we used a multiframe par function, even if I know there are several other ways to plot all my img set( i.g lapply and stack)
plot(dmp2014, col=cl, main="2014")                                                         
plot(dmp2015, col=cl, main="2015")
plot(dmp2016, col=cl, main="2016")
plot(dmp2017, col=cl, main="2017")
plot(dmp2018, col=cl, main="2018")
plot(dmp2019, col=cl, main="2019")
dev.off()

par(mfrow=c(2,1))
plot(dmp2016, col=cl, main="2016") # the driest year
plot(dmp2018, col=cl, main="2018") # the wettest year
dev.off()

# very impressive DMP seasonality changing of the northen hemisphere

##########################DMP300
# Raster Stack (object with the same spatial extent) of DMP images dataset. Resolution 1/3km grid. File format nc period winter-summer.
setwd("D:/exam/DMP300")

rlistdmp300 <- list.files(pattern="DMP300") 
rlistdmp
import <- lapply(rlistdmp300, raster)
dmp300.multitemp <- stack(import)
names(dmp300.multitemp) <- c("Jan 2014","Jul 2014","Jan 2016","Jan 2018","Jul 2018","Jan 2020")
plot(dmp300.multitemp)

# to crop images on my AOI (UPSP; 12°32'15"N lat. and 75°39'46"E long.)
ext <- c(75.10, 75.60, 12.05, 12.55)                                            
dmp300.up <- crop(dmp.multitemp, ext)                                              # up= Uppangala, the Indian study area
plot(dmp300.up)

# Having a comparison of the DMP in north hemisphere with levelplot func. 2019 May 31st/ 10th 
# Loading of my image for the levelplot func. DMP comparison
dmp201905 <- raster("c_gls_DMP300-RT5_QL_201905310000_GLOBE_PROBAV_V1.0.1.TIFF")

# Levelplot function of may-jun DMP2019 images. Having a view of the DMP north hemisphere changing.
levelplot(dmp201905, main="May 31st 2019")
levelplot(dmp2019, main="June 10th 2019")

##########################LAI300
# Vegetation property analyzed is the Leaf Area Index but to have a look on the tree's potential in condensing water vapour
# Let's set my working directory for the project
setwd("D:/exam/LAI300")                            

# Stack and and crop of Global Leaf Area Index images for my analysis. file=nc format, I will import only one raster layer with lapply function to be faster.

rlistlai <- list.files(pattern="LAI300")  
rlistlai
import <- lapply(rlistlai, raster) 
lai.multitemp <- stack(import)
# names(lai.multitemp) <- c("14","15","16","17","18","19") all the images are Jun 20th of multitemp scale 2014-2019
# plot(lai.multitemp)
ext2 <- c(73, 78, 10, 15) 
lai.wg <- crop(lai.multitemp, ext2)                                          
plot(lai.wg)

ext1 <- c(75, 76, 12, 13) 
biotic.pump <- crop(lai.multitemp, ext1)                                          
plot(biotic.pump)

# to crop images on my AOI (UPSP; 12°32'15"N lat. and 75°39'46"E long.)
ext <- c(75.10, 75.60, 12.05, 12.55)                                           # I will use this ext also for my DMP analysis                                                    
lai.up <- crop(lai.multitemp, ext)                                             # up= Uppangala, my Indian study area                              
plot(lai.up)

##########################DMP
## the same procedure of LAI dataset, but on the Dry Matter Productivity and with aggregated pixels images for a very fast and even visive pixels statistical analysis
setwd("D:/exam/DMP")

rlistdmp <- list.files(pattern="DMP") 
rlistdmp
import <- lapply(rlistdmp, raster)
dmp.multitemp <- stack(import)
names(dmp.multitemp) <- c("14","15","16","17","18","19.5","19")
#plot(dmp.multitemp) # already plotted at the begining
# to crop images on my AOI (UPSP; 12°32'15"N lat. and 75°39'46"E long.)
#ext <- c(75.10, 75.60, 12.05, 12.55)                                            
dmp.up <- crop(dmp.multitemp, ext)                                              # up= Uppangala, the Indian study area
plot(dmp.up)

# correlation analysis of my images with ggpairs func. # require (GGally) 
ggpairs(dmp.up)

# raster PCA to analyze all the DMP.UP images in one single component
dmpPCA <- rasterPCA(dmp.up)
summary(dmpPCA$model)     # PC1 model describes 93.7%
plot(dmpPCA$map)      

# High Dry Matter Production areas(pixels) enounce that the DMP is higher in general where It was higher also in the previesly years (where the forest cover is higher and well preserved)
dmp2014.up <- crop(dmp2014, ext)
dmp2015.up <- crop(dmp2015, ext)
dmp2016.up <- crop(dmp2016, ext)
dmp2017.up <- crop(dmp2017, ext)
dmp2018.up <- crop(dmp2018, ext)
dmp2019.up <- crop(dmp2019, ext)

# In Uppangala case study seems that the older part of the forest enounce a better resilience of the DMP Ecosystem Function strectly related at the Biomass Production and Carbon Storage Ecosystem Service
plot(dmp2014.up, main="UP.2014")
plot(dmp2015.up, main="UP.2015")
plot(dmp2016.up, main="UP.2016")
plot(dmp2017.up, main="UP.2017")
plot(dmp2018.up, main="UP.2018")
plot(dmpPCA$map$PC1, main="PC1")  # best predictor

# Differences in dry matter productivity in UP.
plot(dmp2014.up, dmp2019.up, xlab="DMP.UP 2014", ylab="DMP.UP 2019", main="Trend of DMP.UP")
abline(0,1,col="red") # to understand better the cloud point distribution we use ab 45° line that divides in 2 equal parts the cartesian plane

# Another check on DMP trend in UP. Differenzial among pixels. Histogram of the differential.
dif_dmp.up <- dmp2019.up - dmp2014.up
cl <- colorRampPalette(c("orange","white","light green","dark green")) (100)
par(mfrow=c(1,2))
plot(dif_dmp.up, col= cl, main="DMP.UP Differences")
hist(dif_dmp.up, main="histogram")

# to detect the DMP trend year by year using boxplot function
DMP <- stack(dmp2014.up, dmp2015.up, dmp2016.up, dmp2017.up, dmp2018.up, dmp2019.up)
boxplot(DMP,outline=F, horizontal=T, axes=T, names=c("DMP.UP 2014", "DMP.UP 2015", "DMP.UP 2016","DMP.UP 2017","DMP.UP 2018","DMP.UP 2019"), main="Boxplot DMP.UP")       # cancel the outliners 



##################################DMP300
# we want have a comparison with the trend result of previusly analysis (DMP in TIFF format)
# keeping back the stack of the code's top (DMP300)
# Importation of raster layer only of Jan. images

dmp3002014 <- raster("c_gls_DMP300-RT5_201401100000_GLOBE_PROBAV_V1.0.1.NC")            
dmp3002016 <- raster("c_gls_DMP300-RT5_201601100000_GLOBE_PROBAV_V1.0.1.NC")  
dmp3002018 <- raster("c_gls_DMP300-RT5_201801100000_GLOBE_PROBAV_V1.0.1.NC")  
dmp3002020 <- raster("c_gls_DMP300-RT5_202001100000_GLOBE_PROBAV_V1.0.1.NC")  

# cropping on Uppangala
dmp3002014.up <- crop(dmp3002014, ext)
dmp3002016.up <- crop(dmp3002016, ext)
dmp3002018.up <- crop(dmp3002018, ext)
dmp3002020.up <- crop(dmp3002020, ext)

# In Uppangala case study seems that the older part of the forest enounce a better resilience of the DMP Ecosystem Function strectly related at the Biomass Production and Carbon Storage Ecosystem Service
#par(mfrow=c("2,2"))
#plot(dmp3002014.up, main="UP.2014")
#plot(dmp3002016.up, main="UP.2016")
#plot(dmp3002018.up, main="UP.2018")
#plot(dmp3002020.up, main="UP.2020")


# Differences in dry matter productivity in UP.
plot(dmp3002014.up, dmp3002020.up, xlab="DMP.UP 2014", ylab="DMP.UP 2020", main="Trend of DMP.UP")
abline(0,1,col="red") # to understand better the cloud point distribution we use ab 45° line that divides in 2 equal parts the cartesian plane

# Another check on DMP trend in UP.
dif_dmp.up2 <- dmp3002020.up - dmp3002014.up
cl <- colorRampPalette(c("orange","white","light green","dark green")) (100)
par(mfrow=c(1,2))
plot(dif_dmp.up2, col= cl, main="DMP.UP Differences")
hist(dif_dmp.up2, main="histogram")

# to detect the DMP trend year by year
DMP <- stack(dmp3002014.up, dmp3002016.up, dmp3002018.up, dmp3002020.up)
boxplot(DMP,outline=F, horizontal=T, axes=T, names=c("DMP.UP 2014", "DMP.UP 2016","DMP.UP 2018","DMP.UP 2020"), main="Boxplot DMP.UP")
######################################################################################################################################################################
