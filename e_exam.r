## R_code_MyProject_exam.r elisa

# R code for exam project

#library(rgdal)
#library(gdalUtils)
library(raster)
library(RStoolbox)

#set the working directory
setwd("C:/lab/exam/")

# Sentinel-2 L1C. Date 16-10-2019. Localisation: Catalonia, datum=WGS84, coordinates in UTM
# optain GeoTIFF bands using rgdal and gdalUtils packages

#gdal_translate("T30SWG_20191016T110041_B02.jp2", "B02.tif")
#gdal_translate("T30SWG_20191016T110041_B03.jp2", "B03.tif")
#gdal_translate("T30SWG_20191016T110041_B04.jp2", "B04.tif")
#gdal_translate("T30SWG_20191016T110041_B08.jp2", "B08.tif")
#gdal_translate("T30SWG_20191016T110041_B011.jp2", "B011.tif")

#octB02 <- raster("B02.tif")
#octB03 <- raster("B03.tif")
#octB04 <- raster("B04.tif")
#octB08 <- raster("B08.tif")
#octB11 <- raster("B11.tif")

## or elaborate them directly in SNAP 

## next ##
# Creating a RasterStack

# first of all, resampling B11:
b11 <- raster("subset_0_of_T30SWG_20191016T110041_B11.tif") # move in C:/lab/ to not be in conflict with the new file
b11_dis <- disaggregate(b11, fact=2) # B SWIR2 from resolution 20m to 10m as the other bands
writeRaster(b11_dis, "sub_11.tif") # into C:/lab/exam/

# rasterStack of bands 2,3,4,8,11
rlist <- list.files(pattern="sub") # or pattern="B" using gdalUtils
import <- lapply(rlist, raster)
dry_land_2019 <- stack(import)
plotRGB(dry_land_2019, 4,3,2, stretch="lin")

# vegetation indicies: NDVI and NDWI2
# require RStoolbox

ndvi <- spectralIndices(dry_land_2019 , red = "subset_2_of_T30SWG_20191016T110041_B04", nir = "subset_3_of_T30SWG_20191016T110041_B08", indices = "NDVI")
ndwi2 <- spectralIndices(dry_land_2019 , nir = "subset_3_of_T30SWG_20191016T110041_B08", swir2 = "sub_11", indices = "NDWI2")
ndvi_bit <- stretch(ndvi, minv=0, maxv=255) # to compare ndvi and ndwi2 with the same range of pixel values
ndwi2_bit <- stretch(ndwi2, minv=0, maxv=255)

cl <-  colorRampPalette(c("black", "green", "red"))(100)
clb <-  colorRampPalette(c("black", "gold", "blue"))(100)

# displaying on one map the indicies
par(mfrow=c(2,1))
plot(ndvi_bit, col=cl, main="NDVI")
plot(ndwi2_bit, col=clb, main="NDWI2")

# albedo: Black Sky Albedo, satellite MCD43A3, 500m resolution, bands: 1(red),2(NIR),3(blue),4(green), Date: 16-10-2019, coordinates in decimal degrees lon -2.941667, -2.083333, lat 36.95417, 37.85417 
b1 <- raster("MCD43A3.006_Albedo_BSA_Band1_doy2019289_aid0001.tif")
b2 <- raster("MCD43A3.006_Albedo_BSA_Band2_doy2019289_aid0001.tif")
b3 <- raster("MCD43A3.006_Albedo_BSA_Band3_doy2019289_aid0001.tif")
b4 <- raster("MCD43A3.006_Albedo_BSA_Band4_doy2019289_aid0001.tif")

albedo <- stack(b1, b2, b3, b4)
plot(albedo)

# reducing the dimentions for having a rappresentation of all of the bands together
pairs(albedo)# have a fast look at the correlation among bands
albedoPCA2019 <- rasterPCA(albedo)
summary(albedoPCA$model)# PC1 describes 94%
plot(albedoPCA$map)

# checking and displaying on a map NDVI and Albedo
par(mfrow=c(1,2))
plot(albedoPCA2019$map$PC1)
plot(ndvi_bit)

# Applying Linear Regression for albedo and NDVI

# align coordinates systems, extent and ncell
pr <- projectRaster(ndvi, albedoPCA2019$map$PC1)
rm <- stack(pr, albedoPCA2019$map$PC1)
names(rm) <- c("NDVI", "AlbedoPCA")
cla <-  colorRampPalette(c("white", "red", "blue"))(100)
plot(rm, col=cla) #displaying the result and comparing ndvi and albedo

ndvi2019_ex <- extract(rm$NDVI, c(1:1000, 1000, 1000))
albedo2019_ex <- extract(rm$AlbedoPCA, c(1:1000, 1000, 1000))
model <- lm(ndvi2019_ex ~ albedo2019_ex) # y is ndvi, dipendent variable; x is albedo, indipendent variable
plot(albedo2019_ex, ndvi2019_ex)
abline(model, col="red")
summary(model) #R-squared:  0.4976, p-value: < 2.2e-16

# Applying Linear Regression for albedo and NDWI2: as a proxy of soil moisture content
prw <- projectRaster(ndwi2, albedoPCA2019$map$PC1)
rmw <- stack(prw, albedoPCA2019$map$PC1)
names(rm) <- c("NDWI2", "AlbedoPCA")
plot(rmw)

ndwi22019_ex <- extract(rmw$NDWI2, c(1:1000, 1000, 1000))
albedo2019_ex <- extract(rmw$AlbedoPCA, c(1:1000, 1000, 1000))
modelw <- lm( ndwi22019_ex ~ albedo2019_ex)
plot(albedo2019_ex, ndwi22019_ex)
abline(modelw, col="green")
summary(modelw) # R-squared:  0.3429,  p-value: < 2.2e-16

par(mfrow=c(1,2))
plot(albedo2019_ex, ndvi2019_ex, xlab="Albedo 2019", ylab="NDVI 2019", main="NDVI ~ Albedo" )
abline(model, col="red")
plot(albedo2019_ex, ndwi22019_ex, xlab="Albedo 2019", ylab="NDWI2 2019", main="NDWI2 ~ Albedo")
abline(modelw, col="green")

##The Terra Moderate Resolution Imaging Spectroradiometer (MODIS) Vegetation Indices (MOD13Q1) Version 6, 16 days at 250 meter (m) spatial resolution as a Level 3 product
# SDS name: 250m 16 days NDVI
ndvi_2019 <- raster("MOD13Q1_NDVI_2019.tif")
ndvi_2015 <- raster("MOD13Q1_NDVI_2015.tif")
ndvi_2010 <- raster("MOD13Q1_NDVI_2010.tif")

#ndvi layers need to have the same extent of albedo2019 to compare them with albedo layers
extent <- c(-2.941667, -2.083333, 36.95417, 37.85417)
ndvi_2019ext <- crop(ndvi_2019, extent)
ndvi_2015ext <- crop(ndvi_2015, extent)
ndvi_2010ext <- crop(ndvi_2010, extent)

#Difference in biomass 
dif_ndvi <- ndvi_2019ext- ndvi_2010ext
clc <- colorRampPalette(c('red','green','yellow')) (100)
plot(dif_ndvi, col= clc)
hist(dif_ndvi)

# rise or decrease in biomass
plot(ndvi_2010ext, ndvi_2019ext, xlab="NDVI 2010", ylab="NDVI 2019", main="Trend of NDVI")
abline(0,1,col="red") # 45° line that divides the cartesian plane in 2 equal parts

NDVI <- stack(ndvi_2010ext, ndvi_2015ext, ndvi_2019ext)
boxplot(NDVI,outline=F, horizontal=T, axes=T, names=c("NDVI 2010", "NDVI 2015", "NDVI 2019"), main="Boxplot NDVI") 

## The Moderate Resolution Imaging Spectroradiometer (MODIS) MCD43A3 Version 6 Albedo Model dataset, 16 days of Terra and Aqua MODIS data at 500 meter (m) resolution
# Date 16 October 2015 and 2010

# 2015
alb_2015b1 <- raster("MCD43A3.006_Albedo_BSA_Band1_doy2015.tif")
alb_2015b2 <- raster("MCD43A3.006_Albedo_BSA_Band2_doy2015.tif")
alb_2015b3 <- raster("MCD43A3.006_Albedo_BSA_Band3_doy2015.tif")
alb_2015b4 <- raster("MCD43A3.006_Albedo_BSA_Band4_doy2015.tif")

# albedo layers need to have the same extent of albedo2019
extent <- c(-2.941667, -2.083333, 36.95417, 37.85417)
alb_2015b1ex <- crop(alb_2015b1, extent)
alb_2015b2ex <- crop(alb_2015b2, extent)
alb_2015b3ex <- crop(alb_2015b3, extent)
alb_2015b4ex <- crop(alb_2015b4, extent)
albedo_2015 <- stack(alb_2015b1ex, alb_2015b2ex, alb_2015b3ex, alb_2015b4ex)
albedo_2015PCA <- rasterPCA(albedo_2015)
summary(albedoPCA2015$model)

# 2010
alb_2010b1 <- raster("MCD43A3.006_Albedo_BSA_Band1_doy2010.tif")
alb_2010b2 <- raster("MCD43A3.006_Albedo_BSA_Band2_doy2010.tif")
alb_2010b3 <- raster("MCD43A3.006_Albedo_BSA_Band3_doy2010.tif")
alb_2010b4 <- raster("MCD43A3.006_Albedo_BSA_Band4_doy2010.tif")

extent <- c(-2.941667, -2.083333, 36.95417, 37.85417)
alb_2010b1ex <- crop(alb_2010b1, extent)
alb_2010b2ex <- crop(alb_2010b2, extent)
alb_2010b3ex <- crop(alb_2010b3, extent)
alb_2010b4ex <- crop(alb_2010b4, extent)
albedo_2010 <- stack(alb_2010b1ex, alb_2010b2ex, alb_2010b3ex, alb_2010b4ex)
albedo_2010PCA <- rasterPCA(albedo_2010)
summary(albedoPCA2010$model)

# in my opinion faster passage
writeRaster(albedoPCA2019$map$PC1, "albedoPCA2019.tif")
writeRaster(albedo_2010PCA$map$PC1, "albedoPCA2010.tif")
writeRaster(albedo_2015PCA$map$PC1, "albedoPCA2015.tif")

albedoPCA2010 <- raster("albedoPCA2010.tif")
albedoPCA2015 <- raster("albedoPCA2015.tif")
albedoPCA2019 <- raster("albedoPCA2019.tif")

# Difference in albedo 
albedoPCA2019b <- stretch(albedoPCA2019, minv=0, maxv=255)
albedoPCA2015b <- stretch(albedoPCA2015, minv=0, maxv=255)
albedoPCA2010b <- stretch(albedoPCA2010, minv=0, maxv=255)

dif_albedo <- albedoPCA2019b - albedoPCA2010b
plot(dif_albedo)
plot(albedoPCA2010b, albedoPCA2019b, xlab="Albedo 2010", ylab="Albedo 2019", main="Trend of Albedo")
abline(0,1,col="red")

albedobox <- stack(albedoPCA2010b, albedoPCA2015b, albedoPCA2019b)
boxplot(albedobox, outline=F, horizontal=T, axes=T, names=c("Albedo 2010", "Albedo 2015", "Albedo 2019"), main="Boxplot Albedo")



# let's check if with higher amount of data multitemporal analysis of NDVI changes
# satellite data day by day for the month of October

library(MODISTools)
VI <- mt_subset(product = "MOD13Q1", band = "250m_16_days_NDVI", lat = 37, lon = -2, start = "2019-10-01", end = "2019-10-31", km_lr = 100, km_ab = 100, site_name = "testsite", internal = TRUE, progress = FALSE)

V2 <- mt_subset(product = "MOD13Q1",
                band = "250m_16_days_NDVI",
                lat = 37,
                lon = -2,
                start = "2015-10-01",
                end = "2015-10-31",
                km_lr = 100,
                km_ab = 100,
                site_name = "testsite",
                internal = TRUE,
                progress = FALSE)

V3 <- mt_subset(product = "MOD13Q1",
                band = "250m_16_days_NDVI",
                lat = 37,
                lon = -2,
                start = "2010-10-01",
                end = "2010-10-31",
                km_lr = 100,
                km_ab = 100,
                site_name = "testsite",
                internal = TRUE,
                progress = FALSE)

VI_r <- mt_to_raster(df = VI)
V2_r <- mt_to_raster(df = V2)
V3_r <- mt_to_raster(df = V3)

multit_NDVI <- stack(VI_r, V2_r, V3_r)
# multitemporal NDVI needs to have the same coordinates system and extent of previous NDVI data to compare the boxplot results
multit_NDVIex <- projectRaster(multit_NDVI,alb_2015b1ex)
boxplot(multit_NDVIex,outline=F, horizontal=T, axes=T, names=c("NDVI 2010", "NDVI 2015", "NDVI 2019"), main="Boxplot multitemporal NDVI")
