# Species Distrubution Modelling

install.packages("sdm")     # install our packages very useful for ecology and where we can find our dataset of species distribution.
# install.packages("rgdal")
library(sdm)
library(raster)            # library to read the raster and to make use the raster dataset # raster packages are dependent from sp packages, you can solve the problem with " library( rgdal, "dependencies=T"
library(rgdal)             # geodata astract library # use to import the packs of the species data spatial distrubution for this project

# species :  we are looking for a certain species
file <- system.file("external/species.shp", package="sdm")   # automatically importation of the external folder, insiede the folder external there is the species
                                                             # system.file function need to import the file into the sdm package
species <- shapefile(file)                                   # we can use the graphical part of the file using shapefile func., this type of files are very used today for graphical and mapping dataset

species                                                      # looking at the estent and features that are the points in this case
species$Occurrence                                           # let's see the occurence of species, present-absent model (200 features)
plot(species)

plot(species[species$Occurrence == 1,],col='blue',pch=16)   # making a condition inside the "[]" and to make it we need "==" command, we decided only the presence occurences = 1

points(species[species$Occurrence == 0,],col='red',pch=16)  # making the oppost condition= absence with red points (blue points are blue)

# environmental variables
path <- system.file("external", package="sdm")             # making the path inside the external using the same function of before system.file. the folder is called external 

lst <- list.files(path=path,pattern='asc$',full.names = T) # the extention of the file in this case is "asc" (other extention may be tif or png, ect), dollar can be avoided, but to be sure to carry out all data
lst   # have a look of files into the list, into smd pkgs there is external folder, into external there are species and ecological and environmental variables

preds <- stack(lst)                   # we are making a stack of dataset, elevation, temperature, precipitation, vegetation

cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)                  # we can analyze the predictor informations  graphically, we fond that many occorences are at the bottom of the valley, in case of vegetation presence because Brachipodium rupestre love small elevation and grasslands

plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16)  # Brachipodium rupestre love high temperature and low elevation

plot(preds$precipitation, col=cl)                  # Brachipodium rupestre love medium to high precipitation, in the case medium precipitation
points(species[species$Occurrence == 1,], pch=16)

plot(preds$vegetation, col=cl)                     # vegetation index is giving from the red-nir index, Brachipodium rupestre love high ammount vegetation or not. this analysis is called exploration analysis. Brachipodium rupestre love medium vegetation
points(species[species$Occurrence == 1,], pch=16)

                                                   # making the model thanks all this information. one is train = species the other argument of the model is predictors= preds.

d <- sdmData(train=species, predictors=preds)      # class = sdmData, n speci= 1, number or record 200, type of data= presence-absent, plot, abbundance of species but we are looking only prec. absen. is good for the velocity of the sampling efforce.
d                                                  # there is a negative correlation of linear model between species occ. and elevation. In general we use a logistic model that reppresent better ecological functions

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods='glm') # occurence, ~ is equal symble, y~a+bx in this case x is negative related with y in case of elevation. model: y= a+ b elev. + c + temperature + d precipitation + e vegetation, for all the variable we are using a different curve.
                                                   # finally the methods of our model (glm) : generalised linear model is the model that we are going to use.
p1 <- predict(m1, newdata=preds)                   # final prediction, related at our predictors used

plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)  

s1 <- stack(preds, p1)                            # stacking all the predictors together with the final layer added
plot(s1, col=cl)
