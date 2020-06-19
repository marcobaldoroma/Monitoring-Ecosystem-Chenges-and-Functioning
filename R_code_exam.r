# R_code_exam.r

#############################################################################################################################
#############################################################################################################################

# 1. R_code_first.r
# 2. R_code_multipanel.r
# 3. R code spatial.r 
# 4. R_code_point_pattern_analysis.r 
# 5. R_code_multivar.r
# 6. R_code_remote_sensing.r
# 7. R_code_ecosystem_function.r
# 8. R_code_PCA_remote_sensing.r   vedi multivariate analysis remote sensing su iol
# 8. R_code_reflectance.r
# 9. R_code_faPAR.r
# 10. R_code_EBVs.r
# 11. R_code_snow.r
# 12. R_code_NO2.r
# 13. R_code_crop_image.r
# 14. R_code_interpolation.r
# 15. R_code_sdm.r

#############################################################################################################################
#############################################################################################################################

# 1. R_code_first.r
# first test to start to use functions and codes necessary to work in R environment
# sp is a package providing classes and methods for spatial data: points, lines, polygons and grids
# the quotes marks are used to call objects and packages from outside R environment or to delimite them
# Single and double quotes delimit character constants. They can be used interchangeably but double quotes are preferred (and character constants are printed using double quotes), so single quotes are normally only used to delimit character constants containing double quotes
install.packages("sp") 

# call the required and installed packages in our work environment
library(sp)

# data function: Loads specified data sets, or list the available data sets. Load meuse that is a data set inside the sp packages
# meuse dataset: this data set gives locations and topsoil heavy metal concentrations, along with a number of soil and landscape variables at the observation locations, collected in a flood plain of the river Meuse, near the village of Stein (NL). Heavy metal concentrations are from composite samples of an area of approximately 15 m x 15 m
data(meuse)                

# Let's see how meuse dataset is structure:
meuse

# let's look at the first rows of the data set, head of the data set
head(meuse)

# attach function: the database is attached to the R search path. This means that the database is searched by R when evaluating a variable, so objects in the database can be accessed by simply giving their names
# attach meuse dataset
attach(meuse)
# let's see if the zin concentration is realate to the copper concentration
# let's plot the variables together, zinc in the x axis and copper in the y axis
plot(zinc,copper)
# col= "green" to give a colour at the objects and symbols of the plot
plot(zinc,copper,col="green")
# pch function is used to change symbols associated at the variables rappresentated (i.g. pch=19 octagon)
plot(zinc,copper,col="green",pch=19)
# cex= character exaggeration: a function number indicating the amount by which plotting text and symbols should be scaled relative to the default. 1=default, 1.5 is 50% larger, 0.5 is 50% smaller, etc. to change the size of symbols, you can put 0<x<1 to have a smaller symbol size or x>1 to have a bigger symbol size of the rappresented variables position.
plot(zinc,copper,col="green",pch=19,cex=2)

############################################################################################################################
############################################################################################################################

# 2. R_code_multipanel.r

### Multipanel in R: seeing correlation amoung ecological variables

                                                 # install.packages("sp")
install.packages("GGally")                       # the 'GGally' extends 'ggplot2'(R package 'ggplot2' is a plotting system based on the grammar of graphics) by adding several functions to reduce the complexity of combining geometric objects with transformed data. Some of these functions include a pairwise plot matrix, a two group pairwise plot matrix, a parallel coordinates plot, a survival plot, and several functions to plot networks. It is required to use the function ggpairs() we would graphy

                                                 # calling required libraries with library command/function
library(sp)
library(GGally)

data(meuse)                                      # load a data set. There is a dataset available named meuse in sp

                                                 # attach sp function: The database is attached to the R search path. This means that the database is searched by R when evaluating a variable, so objects in the database can be accessed by simply giving their names. Attach function to fix the data set in R environment for the plot analysis
attach(meuse)

                                                 # Exercise: see the names of the variables. Second plot cadmium versus zinc
                                                 # There are two ways to look the names of the variables:
names(meuse)                                     # names functions to get or set the names of an object. Names is a generic accessor function, and names<- is a generic replacement function. The default methods get and set the "names" attribute of a vector (including a list) or pairlist
head(meuse)

meuse                                            # to see all the variables
                                                 # plot cadmium and zinc with symbol 15 of R plot symbol set, character exaggeration to make the size of the symbles more visible
plot(cadmium,zinc,pch=15,col="red",cex=2)

                                                 # Exercise: make all the possible pairs plots of the dataset. to be faster I wrote a generic code where x = copper,lead, zinc, cadmium that changing averytimes depending on the other viriable present
# plot(x,cadmium)
# plot(x,zinc)
# plot(x,lead)
# plot(x,copper)
                                                 # sometimes plot is not a good idea. For example in the pairs analysis
pairs(meuse)                                     # The pairs R function returns a plot matrix, consisting of scatterplots for each variable-combination of a data frame. Pairs is a function stored in GGally and ggplot2. 
                                                 # selection of variables that we want to plot all together from a data set (meuse in our case)
pairs(~ cadmium + copper + lead + zinc, data=meuse)
                                                 # the alternative way could be select a number of rows and collumns from the data set 
pairs(meuse[,3:6])                               # all rows and collumn from 3(cadmium) to 6(zinc)

                                                 # Exercise: prettify this graph
pairs(meuse[,3:6], col="red", pch=18, cex=1.5)

pairs(meuse[,3:6], pch=19)

                                                 # GGally packages necessary to use ggpairs function. It will show infographic even better of pairs function
ggpairs(meuse[,3:6])

############################################################################################################################
############################################################################################################################

# 3. R code spatial

# first real spatial analysis and view of spatial data set gave by us
# install this two packages for spatial analysis "raster" "rgdal"
## R code for spatial analysis and views of data set

install.packages("raster")                                # raster package: reading, writing, manipulating, analyzing and modeling of gridded spa-tial data. The package implements basic and high-level functions. Process-ing of very large files is supported. There is a also support for vector data operations such as in-tersections
install.packages("rgdal")                                 # rgdal package: provides bindings to the 'Geospatial' Data Abstraction Library ('GDAL') and access to projection/transformation operations from the 'PROJ' library. Use is made of classes defined in the 'sp' package. Raster and vector map data can be imported into R, and raster and vector 'sp' objects exported
                                                          # to call the required libraries (sp, raster, rgdal)
library(sp)
library(raster)
library(rgdal)                                            # to load the data set meuse

data(meuse)
                                                          # to have a look at the top lines of the data set
head(meuse)
                                                          # coordinates is the function to sets, converts the objects into sp classes, and sets coordinate reference system to the variables

coordinates(meuse) = ~x+y                                 # alt 126 for the tilde in windows / "sp" function is a R old package. Tilde Operator is used to separate the left- and right-hand sides in a model formula
                                                          # to plot meuse objects
plot(meuse)
                                                          # sp package, spplot function: a wrapper function around spplot (sp package). With spplot it is easy to map several layers with a single legend for all maps. ssplot is itself a wrapper around the levelplot function in the lattice package, and see the help for these functions for additional options
spplot(meuse, "zinc")                                     # simple plot graphs point and vectors, using spplot you can use a variables as functions and graphy the function instead of the variables i.g. concentration 
                                                          # spplot give different colours at the points, ordering the values of the variable in different classes of values
                                                          # Exercise: plot the spatial amount of copper using simple plot
spplot(meuse, "copper")

                                                          # main function: add a label at the graphic
spplot(meuse, "copper", main= "Copper concentration ")
                                                          # just to have a points size directly proportionated at the concentration (variable or object of our analysis) we can use the sp bubble function
                                                          # bubble function: create a bubble plot of spatial data, with options for bicolour residual plots (xyplot wrapper)
bubble(meuse,"zinc", main= "Zinc concentration")          # the bubble function is going to make the different size of points positive relate at the variable

                                                          # Exercise: use bubble funct. for copper, lead and cadmium, in red, blue, gray, rispectively
bubble(meuse, "copper", main="Copper concentration", col="red")
bubble(meuse, "lead", main="Lead concentration", col="blue")
bubble(meuse, "cadmium", main="Cadmium concentration", col="gray")


# We create an excel file covid-19 aggregated in a data. we inserted our data set in the R console and work with our dataset for the first time
# there are three steps to do:  1# download covid_agg.csv from our IOL site from MECF course. 2# create a folder called lab into C:(hard disk) 3# put the covid_agg.csv files into the lab folder
# setting the working directory: lab. remember lab folder without capital letter, R is case sensitive
# setwd function= set the working directory. for Windows users: setwd("C:/lab/") # Mac users: setwd("/Users/yourname/lab/") # Linux users: setwd("~/lab")

setwd("C:/lab/")
                                                           # recalling from lab folder our data set aggregated covid_agg.csv
                                                           # read.table function: reads a file in table format and creates a data frame from it, with cases corresponding to lines and variables to fields in the file
covid <- read.table ("covid_agg.csv", head= T)             # vector name = covid, read.table sp function, covid_agg.csv excel data frame csv is the data format, head= header T=TRUE is true that the first line is the header

head(covid)                                                # have a look at the data frame variables, category, country, cases, lat, lon

attach(covid)                                              # attach covid to the R search path                                    
                                                          
plot(country, cases)                                       # plot only the pair variables country and cases
                                                           # plot (country$country,covid$cases) if you don't attached the covid data set. $ is a good idea when you need to likage objects
                                                           # $ Operators acting on vectors, matrices, arrays and lists to extract or replace parts
                                                           # las sp package: a package providing classes and methods for spatial data: points, lines, polygons and grids. This package provides S4 classes for importing, manipulating and exporting spatial data in R, and for methods including print/show, plot, subset, [, [[, \$, names, dim, summary, and a number of methods specific to spatial data handling
plot(country, cases, las=1)                                # horizzontal labels: las=1

plot(country, cases, las=2)                                # perpendicolar labels: las=2

plot(country, cases, las=0)                                # parellel labels: las=0

plot(country, cases, las=3)                                # vertical labels: las=3

plot(country, cases, las=3, cex.axis=0.5)                  # parellel labels + exaggeration cex.axis to make the trait size (number or names) smaller in axis printed variables

plot(country, cases, las=3, cex.axis=0.7)                  # parellel labels with cex= 0.7 very well balanced trait size for our graphic
                                                           # let's make the spatial plot with ggplot2 packages: ggplot2 book, suggest book for elegant graphics maker in R
                                                           # ggplot2 is a library very usefull for graphs. Create Elegant Data Visualisations Using the Grammar of Graphics
                                                           # ggplot2 packages: Description: A system for 'declaratively' creating graphics, based on "The Grammar of Graphics". You provide the data, tell 'ggplot2' how to map variables to aesthetics, what graphical primitives to use,and it takes care of the details.
                                                           # ggplot2: to make a good graph you need 3 component: 1' the data set, 2' aestetic mapping, 3' tipe of symbol i.g. point
                                                           # let's install.packages ggplot2, you need installed packages "sp"
install.packages("ggplot2")

library(ggplot2)                                            # require ggplot2
                                                            # if ggplot2 not work, you cannot add the coastline at the covid map, following this alternative installation
                                                            # install.packages("devtools")
                                                            # devtools::install_github("tidyverse/ggplot2")
                                                            # save the .RData under the menu file, save with name. In alternative save(list = ls(all.names = TRUE), file = ".RData", envir = .GlobalEnv). It is also what happens with q("yes")
q()                                                         # the function quit or its alias q terminate the current R session

###################### Second Part of R_code_spatial
                                                            # load function: Reload datasets written with the function save. Load the previously saved R work R_code_spatial.RData
                                                            # setting the working directory in lab
setwd ("c:/lab/")
                                                            # load the previews script "R_code_spatial.RData" to recreate the previusly R global environment with all the objects
load("C:/lab/R_code_spatial.RData")
                                                            # have a look of the list of the data sets in my R environment after the loading
ls                                                          # covid
                                                            # if ggplot2 not work, you cannot add the coastline at the covid map, following this alternative installation # install.packages("devtools") # devtools::install_github("tidyverse/ggplot2")
library(ggplot2)                                            # required ggplot2
                                                            # load data set "mpg" present in ggplot2 packages: Fuel economy data from 1999 to 2008 for 38 popular models of cars. This dataset contains a subset of the fuel economy data that the EPA makes available on http://fueleconomy.gov. It contains only models which had a new release every year between 1999 and 2008 - this was used as a proxy for the popularity of the car
data(mpg)
head(mpg)                                                   # ggplot function: ggplot() initializes a ggplot object. It can be used to declare the input data frame for a graphic and to specify the set of plot aesthetics intended to be common throughout all subsequent layers unless specifically overridden
                                                            # key component: data, aesthetic map, geometry. ggplot function
ggplot(mpg, aes(x=displ, y=hwy)) + geom_point()             # ggplot funct., aesthetic mapping x= displ(engine displacement, in litres) y= hwy (highway miles per gallon), geometry of the graph= points
                                                            # plot with points data rappresentation
ggplot(mpg, aes(x=displ, y=hwy)) + geom_line()              # plot with lines data rappresentation
                                                            # plot with poligones data rappresentation
ggplot(mpg, aes(x=displ, y=hwy)) + geom_polygon()
                                                            # let's graphy covid(dataset) with the ggplot function, putting x=lat and y=lon in funcion of the variable cases and country(aes), and the geom_point(geometry) = plot with points data rappresentation
head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point()
ggplot(covid, aes(x=lon, y=lat, size=cases, col= "red")) + geom_point() # giving a different colour= "red" at the objects(points) in my plot

############################################################################################################################
############################################################################################################################

# 4. R_code_multivar.r

# R code for multivariate analysis

setwd ("C:/lab/")

install.packages("vegan")
library (vegan)                                                          #packages for vegetation analysis installed

biomes <- read.table ("biomes.csv", header=TRUE, sep = ",")              #give a name, linkage, give the command to read exel or grid dataset files, give the header name "first line" of the collums "intestazione", in the end the separator of the collums = comma
                                                                         #biomes   <- read.table     biomes.csv csv #tipe of the data file   header = T #first line of my table                                  sep = ","          
head (biomes)                                                            #have a look at the dataset

                                                                         #for multivariate analysis in many programs is a little difficult built multivariate function, in R scientist already create a function to facilities it, names DECORANA
                                                                         #DEtreanded COrrespondence ANAlysis = DECORANA
multivar <- decorana (biomes)

plot (multivar)

                                                                         #plot with exageration
plot (multivar, cex=1.2)

plot(multivar)

multivar                                                                 #we can see in the R macchine the percentage of the reality in 2 dimentions. our algorithm give us the 82% of our data just in 2 dimantion, this is not bad 

#               DCA1   DCA2    DCA3    DCA4
#Eigenvalues     0.5117 0.3036 0.12125 0.14267
#Decorana values 0.5360 0.2869 0.08136 0.04814
#Axis lengths    3.7004 3.1166 1.30055 1.47888

                                                                         #let to link same biomes for the next plot
biomes_types <- read.table ("biomes_types.csv", header=TRUE, sep = ",")
head ("biomes_types")                                                    #biomes_types= new data set
                                                                         #attach the dataset 

attach (biomes_types)

ordiellipse (multivar, type, col = 1:4, kind = "ehull", lwd = 3)         #to order inside an ellipse our biomes, the variables, collum: types, color from 1 to 4,form of the shape "ehull",  size of the line: 3

ordispider(multivar, type, col=1:4, label = T)                           #we use ordispider to connect the species at the ellipse, this can be a way to add a dimention at our graph

############################################################################################################################
############################################################################################################################

# 5. R_code_remote sensing.r

# R code for  Remote Sensing 

setwd ("C:/lab/")

library (raster)                                                         # install.packages ("raster")

install.packages("RStoolbox")                                            #packages to make our analysis, the book that the Prof. wrote is about the algorith of this R function that he created


library (RStoolbox)                                                      #to be faster in the istallation of pkgs: install.packages(c("raster", "RStoolbox")

                                                                         #load images from the lab folder
p224r63_2011 <- brick ("p224r63_2011_masked.grd")
plot(p224r63_2011)

#landsat resolution of 30m (each pixel)...B1 blu, B2 green, B3 red, B4 NIR ecc
# B1: blue
# B2: green
# B3: red
# B4: near infrared (nir) 
# B5: medium infrared
# B6: thermal infrared
# B7: medium infrared

ls()

cl <- colorRampPalette(c('black','grey','light grey'))(100)               #changing our scaled colors for graphs (gray)

plot(p224r63_2011, col=cl)                                                #Exercise: plot the image with scale color ramp palette built

                                     
cllow <- colorRampPalette(c('black','grey','light grey'))(5)              # grey scaled low amount of colours
plot(p224r63_2011, col=cllow)

names(p224r63_2011)                                                       # [1] "B1_sre" "B2_sre" "B3_sre" "B4_sre" "B5_sre" "B6_bt"  "B7_sre"

clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)          ##changing our scaled colors for graphs (blue)
plot(p224r63_2011$B1_sre, col=clb)                                        #symble $ links together collum with the data I want  

clnir <- colorRampPalette(c('red','orange','yellow'))(100)                # Exercise: plot the near infrared band with:
plot(p224r63_2011$B4_sre, col=clnir)                                      #colorRampPalette color that go from red, to orange, until yellow

dev.off ()

                                                                          # multiframe to plot more graphs together, we say at the console how many row and collums we want in the Plot Space: par(mfrow=c(2,2)) 
par(mfrow=c(2,2)) 

                                                                          # B1 blue band, we want set colors intensity in blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)                                        # func. par, multivariateframe  =  col = c finally names the colors you would 
                                                                          # plot the immages with colors calls clb 
                                                                          #dev.off to close all the windows if I need

                                                                          # Exercise: do the same for the B2_sre green band
clg <- colorRampPalette(c('dark green','green','light green'))(100)
plot(p224r63_2011$B1_sre, col=clg)

                                                                          #B3_sre red band
clr <- colorRampPalette(c('dark red','red','pink'))(100)
plot(p224r63_2011$B1_sre, col=clr)

                                                                          #B4_sre near infrared band
cln <- colorRampPalette(c('red','orange','yellow'))(100)
plot(p224r63_2011$B1_sre, col=cln)

                                                                          #change the layout of our graph in 4 rows and 1 collum and to make it different and stretched, could be really usefull for data analysis

dev.off ()
par(mfrow=c(4,1))

# B1: blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) 
plot(p224r63_2011$B1_sre, col=clb)

# B2 green
# Exercise: do the same for the green band B2_sre
clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)

# B3 red
clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)

# B4 NIR
cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)

# natural colours
# 3 components: R G B
# 3 bands: R = red band, G = green band, B = blue band
# plotRGB(p224r63_2011,r=3,g=2,b=1)

# B1: blue - 1
# B2: green - 2
# B3: red - 3
# B4: near infrared (nir) - 4 

dev.off ()                                                                  # close others plots

                                                                            # plotRGB function RGB colors

plotRGB (p224r63_2011, r=3, g=2, b=1, stretch="Lin")
                                                                            # stretching as much possible the colors

                                                                            # substitute red component with our NIR # in this way you can discriminete better the vegetation in the images, dark green is the forest, red and pink in general the pitch 
plotRGB (p224r63_2011, r=4, g=3, b=2, stretch="Lin")

                                                                            # EX. NIR on top of the G component # NIR false colours
plotRGB (p224r63_2011, r=3, g=4, b=2, stretch="Lin")
                                                                            # dark forest is forest with high amount of water

                                                                            # yellow area are bare in this way to set the colors
plotRGB (p224r63_2011, r=3, g=2, b=4, stretch="Lin")

################ day 2

#setwd("~/lab/") # linux
setwd("C:/lab/") # windows
# setwd("/Users/nome/Desktop/lab") # mac

load("R_code_remote_sensing.RData")

library(raster)                                                             # list
library (RStoolbox)

ls()

p224r63_1988 <- brick("p224r63_1988_masked.grd")

plot(p224r63_1988)

p224r63_2011 <- brick("p224r63_2011_masked.grd")

plot(p22r63_2011)
                                                                            # Exercise: plot the image using the nir on the "r" component in the RGB space
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")

                                                                            # Exrecise: plot in visible RGB 321 both images
par(mfrow=c(2,1))
  plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
  plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")
  
# plot together two images, 1988 e 2011
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin", main="1988")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin", main="2011")

# plot it in the way to see water, so with different colors relate at the NIR= water= 0
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")


# enhance the noise!
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist")
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist")
                                                                           # spectral indices
                                                                           # dvi1988 = nir1988-red1988

dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre
plot(dvi1988)

                                                                           # Exercise: calculate dvi for 2011
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre
plot(dvi2011)

dev.off()

                                                                           #giving colors at the div looking for the best grafic information
cl <- colorRampPalette(c("darkorchid3","light blue","lightpink4"))(100) 
plot(dvi2011, col=cl)

#or

cldvi <- colorRampPalette(c("light blue","light green","green"))(100)      # not the best for my purpose
plot(dvi2011, col=cldvi)

#or

cl2 <- colorRampPalette(c("yellow",'light blue','lightpink4'))(100)
plot(dvi2011,col=cl2)

                                                                           # Exercise: dvi for 1988

par(mfrow= c(2,1))
plot(dvi1988, col=cl, ylab="1988")
plot(dvi2011, col=cl, ylab="2011")


par(mfrow= c(2,1))
                                                                           # differances between 1988-2011 
plot(dvi1988,col=cl3)
plot(dvi2011, col=cl)

dev.off()

diff <- dvi2011 - dvi1988
plot (diff)

                                                                            # aggregate pixels 
                                                                            # changing the grain and resampling with factor=10 and 100 to minimize the weight
p224r63_2011res10 <- aggregate(p224r63_2011, fact=10)
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)
                                                                            #EX. plot all together images with diff grain
par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res10, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin") 

                                                                            #information about image p224r63_2011
p224r63_2011

############################################################################################################################
############################################################################################################################

# 6. R_code_PCA_remote_sensing

setwd("C:/lab/")

    #Band 1 Visible (0.45 - 0.52 µm) 30 m b
    #Band 2 Visible (0.52 - 0.60 µm) 30 m g
    #Band 3 Visible (0.63 - 0.69 µm) 30 m r
    #Band 4 Near-Infrared (0.77 - 0.90 µm) 30 m
    #Band 5 Near-Infrared (1.55 - 1.75 µm) 30 m
    #Band 6 Thermal (10.40 - 12.50 µm) 60 m Low Gain / High Gain, hydrologic analysis
    #Band 7 Mid-Infrared (2.08 - 2.35 µm) 30 m
    #Band 8 Panchromatic (PAN) (0.52 - 0.90 µm) 15 m

library (raster)
library (RStoolbox)
library (ggplot2)

p224r63_2011 <- brick("p224r63_2011_masked.grd")                    # brick import all the data of satellite image. the extantion is grd file, graphicic file
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")
ggRGB (p224r63_2011, 5, 4, 3)                                       #show the data with ggplot2


p224r63_1988 <- brick("p224r63_1988_masked.grd")                    # EX do the same, with 1988 image!
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
ggRGB (p224r63_1988, 5, 4, 3)

par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

names(p224r63_2011)

dev.off()                                                           # bet correlation between B1 and B3  # correlarion coefficient= R+= 1  Rneutral = 0  R-=-1


plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)                      # in vegetation and in our example for 2011 img. B1 and B3 have R=0.9 so very high correlation

                                                                    #PCA analysis!
                                                                    # resolution
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)
 
                                                                    # RStoolbox is now needed
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)
plot(p224r63_2011_pca$map)                                          # plot img linked map
                                                                    # coovariance matrix analysis= how much one component is far from the other 

cl <- colorRampPalette(c('dark grey','grey','light grey'))(100)     # good color for our graphic analysis
plot (p224r63_2011_pca$map, col=cl)                                   # in Principal Component analysis we have three component parts= call, model, map with $ we can link imag at the function

summary(p224r63_2011_pca$model)                                     # have a look how much is it the correlation
                                                                    # we can now see the cumulative proportion= in this case is huge correlate 99.83 
                                                                    
pairs(p224r63_2011)                                                 # with function pairs we can see the correlation between the two PCs in a graphic matrix
                                                                    # B1,B2,B3 are very correlate so the rs analysis put them in one PC
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin")

p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res) 
plot(p224r63_1988_pca$map, col=cl)                                  # we want decrease the number of variances for a 10 factor 

summary(p224r63_1988_pca$model)                                     # 99.56% of correlation in 1988 image too
                                                                
pairs(p224r63_1988)                                                 # seeing it in a graphic matric with pairs we see less correlation comparising at 2011 but very high the same!! 
                                                                    # questa cosa sarà troppo importante per vedere i cambiamenti in matematico modo
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(difpca)                                                        # differences amoung components
plot(difpca$PC1)                                                    #the PC1 have the 99% of the component
                                                                    
cldif <- colorRampPalette(c('blue','black','yellow'))(100)          #final plot:
plot(difpca$PC1,col=cldif)                                          #difpca linked with pc1 
                                                                    #with this function we can see the highest dvi ration between dvi of 2011 and 1988

plot(difpca$PC1)

############################################################################################################################
############################################################################################################################


#############################################################################################################################
#############################################################################################################################


#############################################################################################################################
#############################################################################################################################


#############################################################################################################################
#############################################################################################################################


#############################################################################################################################
#############################################################################################################################


#############################################################################################################################
#############################################################################################################################



#############################################################################################################################
#############################################################################################################################



#############################################################################################################################
#############################################################################################################################



#############################################################################################################################
#############################################################################################################################




#############################################################################################################################
#############################################################################################################################




#############################################################################################################################
#############################################################################################################################



#############################################################################################################################
#############################################################################################################################









