# R_code_exam.r

#############################################################################################################################
#############################################################################################################################

# 1. R_code_first.r
# 2. R_code_multipanel.r
# 3. R code spatial.r 
# 4. R_code_point_pattern_analysis.r 
# 5. R_code_multivariate_analysis.r
# 6. R_code_remote_sensing.r
# 7. R_code_ecosystem_function.r
# 8. R_code_remote_sensing_multivariate_analysis.r
# 9. R_code_ecosystem's_reflectance.r
# 10. R_code_faPAR.r
# 11. R_code_EBVs.r
# 12. R_code_snow.r
# 13. R_code_monitoring_air_pollution_no2.r
# 14. R_code_crop_image.r
# 15. R_code_interpolation.r
# 16. R_code_sdm.r
# 17. R_code_myproject_exam.r

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
                                                            # if ggplot2 not work, you can try to install the required packages ggplo2 in the following way: 
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
ls()                                                        # ls function: list of objects in this case covid. ls and objects return a vector of character strings giving the names of the objects in the specified environment. When invoked with no argument at the top level prompt, ls shows what data sets and functions a user has defined. When invoked with no argument inside a function, ls returns the names of the function's local variables: this is useful in conjunction with browser
                                                            # if ggplot2 not work, you can try to install the required packages ggplo2 in the following way: # install.packages("devtools") # devtools::install_github("tidyverse/ggplot2")
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

# 4. R_code_point_pattern_analysis.r 

# R_code_point_pattern_analysis: Density analysis using covid-19 number of cases per nation

install.packages("spatstat")                                          # spatstat packages: Spatial Point Pattern Analysis, Model-Fitting, Simulation, Tests: 
library(spatstat)                                                     # spatstat pack. description: Comprehensive open-source toolbox for analysing Spatial Point Patterns. Focused mainly on two-dimensional point patterns, including multitype/marked points, in any spatial region. Also supports three-dimensional point patterns, space-time point patterns in any number of dimensions, point patterns on a linear network, and patterns of other geometrical objects. Supports spatial covariate data such as pixel images. Contains over 2000 functions for plotting spatial data, exploratory data analysis, model-fitting, simulation, spatial sampling, model diagnostics, and formal inference 

setwd("C:/lab/")

covid <- read.table("covid_agg.csv", head=T)                          # import and read covid data set
 
attach(covid)
head(covid)
                                                                       # set the coordinates of the vectors in covid in relation to the coordinate reference system of the global map. For the lat and lon we used min and max values range
covids <- ppp(lon, lat, c(-180,180), c(-90,90))                        # "c" is used to clastering all the variables/numeber together. i.g. dead <- 12,34,55,66,77,88,89. This is a generic function which combines its arguments. The default method combines its arguments to form a vector. All arguments are coerced to a common type which is the type of the returned value, and all attributes except names are removed
                                                                       # ppp means panel point pattern. The general form is ppp(x.coordinates, y.coordinates, x.range, y.range) it creates a point pattern dataset in the two-dimensional plane for the range of values of lat and long.
                                                                       # attach covid at the coordinate : covids <- ppp(covid$lon, covid$lat, c(-180,180), c(-90,90)) alternative way without attach covid dataset to the R search path

d <- density(covids)                                                   # create a variable names d = density
                                                                       # density function: compute a kernel smoothed intensity function from a point pattern
plot (d)                                                               # plot this density map
                                                                       # points function: Add Points to a Plot
points(covids)                                                         # add points at the density map plot

#################################################### second part

setwd("C:/lab/")                                                       # set work directory

load("point_pattern_analysis.RData")                                   # load previously saved work from the lab folder ("point_pattern_analysis.RData")

library(spatstat)                                                      # call libraries requested
library(rgdal)                                                         # costline cannot uploaded without rgdal library

ls()                                                                   # list of objects I have in my R environment
                                                                       # plot density + adding points georeferantiated for the covids vector
plot(d)
points(covids)
                                                                       # import in lab folder and then in R "ne_10m_coastline.shp" image. shp is a data format named shape file format
                                                                       # readORG rgdal packages function: read OGR vector maps into Spatial objects. the function reads an OGR data source and layer into a suitable Spatial vector object. It can only handle layers with conformable geometry features (not mixtures of points, lines, or polygons in a single layer). It will set the spatial reference system if the layer has such metadata
                                                                       # let’s input vector lines (x0y0, x1y1, x2y2..)
coastlines <- readOGR("ne_10m_coastline.shp")                          # reading the shape file and giving a vector to the global coastlines image necessary to be added at the plot of the covid density map
                                                                       # install additional packages, alternative way to plot the coastlines at the covid density map
                                                                       # install.packages("rnaturalearth")
                                                                       # coastlines <- rnaturalearth::ne_download(scale = 10, type = 'coastline', category = 'physical')
plot(d)
points(covids)
                                                                       # plot density maps (covid-19 num. of cases per country) + add=TRUE the world coastlines map at the plot of density map and covids points
plot(coastlines, add=T)
                                                                       # colour interpolation: we can create a vector for our colour scale. 
cl <-colorRampPalette(c("yellow","orange","red")) (100)                # change the colour scale and make the graph more elegant with colorRampPalette func. We build a new object names "cl". (100) are the intermediate number of colours
                                                                       # colorRampPalette function: These functions return functions that interpolate a set of given colors to create new color palettes and color ramps, functions that map the interval [0, 1] to colors (like grey)
                                                                       
plot(d, col=cl, main="Densities of covid-19")                          # plot the density map with the new colours scale and add at the plot a label
points(covids)
plot(coastlines, add=T)                                                # add=T means that is true that we are adding the coastlines vector at the plot

                                                                       # export your Density map of covid-19 plot in pdf. pdf function: pdf starts the graphics device driver for producing PDF graphics
pdf("covid_density.pdf")

                                                                       # number of colours: abrupt change of colours!!! making just an example
cll <- colorRampPalette(c("light green", "yellow","orange","violet")) (5)
plot(d, col=cll, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

                                                                       # Exercise: create a new colourRampPalette colours scale for the same density map and analysis
clr <-colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

pdf("covid_density2.pdf")                                              # higher number of intermediate colours (1000) to have more colours quality. The best option is a balance between the heavyness of the file and the quality of the colours
clrr <-colorRampPalette(c("light green", "yellow","orange","violet")) (1000)
plot(d, col=clrr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

dev.off()                                                              # dev.off function: Control Multiple Devices: these functions provide control over multiple graphics devices, with the command off, we say at the console to close all the plots in my plot environment

pdf("covid_density3.pdf")                                              # this command is a wrapper for the pdf function

png("covid_density.png")                                               # or export in jpeg data format. png function: graphics devices for BMP, JPEG, PNG and TIF format bitmap files

clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

dev.off()                                                               # to close all open plots

############################################################################################################################
############################################################################################################################

# 5. R_code_multivar.r

# R code for multivariate analysis

setwd ("C:/lab/")
                                                                         # vegan packages: Community Ecology Package: Ordination, Diversity and Dissimilarities: the vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. Most of its multivariate tools can be used for other data types as well
install.packages ("vegan")                                               # packages for vegetation analysis installed
library (vegan)                                                           

biomes <- read.table ("biomes.csv", header=TRUE, sep = ",")              # give a name at the new vector,    linkage it with data set,     command necessary to read excel or grid data files,   tell at the console that the "first line" is the "header" of the data frame,      separator of the collums = comma
                                                                         # biomes                            <-                            read.table (biomes.csv) csv=tipe of the data file     header = T #first line is the header                                              sep = ","          
head (biomes)                                                            # have a look at the head of the dataset biomes

                                                                         # multivariate statistic analysis in many softwares could be a difficult task. To build a multivariate function in a  script could be really complex, for those reasons in R, scientist created a function to facilities the multivariate analysis. The name of M.A. function in R vegan packages is DECORANA
                                                                         # DEtreanded COrrespondence ANAlysis = DECORANA. Detrended Correspondence Analysis and Basic Reciprocal Averaging. Performs detrended correspondence analysis and basic reciprocal averaging or orthogonal correspondence analysis.
multivar <- decorana (biomes)                                            # for multivariate stat. analysis of biomes data set, we need just to process biomes with DECORANA vegan pkgs function.

plot (multivar)
                                                                         # plot multivariate analysis with a little exageration
plot (multivar, cex=1.2)

multivar                                                                 # we can see in the R macchine the percentage of the reality in 2 dimentions. Our algorithm give us the 82% of our data set rappresented just in 2 dimantion. This is a good rappresentation of the real data set information reducing as much as I can the real dimentions, loosing only the 18 % of our initial information
#               DCA1   DCA2    DCA3    DCA4
#Eigenvalues     0.5117 0.3036 0.12125 0.14267
#Decorana values 0.5360 0.2869 0.08136 0.04814
#Axis lengths    3.7004 3.1166 1.30055 1.47888

                                                                         # let to link same biomes types together to see correlations among the species distribution in my multivariate analysis rappresentation 
biomes_types <- read.table ("biomes_types.csv", header=TRUE, sep = ",")  # import a new data set with bioms_types. type= temperate, tropical, conifer forest, boreal
head ("biomes_types")                                                    # biomes_types = new data set
                                                                         # attach the dataset at the R search path
attach (biomes_types)

ordiellipse (multivar, type, col = 1:4, kind = "ehull", lwd = 3)         # to order inside an ellipse our biomes, the variables = names of the species associated at the biomes(multivar), collum interested as second data set = types = second variable biomes types, colours from 1 to 4, form of the shape "ehull",  size of the line: 3
                                                                         # ordiellipse and ordispider function: Display Groups or Factor Levels in Ordination Diagrams. description: functions to add convex hulls, “spider” graphs, ellipses or cluster dendrogram to ordination diagrams. The ordination diagrams can be produced by vegan plot.cca, plot.decorana or ordiplot.
ordispider (multivar, type, col=1:4, label = T)                          # we use ordispider to connect the species at the ellipse, this can be a way to add a rappresentation of another dimention at our graphical analysis

############################################################################################################################
############################################################################################################################

# 6. R_code_remote_sensing.r

# R code for  Remote Sensing 

setwd ("C:/lab/")                                                        # install.packages ("raster")

install.packages("RStoolbox")                                            # packages to make our analysis described into the professor's book. 
                                                                         # RStoolbox packages: Tools for Remote Sensing Data Analysis: description: Toolbox for remote sensing image processing and analysis such as calculating spectral indices, principal component transformation, unsupervised and supervised classification or fractional cover analysis
library (raster) 
library (RStoolbox)                                                      # to be faster in the installation of pkgs: install.packages(c("raster", "RStoolbox")

                                                                         # load image "p224r63_2011_masked.grd" from the lab folder
p224r63_2011 <- brick ("p224r63_2011_masked.grd")                        # brick RStoolbox function: Create a RasterBrick object: a RasterBrick is a multi-layer raster object. They are typically created from a multi-layer (band) file; but they can also exist entirely in memory. They are similar to a RasterStack (that can be created with stack), but processing time should be shorter when using a RasterBrick. Yet they are less flexible as they can only point to a single file. A RasterBrick can be created from RasterLayer objects, from a RasterStack, or from a (multi-layer) file. The can also be created from SpatialPixels*, SpatialGrid*, and Extent objects, and from a three-dimensional array
plot(p224r63_2011)

#landsat resolution of 30m (each pixel)...B1 blu, B2 green, B3 red, B4 NIR ecc
# B1: blue
# B2: green
# B3: red
# B4: near infrared (nir) 
# B5: medium infrared
# B6: thermal infrared
# B7: medium infrared

ls()                                                                      # list of objects

cl <- colorRampPalette(c('black','grey','light grey'))(100)               # changing our colours scale for the plot (gray)

plot(p224r63_2011, col=cl)                                                # Exercise: plot the image with colorRampPalette scale built = cl

                                     
cllow <- colorRampPalette(c('black','grey','light grey'))(5)              # grey scaled lower amount of intermediate colours = 5
plot(p224r63_2011, col=cllow)
                                                                          # names function: functions to get or set the names of an object
names(p224r63_2011)                                                       # [1] "B1_sre" "B2_sre" "B3_sre" "B4_sre" "B5_sre" "B6_bt"  "B7_sre"

clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)          # changing our colours scale for graphs in blue
plot(p224r63_2011$B1_sre, col=clb)                                        # symbol $ links together collumns with the data set I want, in this case the light band 1 (blue)
                                                                          # plot image with band layer that I decided to extract and focus on the image
clnir <- colorRampPalette(c('red','orange','yellow'))(100)                # Exercise: plot the near infrared band with:
plot(p224r63_2011$B4_sre, col=clnir)                                      # colorRampPalette colours that go from red, to orange, until yellow. symbol $ links together collumns with the data set I want, in this case the light band 4 (near infrared)

dev.off ()
                                                                          # multiframe to graphy more plots together, we say at the console how many rows and collums we want in the Plot Space: par(mfrow=c(2,2)) so 2x2 plots.
par(mfrow=c(2,2))                                                         # par (parameters) function: Set or Query Graphical Parameters: par can be used to set or query graphical parameters. Parameters can be set by specifying them as arguments to par in tag = value form, or by passing them as a list of tagged values
                                                                          # B1 blue band, we want set colors intensity in blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)                                        # func. par, multivariateframe  =  col = clb finally names the colors you prefer
                                                                          # plot the immages with colors calls clb 

                                                                          # Exercise: do the same for the B2_sre green light band
clg <- colorRampPalette(c('dark green','green','light green'))(100)
plot(p224r63_2011$B1_sre, col=clg)
                                                                          # B3_sre red light band
clr <- colorRampPalette(c('dark red','red','pink'))(100)
plot(p224r63_2011$B1_sre, col=clr)
                                                                          # B4_sre near infrared light band
cln <- colorRampPalette(c('red','orange','yellow'))(100)
plot(p224r63_2011$B1_sre, col=cln)
                                                                          # change the layout of our graph in 4 rows and 1 collum, could be really usefull for data analysis have parallel plots

dev.off ()

par(mfrow=c(4,1))                                                         # multiframe of 4x1 plots matrix so we could have 4 rows in a single collumn

# B1: blue
#clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) 
plot(p224r63_2011$B1_sre, col=clb)

# B2 green
# Exercise: do the same for the green band B2_sre
#clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)

# B3 red
#clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)

# B4 NIR
#cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)

# primary natural colours # Make a Red-Green-Blue plot based on three layers (in a RasterBrick or RasterStack). Three layers (sometimes referred to as "bands" because they may represent different bandwidths in the electromagnetic spectrum) are combined such that they represent the red, green and blue channel. This function can be used to make 'true (or false) color images' from Landsat and other multi-band satellite images
# 3 components: R G B
# 3 bands: R = red band, G = green band, B = blue band. r integer index of the Red channel, g integer index of the Green channel, b integer index of the Blue channel.
# plotRGB(p224r63_2011,r=3,g=2,b=1)

# B1: blue - 1
# B2: green - 2
# B3: red - 3
# B4: near infrared (nir) - 4 

dev.off ()                                                                  # close plots

                                                                            # RGB colours: The RGB function describes a color giving the intensity of the 3 primary colors: red, green and blue. This function creates colors of class RGB
                                                                            # plotRGB raster packages function: Red-Green-Blue plot of a multi-layered Raster object: make a Red-Green-Blue plot based on three layers (in a RasterBrick or RasterStack). Three layers (sometimes referred to as "bands" because they may represent different bandwidths in the electromagnetic spectrum) are combined such that they represent the red, green and blue channel. This function can be used to make 'true (or false) color images' from Landsat and other multi-band satellite images
plotRGB (p224r63_2011, r=3, g=2, b=1, stretch="Lin")                        # stretch function: Linear stretch of values in a Raster object. Provide the desired output range (minv and maxv) and the lower and upper bounds in the original data, either as quintiles (if minq=0 and maxq=1 you use the minimum and maximum cell values), or as actual values (smin and smax; e.g. precomputed quantile values). If smin and smax are both not NA, minq and maxq are ignored
                                                                            # normal band layers, red=3, green=2, blue= 1. stretch= stretching as much possible the colours

                                                                            # substitute red component with our NIR band layer # in this way you can discriminete better the vegetation in the images, dark green is the forest, red and pink in general the agriculture fields.
plotRGB (p224r63_2011, r=4, g=3, b=2, stretch="Lin")                        # putting NIR light band on top of the image layers, red at the place of green and green at the place of blue

                                                                            # EX. NIR on top of the G component # NIR false colours
plotRGB (p224r63_2011, r=3, g=4, b=2, stretch="Lin")
                                                                            # dark forest is forest with high amount of water. This time we put at the top layer the red light band reflectance, later the NIR and for the last the green band

                                                                            # yellow areas are the not vegetated areas or poorly vegetated in this chromatic way to set the colours
plotRGB (p224r63_2011, r=3, g=2, b=4, stretch="Lin")

################  second part 

#setwd("~/lab/") # linux
setwd("C:/lab/") # windows
# setwd("/Users/nome/Desktop/lab") # mac

load("R_code_remote_sensing.RData")                                         # load the saved "R_code_remote_sensing.RData"

library(raster)                                                             
library (RStoolbox)

ls()                                                                        # list of objects

p224r63_1988 <- brick("p224r63_1988_masked.grd")                            # brick the image ("p224r63_1988_masked.grd") of the same place, but 20 years older 1988-2011

plot(p224r63_1988)

p224r63_2011 <- brick("p224r63_2011_masked.grd")                            # making a comparation in the vegetated area in 1988 and in 2011 in a specific amazon forest locality

plot(p224r63_2011)
                                                                            # Exercise: plot the image using the "nir" on the "r" component in the RGB space
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")

                                                                            # Exrecise: plot in visible RGB space r=3 g=2 b=1, both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")                         # multiframe par function of 2 parallel plots, 1988 image vs 2011 image, to have a clear evidence of eventual land cover modifications
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")
  
# plot together two images, 1988 e 2011                                     # making the same multiframe analysis, with the plotRGB function, primary colours, RGB space but putting the Near InfraRed band layer over the Red band layer and the Red over the Green band layer
#par(mfrow=c(2,1))                                                          # putting the labels at both the plots
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin", main="1988")            # r=4=nir; g=3=red; b=2=green layers
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin", main="2011")

# plot it in the way to see water, so with different colours relate at the NIR= water= 0 reflectance
#par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")                         # we want have the same graphical analysis and RGB layers order of the previous multiframe, but without labels to have the maximum image extention
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")


                                                                           # enhance the noise, probably high ammount of water in the air, stretching these images with histogram equalization of band information result to lose the images quality               
#par(mfrow=c(2,1))                                                         # hist: calcutation of the area under the integral that describes the shock of enhancing the colour
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist")                       # same multiframe and RGB order of layers. 
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist")                       # stretch="hist": character. Option to stretch the values to increase the contrast of the image: "lin" or "hist". Linear equalization or Histogram equalization.
                                                                           # spectral indices
dev.off()

# to calculate DVI of 1988                                                 # NDVI: Calculation of the Normalized Difference Vegetation Index (NDVI). Chlorophyll (a health indicator) strongly absorbs visible light, and the cellular structure of the leaves strongly reflect near-infrared light. When the plant becomes dehydrated, sick, afflicted with disease, etc., the spongy layer deteriorates, and the plant absorbs more of the near-infrared light, rather than reflecting it. Thus, observing how NIR changes compared to red light provides an accurate indication of the presence of chlorophyll, which correlates with plant health
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre                       # mesuring the NDVI: Normalized Difference Vegetation Index. In Landsat 4-7, NDVI = (Band 4 – Band 3) / (Band 4 + Band 3). I.G. DVI1988 = (nir1988-red1988)/(nir1988+red1988). Higher the index and healthier should be the vegetation
plot(dvi1988)                                                              # plot the normalized difference vegetation index (NDVI) values derivated from the previous band NIR-Red analysis
                                                                           # we linked together image with the band interested for our function. In this way we isolated values of the 1988 image for the NIR and Red wavelenght bands, neccessary for the DVI equation
# to calculate DVI of 2011                                                 # Exercise: calculate the DVI in 2011 image
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre
plot(dvi2011)

diff <- dvi2011 - dvi1988                                                  # consider the difference between DVI2011 and DVI1988                             
plot (diff)

dev.off()
                                                                           # trying to give colours contrast and gradation scale at the DVI images, looking for the best grafic information
cl <- colorRampPalette(c("darkorchid3","light blue","lightpink4"))(100) 
plot(dvi2011, col=cl)

#or

cldvi <- colorRampPalette(c("light blue","light green","green"))(100)      # not the best for my purpose
plot(dvi2011, col=cldvi)

#or

cl2 <- colorRampPalette(c("yellow",'light blue','red'))(100)              # well balanced colorRampPalette colours vector
plot(dvi2011,col=cl2)

                                                                          # making two final comparisons in a multiframe of 2x1 plot space

par(mfrow= c(2,1))
plot(dvi1988, col=cl, ylab="1988")
plot(dvi2011, col=cl, ylab="2011")


par(mfrow= c(2,1))                                                        # differances between 1988-2011 without labels. Yellow enhance bad vegetation health or cover, or high ammount of water (rivers)
plot(dvi1988,col=cl2)
plot(dvi2011, col=cl2)

dev.off()
                                                                          # aggregation of pixels. "aggregate" sp function: spatial aggregation of thematic information in spatial objects
                                                                          # changing the grains (dimention of pixel), aggregating them with factor=10 and =100. 1 grain that rappresent 10 or 100 pixels instead than 1. In this way you can obtain a file lighter
p224r63_2011res10 <- aggregate(p224r63_2011, fact=10)
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)
                                                                          # Exercise plot all together images with diff. grains aggregation. In this way you will understand positive and negative aspects of using this procedure
par(mfrow=c(3,1))                           
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res10, r=4, g=3, b=2, stretch="Lin")                  # sometime is a good idea reduce the resolution.     DOI: 10.1007/s10531-008-9479-0
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")                 # we lost a lot the quality of the image

                                                                          # information about image p224r63_2011
p224r63_2011

############################################################################################################################
############################################################################################################################

# 7. R_code_ecosystem_functions.r

# R_code_ecosystem_functions.r

# R code to view biomass over the world and calculate changes in ecosystem functions
# energy
# chemical cycling
# proxies
                                                                                       # rasterdiv packages: Diversity Indices for Numerical Matrices: Rasterdiv basics. Derive indices of diversity from NDVI.Providing functions to calculate indices of diversity on numerical matrices based on information theory. The rationale behind the package is described in Rocchini, Marcantonio and Ricotta (2017) doi:10.1016/j.ecolind.2016.07.039
install.packages("rasterdiv")                                                          # packages to assess Rao variations in different ecosystems
library(rasterdiv)
install.packages("rasterVis")                                                          # rasterVis packages: Visualization Methods for Raster Data: Methods for enhanced visualization and interaction with raster data. It implements visual-ization methods for quantitative data and categorical data, both for univariate and multivari-ate rasters. It also provides methods to display spatiotemporal rasters, and vec-tor fields. See the website for examples
library (rasterVis)                                                                    ## Loading required package:  raster 
library(ratser)                                                                        ## Loading required package:  sp
                                                                                       ## Loading required package:  lattice
                                                                                       ## Loading required package:  latticeExtra

data(copNDVI)                                                                          # load copNDVI dataset (inside rasterdiv library): a RasterLayer (EPSG: 4326) of the global average NDVI value per pixel for the 21st of June over the period 1999-2017
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253, 255, NA), right=TRUE)                        # removing water pixels using cbind argument for 253,254, 255 values= Not Assigned
levelplot(copNDVI)                                                                     # plot the same image with pixels aggregation of factors=10 and factors=100
                                                                                       # reclassify function: Reclassify values of a Raster* object. The function (re)classifies groups of values to other values. For example, all values between 1 and 10 become 1, and all values between 11 and 15 become 2 
copNDVI10 <- aggregate(copNDVI, fact=10)                                               # levelplot function: #plot a raster object and displaying the level of the values along x and t axis. Level plots and contour plots. Draws false color level plots and contour plots
levelplot(copNDVI10)

copNDVI100 <- aggregate(copNDVI, fact=100)                                             # sometimes is better aggregate pixels for the raster analysis. High resolution satellite imagery for tropical biodiversity studies: The devil is in the detail  DOI: 10.1007/s10531-008-9479-0
levelplot(copNDVI100)

library(ggplot2)
library (RStoolbox)

myPalette <- colorRampPalette(c('white','green','dark green'))                         # scale_colour_gradientn: gradient colour scales: 
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))                 # to plot NDVI score on the map, very impressive map

ggR(copNDVI, geom_raster = TRUE)+
scale_fill_gradientn(name = "NDVI", colours = myPalette(100))+
labs(x="Longitude",y="Latitude", fill="")+
theme(legend.position = "bottom") +
  NULL +
ggtitle("NDVI")

## deforestation project

setwd("C:/lab/")
defor1 <- brick("defor1_.jpg")                                                          # to import images and link them with bands
defor2 <- brick("defor2_.jpg")                                                          # band1=NIR, Band2=red, band3=green defor1_.1 or defor2_.1 (band1)

                                                                                        # plotRGB: plot using RGB colours: Make a Red-Green-Blue plot based on three layers (in a RasterBrick or RasterStack). Three layers (sometimes referred to as "bands" because they may represent different bandwidths in the electromagnetic spectrum) are combined such that they represent the red, green and blue channel. This function can be used to make 'true (or false) color images' from Landsat and other multi-band satellite images 
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")                                           # par function to make a multiframe plot to compare the same image in two different years 1x2
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")                                           # default association of RGB colours


par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

                                                                                        # calculating the DVI for both images NIR-Red bands of the associated at the two images
                                                                                        # $ symbol to link the layers at the image, each layer rappresents a band
dvi1 <- defor1$defor1_.1 - defor1$defor1_.2                                             
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2                                             # to plot NDVI score on the map, very impressive map, we resampled the colorRampPalette coloration
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100)
par(mfrow=c(2,2))                                                                       # multiframe of plot matrix 4x4 to have a visive comparation of the two image(same place different time), adding difference in NDVI index, graphically and through an histogram (diff. DVI)
plot(dvi1, col=cl)
plot(dvi2, col=cl)
                                                                                        # make a difference between DVI of 1' image and 2' image
difdvi <- dvi1 - dvi2                                                                   # Warning in dvi1 - dvi2:  Raster objects have different extents.  Result for their intersectionis returned
cld <- colorRampPalette(c('blue','white','red'))(100) 

plot(difdvi, col=cld)                                                                   # to see the lost of ecosystem services, related at the lost of the vegetation. Healthness of the vegetation in the area is lower= DVI is lower
hist(difdvi)                                                                            # high lost in primary production ecosystem services and biomass production
                                                                                        # hist=histogram: function: The generic function hist computes a histogram of the given data values. If plot = TRUE, the resulting object of class "histogram" is plotted by plot.histogram, before it is returned. Have a look also of hist.im: Histogram of Pixel Values in an Image
sessionInfo ()                                                                          # R session information, matrix, packages

############################################################################################################################
############################################################################################################################


# 8. R_code_remote_sensing_multivariate_analysis.r
# Principal Component Analysis (PCA) in Remote Sensing

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

p224r63_2011 <- brick("p224r63_2011_masked.grd")                    # brick to import all the data of satellite image (rasterlayers related at the bands). the extantion is grd file format, graphic files
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")                 # in the RGB layers space with resampled the bands related at r=5= NIR(over 1.55 µm),  g=4= NIR(lower 0.90 µm),  b=3=Red
ggRGB (p224r63_2011, 5, 4, 3)                                       # ggRGB function in ggplot2 packages: create ggplot2 Raster Plots with RGB colours from 3 RasterLayers: calculates RGB color composite raster for plotting with ggplot2. Optional values for clipping and and stretching can be used to enhance the imagery


p224r63_1988 <- brick("p224r63_1988_masked.grd")                    # Exercise do the same procedure and plot, with 1988 image
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
ggRGB (p224r63_1988, 5, 4, 3)

par(mfrow=c(1,2))                                                   # plotting the two images together with par function(multiframe 1x2) with RGB focalized on NIR and Red bands
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")                 # is quite impressive the discontinuity of the image 2011 in plant cover and reflectance
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

names(p224r63_2011)                                                 # get the names of objects

dev.off()                                                           # betting the correlation between B1 and B3  # correlarion coefficient= R+= 1  Rneutral = 0  R-=-1

                                                                    # PCA, multivariate analysis on the two wavelenght bands blue and red
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)                      # in vegetation and in our example for 2011 img. B1 and B3 have R=0.9 so very high positive correlation of the two bands values (Blue and Red wavelenght bands), this means that the information on reflectance comes from vegetation(assorbing in a good and close ratio B1 and B3 blue and red)                                                                  
                                                                    # changing the resolution
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)
                                                                    # RStoolbox is now needed: rasterPCA function: Principal Component Analysis for Rasters: calculates R-mode PCA for RasterBricks or RasterStacks and returns a RasterBrick with multiple layers of PCA scores
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)                     # PCA, raster analysis of principal component
plot(p224r63_2011_pca$map)                                          # plot image linked to the map

                                                                    # coovariance matrix analysis= how much one component is far from another component
cl <- colorRampPalette(c('dark grey','grey','light grey'))(100)     # grey gradation for our graphic analysis
plot (p224r63_2011_pca$map, col=cl)                                 # in a principal component analysis: we have three components = call, model and map with $ we can link component at the model
                                                                    # summary function: summary is a generic function used to produce result summaries of the results of various model fitting functions. The function invokes particular methods which depend on the class of the first argument
summary(p224r63_2011_pca$model)                                    
                                                                    # have a look how much is it the correlation amoung components
                                                                    # we can now see the cumulative proportion and standard deviation of our principals component. In this case there is a huge correlation = 0.98 
                                                                    
pairs(p224r63_2011)                                                 # with function pairs we can see the correlation between the two PCs in a graphic matrix
                                                                    # B1,B2,B3 are very correlate so the rs analysis putting them almost in a single PC
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin")         # plotting in RGB plot space

p224r63_1988_res <- aggregate(p224r63_1988, fact=10)                # same PCA and remote sensing analysis, but with the 1988 image
p224r63_1988_pca <- rasterPCA(p224r63_1988_res) 
plot(p224r63_1988_pca$map, col=cl)                                  # we want decrease the number of variances for a 10 factor 

summary(p224r63_1988_pca$model)                                     # model function: Automatically selects and then provides an analysis of a linear model
                                                                
pairs(p224r63_1988)                                                 # seeing it in a graphic matric with pairs we see less correlation comparising at 2011 but the correlation amoung variables is high as well
                                                                    # very good matematic and statistical analysis to have evidance in environmnental changes
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(difpca)                                                        # differences amoung components
plot(difpca$PC1)                                                    # PCA correlation is very high in 1988 image 0.82
                                                                    
cldif <- colorRampPalette(c('blue','black','yellow'))(100)          # final plots, plotting the dvi difference
plot(difpca$PC1,col=cldif)                                          # differance between components.
                                                                    # with this function we can see best dvi differences between dvi of 2011 and 1988

plot(difpca$PC1)

############################################################################################################################
############################################################################################################################

# 9. R_code_ecosystem's_reflectance.r

# R_code_radiance

library(raster) 
                                                                       # raster function in raster packages: Create a RasterLayer object: Methods to create a RasterLayer object. RasterLayer objects can be created from scratch, a file, an Extent object, a matrix, an 'image' object, or from a Raster*, Spatial*, im (spatstat) asc, kasc (adehabitat*), grf (geoR) or kde object. In many cases, e.g. when a RasterLayer is created from a file, it does (initially) not contain any cell (pixel) values in (RAM) memory, it only has the parameters that describe the RasterLayer. You can access cell-values with getValues, extract and related functions. You can assign new values with setValues and with replacement
toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)              # create a new raster dataset: you have to specify ncollumns, number of rows, x and y minimum and maximum value
values(toy) <- c(1.13,1.44,1.55,3.4)                                   

plot(toy)                                                              # digits function: Return the digits that make up an integer. Takes an integer or vector of integers and returns a vector, list, or matrix of the individual digits (decimal) that make up that number
text(toy, digits=2)                                                    # text function: text draws the strings given in the vector labels at the coordinates given by x and y. y may be missing since xy.coords(x, y) is used for construction of the coordinates. Put the informations in the pixel
                                                                       # stretching the values for toy dataset
toy2bits <- stretch(toy,minv=0,maxv=3)                                 # 2bits= 2^2 informations = numbers. we can traslate the informations (radiance) in bits informations and the opposite (traslating not integer numbers, in integer numbers)into my dataset
                                                                       # storage.mode function: Retrieve or set storage mode for eSets. These generic functions report or change the storage mode in a dataset
storage.mode(toy2bits[]) = "integer"                                   # the way to storage our informations = integer

plot(toy2bits)
text(toy2bits, digits=2)                                               # plot the new frame with text on the pixels

toy4bits <- stretch(toy,minv=0,maxv=15)                                # we can change the information from 2bits to 4bits= 16 integers 2^4
storage.mode(toy4bits[]) = "integer"
                                                                       # storage and plot the new informations
plot(toy4bits)    
text(toy4bits, digits=2)
                                                                       # passing from 4 bits to 8 bits, 2^8= 256 integers                              
toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer"

plot(toy8bits)
text(toy8bits, digits=2)                
                                                                       # par multiframe of 1 row and 4 columns, where i can plot all the 4 examples made in the R_code_ecosystem's_reflectance                 
par(mfrow=c(1,4))  

plot(toy)
text(toy, digits=2)
plot(toy2bits)
text(toy2bits, digits=2)
plot(toy4bits)
text(toy4bits, digits=2)
plot(toy8bits)
text(toy8bits, digits=2)

#############################################################################################################################
#############################################################################################################################

# 10. R_code_faPAR.r

# how to look energy flow and CO2 cycle from satellite: faPAR fraction of absorbed photosynthetically active radiation. The faPAR vegetation index is really related at the biomass producted from the primary prodation of the vegetation. faPAR=ammount of photosynthesis.

library(raster)                                             # install.packages("raster")
library(rasterVis)                                          # needed for levelplot
library(rasterdiv)

                                                            # copDVI= copernicus DIV = World wide analysis DVI different variation index
setwd("C:/lab/")
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA))          # reclissifaty coppenNDVI   remuving data from 255-255 and putting No assigned= NA. to remuve water values
levelplot(copNDVI)                                          # levelplot: Draw Level Plots and Contour plots arround the axis
                                                            # NDVI: A RasterLayer (EPSG: 4326) of the global average NDVI value per pixel for the 21st of June over the period 1999-2017. Normalized difference vegetation index. 
faPAR10 <- raster("faPAR10.tif")                            # named faPAR10 because pixels aggregation factor of 10 ,   raster function that import one layer at time,  faPAR10= name of image
levelplot(faPAR10)                                          # levelplot(copNDVI) makes an avarage of the pixel variable value, very high NDVI are in tropical forests and temperate northern forests, NDVI index in the two biomes is similarly high 
                                                            # levelplot(faPAR10): faPAR index put in evidence the highest capability of tropical forest (much more of temperate) to use the light photosinthetically active ratio and produce biomass
dev.off()

pdf("copNDVI.pdf")                                          # make a pdf for NDVI and faPAR plot indicies
levelplot(copNDVI)
dev.off()

pdf("faPAR.pdf")
levelplot(faPAR10)

dev.off()

######################################### second part 

setwd("C:/lab/")

load("faPAR.RData")                                          # load the R_code_faPAR
                                                             # the original faPAR image from copernicus is 2gigabits
                                                             # let's see how much space needs for 8-bits image
ls()                                                         # lookin for list of object in faPAR10
faPAR10                                                      # have a look on the informations related to the raster layer
                                                             # writeRaster function in raster pkgs: Write raster data to a file. Write an entire Raster* object to a file, using one of the many supported formats. See writeValues for writing in chunks (e.g. by row)
writeRaster(copNDVI, "copNDVI.tif")                          # write the data copNDVI in .tif, 5.3MB

faPAR <- stretch(faPAR10,minv=0,maxv=250)                    # faPAR10 in bits from 0 to 0.93. Change from 0 to 255

writeRaster (faPAR, "faPAR.tif")
faPAR

levelplot(faPAR)                                             # Exercise make a levelplot for this dataset

#### third part

setwd("C:/lab/")
load("faPAR.RData")                                          # load the R_code_faPAR                                                             # regression model between faPAR and NDVI
                                                             # make before a general example with erosion and heavy metal
erosion <- c(12,14,16,24,26, 40,55,67)                       # to make a dataset of values. I.G. % of erosion or concentration of heavy metals
                                                             # heavy metal in ppm => "hm"
hm <- c(30, 100, 150, 200, 260, 340, 460, 600)

plot(erosion, hm, col="red", pch=19, xlab = "Erosion", ylab= "Heavy Metal")   # xlabel, ylabel, make a specific axis label
                                                                              # octagon character is the symble pch=19 
                                                             
model1 <- lm(hm ~ erosion)                                   # linear model function "lm", Fitting Linear Models: lm is used to fit linear models. It can be used to carry out regression, single stratum analysis of variance and analysis of covariance (although aov may provide a more convenient interface for these)
summary(model1)
                                                             # the line describe by the ab ( slope and intercept)
abline(model1)                                               # abline function: Add Straight Lines to a Plot. This function adds one or more straight lines through the current plot
                                                             # make the linear model for the relationship between biomass and ecosystem function # we want the same analysis for faPAR and NDVI vegetation indecies
#### fourth part

setwd("C:/lab/")

faPAR10 <- raster("faPAR10.tif")                             # import faPAR10 agg. img raster layer
faPAR10                                                      # 56 milions of pixels/cells
plot(faPAR10)                                                # plot fapar and ndvi
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)   # use reclassify and cbind functions to exclude values of radiation between 253-255(water)

install.packages ("sf")                                      # sf packages: Support for simple features, a standardized way to encode spatial vector data. Binds to 'GDAL' for reading and writing data, to 'GEOS' for geometrical operations, and to 'PROJ' for projection conversions and datum transformations.
library (sf)                                                 # to call st_* functions you need sf packages
                                                             # create a random point algorithm
random.points <- function(x,n)                               # algorithm to select the point on our raster image randomly
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union function: to dissolve geometries: Combine several feature geometries into one, without unioning or resolving internal boundaries
pts <- spsample(pol[1,], n, type = 'random')
}
pts <- random.points( faPAR10, 1000)                         # then select 1000 points randomly from faPAR10 image raster layer
plot(faPAR10$pts)                                            # random.points is our new function/algorithm created
                                                             # having a view of this 1000 points. we use 1000 points to rappresent the images
copNDVIp <- extract(copNDVI, pts)                            # extract function: Operators acting on vectors, matrices, arrays and lists to extract or replace parts
faPAR10p <- extract(faPAR10, pts)
                                                             # photosinthesis vs biomass production linear model
model2 <- lm(faPAR10p ~ copNDVIp)
summary(model2)                                              # having a look of the correlation of the two variables

plot(copNDVIp, faPAR10p, col="green", xlab="biomass", ylab="photosynthesis")
abline(model2, col="red")                                    # second linear model for my function, phot. and biom. are related                   

levelplot(copNDVI)
levelplot(faPAR10)          


#############################################################################################################################
#############################################################################################################################

# 11. R_code_EBVs.r

# R code Essential Biodiversity Variable
# Understanding Heterogeneity 

setwd("C:/lab/")   # set always the working directory in lab folder!
library (raster)

# ratser function can import a single layer band.
# brick func. can import all the layer bands of my file/image.
snt <- brick("sentinel.tif")

snt  # our new object

# a bits theory remind
# 2^8 # bits
#[1] 256
#> 2^12
#[1] 4096
#> 2^13
#[1] 8192

# Boa Vista Brazil (our location)
# plot this images
plot(snt)

#B1 blue
#B2 green
#B3 red
#B4 NIR

# R G B
plotRGB(snt,3,2,1, stretch="lin")
plotRGB(snt,4,3,2, stretch="lin") #putting the nir band on the top layer

# how the different layers(bands) are related... for a better understanding of PCA analysis
# standard deviation is possible to calculate only in one layer at time.
# we can do it with the multivariate analysis, in the way to have only 1 dimention instead than one PCA
# have a look of the relationship of the bands
pairs(snt)

# PCA analysis

# install.packages("RStoolbox")
library(RStoolbox) # essential for the PCA

sntpca <- rasterPCA(snt)  # to select the raster PCA in which we are interested
sntpca                    # having a look of the informations of my PCA vector

# we want link the model at our image
# information about the output of the model, we are looking for the percentage of variance related to the components
summary(sntpca$model) 
# Cumulative Proportion    0.7015076  means that we are using 70% of the original information

plot(sntpca$map)    # plot sntpca linked to the raster map

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")  # in this way we can have different colorations of the area based on the bands RGB

# set the moving window, we projected 5x5 pixels at time (# create the moving window: it is a matrix), we need that because we are mesuring the standard deviation of window
# all the values are set to 1, empty window

window <- matrix(1, nrow= 5, ncol = 5)

# focal function for sd, in this case, works only for RasterLayer
# focal function: Calculate focal ("moving window") values for the neighborhood of focal cells using a matrix of weights, perhaps in combination with a function
# function focal is using for the calculation of our SD of the moving window, see the description on CRAN for Focal Function to know mathematic analysis that we can do with focal function.

sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)         # meusuring standard deviation of the PCA in built moving windows on my raster layer map


cl <- colorRampPalette(c('dark blue','green','orange','red'))(100)  
plot(sd_snt, col=cl)

# we want a multiframe to plot RGB image and PCA sd image together

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
#plotRGB(clad, 4,2,3, stretch = "lin")
#plotRGB(clad, 3,2,1, stretch = "lin")
#plotRGB(clad, 1,3,2, stretch = "lin")
#plotRGB(clad, 2,3,1, stretch = "lin")
# the goal is to apply standard deviation to the principal component analysis of the clad image
# we have to decide our moving window, to calculate the SD and report the SD on our pixel. "1" is to set the value of the moving window and it doesn't influence the calculation.
window <- matrix(1, nrow= 3, ncol= 3)
window

# focal function have the capacity to make the calcolation of SD
cladpca <- rasterPCA(clad)
cladpca  # we can see the output(informations) for cladpca

summary(cladpca$model)
# 98% of the images informations is showing from the model

# link the cladpca with map and with PC1 # folcal function can be applied also to a image directly taken in field: cladonia.jpg. function=standard deviation
sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd)

# we want aggregate the image's pixels to make the sd analysis faster

PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)

# relate powerful colours for our function and use par to plot a frame of plots
# plot the calculation both sd and sd aggregated

par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl)                                 # less resolution of the image 

# we want see the original image and secon the variability of the image. to see the variability of each individual, so the inner variation with sd calculation
# could be a biodiversity analysis, just using the standard deviation variability of the reflectance colours bands
par(mfrow=c(1,2))
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plotRGB(clad, 1,2,3, stretch = "lin")
plot(sd_clad, col=cl)  # less accurancy

#############################################################################################################################
#############################################################################################################################

# 12. R_code_snow.r

# R_code_snow.r 

# setwd("~/lab/") #linux
# setwd("/Users/utente/lab") #mac
setwd("C:/lab/") 
                                                                      # ncdf4 packages: returns a string that is the version number of the ncdf4 package: Interface to Unidata netCDF Format Data Files
install.packages("ncdf4")                                             # new library to read a different format data files
library(ncdf4)                                                        
library(raster)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.NC")     # giving the name for visualization process, raster function to import a rasterlayer
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)        # we make a colorramppalette for snow changes blue to white

# Exercise: plot snow cover with the cl palette
plot(snowmay,col=cl)  

                                                                      # import snow data # for several layers all together in one time 
                                                                      # he aggregates and clumped several pixels to make the images less heavy
                                                                      # you can create a new folder i.e. snow # needed to stack multitemporal images in few command with lapply and stack functions
                                                                      # slow manner to import dataset
                                                                      # new working directory setwd("C:/lab/snow")
# setwd("~/lab/snow") #linux
# setwd("/Users/utente/lab/snow") #mac

setwd("C:/lab/snow") # windows                                        # we want to show the multitemporal scale (20years) of snow cover in a specific geografic area

snow2000 <- raster("snow2000r.tif")                                   # import all the temporal images raster layer of interest
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))                                                     # we want plot all images together we used a multiframe par function
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

                                                                      # import the images set in faster way
                                                                      # lapply function: you can repet the same fuction for the whole set (i.g. raster func.): Apply a Function over a List or Vector, lapply returns a list of the same length as X, each element of which is the result of applying FUN to the corresponding element of X
                                                                      # you can list files of the all folder. i.e. snow_2000, snow_2005, ..... snow_2020
                                                                      # create a list of files, list.file function: List the Files in a Directory/Folder. you must say the pattern at the console where it can keep the files
rlist <- list.files(pattern="snow")                                   # lapply function repet the raster fuction to import the interest layer of the all set snow images
rlist                                               
import <- lapply(rlist, raster)                                       # we import the layer (raster func.) for all the single files (images) with the function lapply
                                                                      # this work procedure is call stack or raster stack
                                                                      # stack function: stack of all the dataset that we import: Stack or Unstack Vectors from a Data Frame or List: stacking vectors concatenates multiple vectors into a single vector along with a factor indicating where each observation originated. Unstacking reverses this operation
snow.multitemp <- stack(import)                                       # we give the proper name at the vector (snow multitemp, because we are analysing snow cover in a multitemporal scale)
 
plot(snow.multitemp, col=cl)                                          # in ecology as well as in science is important to make several analysis all together at the same time with much powerful information range and in faster way.
                                                                      # less codes and commands, less time involved, same result
#########################################################################################################################################
# make a prediction of snow cover, making a graph with time(x) and snow cover %(y)
# with a regression functions we can fit our set to predict snow cover for the next 5 years using the trend of the last 20 years
# Big Data structure of complex data, we want see how make the importation of data.
# lets look at the function "prediction"
# let make a prediction of snow cover with a simple line code
                                                                      # prediction is a R code already prepared and added in the working directory folder (snow)
source("prediction.r")                                                # with source function you can select a file like a R script and it will work by themself for all my analysis!!!!
                                 
############################################################ what there are inside the prediction function that we created?
# require(raster)
# require(rgdal)
# define the extent
# ext <- c(-180, 180, -90, 90)
# extension <- crop(snow.multitemp, ext)   
# make a time variable (to be used in regression)
# time <- 1:nlayers(snow.multitemp)
# run the regression
# fun <- function(x) {if (is.na(x[1])){ NA } else {lm(x ~ time)$coefficients[2] }} 
# predicted.snow.2025 <- calc(extension, fun) # time consuming: make a pause!
# predicted.snow.2025.norm <- predicted.snow.2025*255/53.90828

##################### day second
                                                                       # Read R Code from a File, a Connection or Expressions
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


load("R_code_snow.r.RData")                                            # load preview work/script

prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)

# the idea is to make a regretion model to understand the future snow cover in continuity with the trend until today
# export the output
# you made the calculation and you want to send the output to a collegue, you can use writeRaster func.

writeRaster(prediction, "final.tif")  # is good func for transform images in tif, It create new data 1. we import data (5images), 2. we stack all that ims and added prediction, 3. with writeRaster
# writeRaster is the oppost of Raster function!! that import data into R. at the opp. writeRaster you can export the image from R to your PC folder.
# final stack (to see if the extent of the image is different)
# https://gdal.org/ to have a look in which format will work my function


final.stack <- stack(snow.multitemp, prediction)   # we are going to stack the all multitemp+ prediction of snow cover
plot(final.stack, col=cl)
# export the R graph for a thesis in pdf or png format

pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()

png("my_final_exciting_graph.png")
plot(final.stack, col=cl)
dev.off()

#############################################################################################################################
#############################################################################################################################

# 13. R_code_monitoring_air_pollution_no2.r

# R_code_monitoring_air_pollution_no2.r
# create a new folder to make lapply and stack dataset for this new set of images relating at the NO2 concentration evolution before and during the covid pandemic

setwd("C:/lab/NO2/")  

# R_code_no2.r
library(raster)

setwd("C:/lab/no2/")
# create RasterStack
rlist <- list.files(pattern="EN")                      # create a list of files. EN = environmental NOxs
import <- lapply(rlist, raster)                        # using lapply to apply the interested function in a multifiles list of images multitemporal scaled
EN <- stack(import)                                    # stack import dataset
cl <- colorRampPalette(c('red','orange','yellow'))(100)
plot(EN, col=cl)

par(mfrow=c(1,2))                                      # to see the change during Sar-CoV-2 pandemic time in UE. Two periods before and during pandemic
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)
                                                       # selection of 3 layers for RGB plot from EN images. See where and when Europe was more affected by NO2 pollution
plotRGB(EN, r=1, g=7, b=13, stretch="lin")
                                                       # difference map between 2 situations
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100)  
plot(dif, col=cld) 

                                                       # quantitative decrease of no2 
boxplot(EN)                                            # five-objects summary is the minimum, first quartile, median, third quartile, and maximum. The first quartile is the median of the data points to the left of the median.
boxplot(EN,outline=F)                                  # without outlines. making a box plot: boxplot function: produce box-and-whisker plot(s) of the given (grouped) values
boxplot(EN,outline=F, horizontal=T)                    # bars in horizontal position
boxplot(EN,outline=F, horizontal=T, axes=T)            # with axis informations
                                                       # multivariate analysis
plot(EN$EN_0001, EN$EN_0013)                           # plot EN stacked images list of NO2 pollution liked with the farest temporal images
abline(0,1,col="red")                                  # to see if the values are under the line, this means a decrease in NO2 concentration!! 45 degree line dividing x and y plan in equal parts 
                                                       # abline: fit the multivariate linear model: most of [NO2] values are under the line used to fit the dots cloud. means that in the first period(no lockdown) NO2 concentration values were higher

                                                       # fast version of import and plot of many data and stat. analysis for lazy people :)
rlist <- list.files(pattern="snow")
import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1) # most of the value under the curve
plot(snow.multitemp$snow2000r, snow.multitemp$snow2020r) 
abline(0,1,col="red")
                                                       # you can apply a code like this to NO2 analysis

#############################################################################################################################
#############################################################################################################################

# 14. R_code_crop_image.r

setwd("C:/lab/")
library(raster)
library(ncdf4)                                                        # required to read our data format ".nc"

snow <- raster( "c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)        # colorRamp with the meanest colours
plot(snow, col=cl)
                                     # making a crop of our image i.g. Italy 
ext <- c(0, 20, 35, 50)              # this is the extent(polygon) in which we want to zoom on
                                     # extent function: Objects of class Extent are used to define the spatial extent (extremes) of objects of the BasicRaster and Raster* classes. we have to say 1' the image vector and 2' the extent 
zoom(snow, ext=ext)                  # zoom function: Zoom in on a map (plot) by providing a new extent, by default this is done by clicking twice on the map
                                     # but we can also cutting the image with crop function 
                                     # in this case the extent is direct related at the previous one
                                     # this is very useful to crop the area of interest
snowitaly <- crop(snow, ext)         # crop function: crop returns a geographic subset of an object as specified by an Extent object (or object from which an extent object can be extracted/created). If x is a Raster* object, the Extent is aligned to x. Areas included in y but outside the extent of x are ignored (see extend if you want a larger area)
plot(snowitaly, col=cl)             
                                     # you can drow one specific area in other way. drowing the extent                       
                                     # to give a geometry for example a rectangle # () this command to say we are working in the image as a vector 
zoom(snow, ext=drawExtent())         # drawExtent: Create an Extent object by drawing on a map, click on two points of a plot (map) to obtain an object of class Extent ('bounding box'). draw and zoom can be done with zoom and drawExtent function

#############################################################################################################################
#############################################################################################################################

# 15. R_code_interpolation.r

# R_code_interpolation.r
# interpolation field data
# Interpolation of the data thanks spatstat library and species distribution modelling to understand the position and structure of the forest 
# Another nice analysis is about diameter and height of the forest
# Interpolation: spatstat library
# library(dbmss): Distance-Based Measures of Spatial Structures: simple computation of spatial statistic functions of distance to characterize the spa-tial structures of mapped objects, following Marcon, Trais-sac, Puech, and Lang (2015) <doi:10.18637/jss.v067.c03>.

setwd("C:/lab/")
library(spatstat)
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T)                 # import dataset/ data table of excel, ";" is the separator of dataset col. and header is true (first line is the header of the dataframe)
head(inp)                                                                   # look the head of dataset
attach(inp)                                                                 # attach the dataset inp
                                                                            # estimate canopy cover
plot(X,Y)                                                                   # XY geospatial coordinates
summary(inp)                                                                # let's see the minumum and maximum of X and Y, in order to give an extent to spatstat
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000))                 # range of min and max of X and Y, assign the coordinates to spatstat
                                                                            # give information about the variable: lables
names(inp)                                                                  # names of the variables: Canopy.cov (cover)
marks(inppp) <- Canopy.cov                                                  # marks function: Marks of a Point Pattern: extract or change the marks attached to a point pattern dataset.mark the variable with coordinates
                                                                            # Smooth function(spatstat): spatial smoothing of data: generic function to perform spatial smoothing of spatial data. Is an implementation of running median smoothers (algorithm proposed by Tukey) DOI: 10.17485/ijst/2016/v9i28/97354
canopy <- Smooth(inppp)                                                     # visualize the objects where they are not been measured in pixels
                                                                            # list validation distance of value from the line that record the values: means measured the amount of error
plot(canopy)                                                                # density of the vegetation
points(inppp, col="green")                                                  # adding points recorded on the map
 
                                                                            # adding the lichens for detecting the quality of air
marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)                                                      # same procedure used for canopy variable
plot(lichs)
points(inppp)                                                               # no congruence with canopy cover and lichens

par(mfrow=c(1,2))                                                           # par plot function to have a look of both canopy and lichens density distribution
plot(canopy)
points(inppp)                                                               # point refered at the inppp coordinates on canopy density graph
plot(lichs)
points(inppp)                                                              

par(mfrow=c(1,3))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)
plot(Canopy.cov, cop.lich.mean, col="red", pch=19, cex=2)                   # to see the negative correlation between canopy cover and lichens cover

# Dati psammofile from Giacomo                                              # interpolation analysis

inp.psam <- read.table("dati_psammofile.csv", sep=";", head=T)
attach(inp.psam)
head(inp.psam)
 
plot(E,N)                                                                  # clumped distribution North/East
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))           # giving coordinates
marks(inp.psam.ppp) <- C_org                                               # ecological data
C <- Smooth(inp.psam.ppp)                                                  # warning message:lower amount of data (or no point) for some part because of clumped set (Numerical underflow detected: sigma is probably too small)
plot(C)
points(inp.psam.ppp)                                                       # solution: mean value for each clumped zone or select the zone of the graph and zoom on the top of them: separation of the main graph in several graphs

#############################################################################################################################
#############################################################################################################################

# 16. R_code_sdm.r

# Species Distrubution Modelling

install.packages("sdm")     # install.packages sdm: Species Distribution Modelling. Description: An extensible framework for developing species distribution models using individual and community-based approaches, generate ensembles of models, evaluate the models, and predict species potential distributions in space and time. For more information, please check the following paper: Naimi, B., Araujo, M.B. (2016) <doi:10.1111/ecog.01881>. our packages very useful for ecology and where we can find our dataset of species distribution
# install.packages("rgdal")
library(sdm)
library(raster)             # library to read the raster files and to make use the raster dataset # raster packages are dependent from sp packages, you can solve the problem with " library( rgdal, "dependencies=T")
library(rgdal)              # geodata astract library # use to import the packs of the species data spatial distrubution for this project

# species :  we are looking for a species present-absence model (Brachypodium rupestre)
file <- system.file("external/species.shp", package="sdm")   # the location of the species data, automatically importation of the external folder, insiede the folder external there is the species
                                                             # system.file function: Find Names of R System Files. Finds the full file names of files in packages etc. Needed to import the file from the sdm package
species <- shapefile(file)                                   # shapefile raster function: read or write a shapefile: we can use the graphical part of the file using shapefile func., this type of files are very used today for graphical and mapping dataset
                                                             
species                                                      # looking at the estent and features that are the points in this case
species$Occurrence                                           # let's see the occurence of species, present-absent model (200 features)
plot(species)

plot(species[species$Occurrence == 1,],col='blue',pch=16)    # making a condition inside the "[]" and to make it we need "==" command, we decided to show only the presence occurences = 1
                                                             # point function: plots points on a map or plot previuosly processed
points(species[species$Occurrence == 0,],col='red',pch=16)   # making the oppost condition= absence with red points (blue points are blue)

# path to the folder contains the data
path <- system.file("external", package="sdm")               # making the path inside the external using the same function of before system.file. the folder is called external 

lst <- list.files(path=path,pattern='asc$',full.names = T)   # creating a list of objects. the extention of the files needed in this case is "asc" (other extention may be tif or png, ect), dollar can be avoided, but to be sure to carry out all data
lst                                                          # have a look of files into the list, into smd pkgs there is external folder, into external there are species distribution + ecological and environmental variables

preds <- stack(lst)                                          # making a raster objec (stack function). We are making a stack of dataset, elevation, temperature, precipitation, vegetation

cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)                                          # we can analyze the predictor informations  graphically, we found that many occorences are at the bottom of the valley, in case of vegetation presence, because Brachypodium rupestre love small elevation and grasslands
#par(mfrow=c(2,2))                                           # ecologicic factors influencing the species distribution, ecological analysis
plot(preds$elevation, col=cl)                                
points(species[species$Occurrence == 1,], pch=16)            # plot species occurrence with elevetion map

plot(preds$temperature, col=cl)                              # plot species occurrence with temperature map
points(species[species$Occurrence == 1,], pch=16)            # Brachypodium rupestre love high temperature and low elevation

plot(preds$precipitation, col=cl)                            # Brachypodium rupestre love medium to high precipitation, in the case medium precipitation
points(species[species$Occurrence == 1,], pch=16)            # plot species occurrence with precipitation map
                     
plot(preds$vegetation, col=cl)                               # vegetation index is gave from the NDVI theory, Brachipodium rupestre love high ammount vegetation or not? This analysis is called exploration analysis. Brachipodium rupestre love medium vegetation
points(species[species$Occurrence == 1,], pch=16)            # plot species occurrence with vegetation map
                                                             # making the model thanks all this information. One is train = species the other argument of the model is predictors= preds.
 
d <- sdmData(train=species, predictors=preds)                # sdmData function: creating sdm Data object: Creates a sdmdata objects that holds species (single or multiple) and explanatory variates. In addition, more information such as spatial coordinates, time, grouping variables, and metadata (e.g., author, date, reference, etc.) can be included. class = sdmData, n species= 1, number or record 200, type of data= presence-absent, plot, abundance of species but we are looking only prec. absen. is good for the velocity of the sampling efforce.
d                                                            # there is a negative correlation of linear model between species occ. and elevation. In general we use a logistic model that reppresent better ecological functions
                                                             # sdm function: Fit and evaluate species distribution modelsmodel: Fits sdm for single or multiple species using single or multiple methods specified by a user in methods argument, and evaluates their performance
m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods='glm') # occurence, ~ symble is like an =, y~a+bx in this case x is negative related with y in case of elevation. model: y= a+ b elev. + c + temperature + d precipitation + e vegetation, for all the variable we are using a different curve.
                                                             # finally the methods of our model (glm) : generalised linear model is the model that we are going to use.
p1 <- predict(m1, newdata=preds)                             # final prediction, related at our predictors used
                                                             # predict function: model predictions: predict is a generic function for predictions from the results of various model fitting functions. The function invokes particular methods which depend on the class of the first argument
plot(p1, col=cl)                                             # plotting the prediction model
points(species[species$Occurrence == 1,], pch=16)  

s1 <- stack(preds, p1)                                       # stacking all the variables maps + prediction model of species presence 

#############################################################################################################################
#############################################################################################################################




















# 17. R_code_myproject_exam.r

setwd("C:/upsp/")

library(raster)
library(rgdal)
library(RStoolbox)
library(rasterdiv)
library(rasterVis)
library(ggplot2)
library(ncdf4)                                                        # required to read our data format ".nc"
library(spatstat)
library (sf)  

setwd("C:/upsp/DMP")

# stack images dmp
rlistdmp <- list.files(pattern="DMP")
import <- lapply(rlistdmp, brick)
dmp.multitemp <- stack(import)
cl <- colorRampPalette(c('yellow','light green','dark green'))(100)
plot(dmp.multitemp, col=cl)    

#dmpitaly <- crop(dmp.multitemp, ext)         # crop function: crop returns a geographic subset of an object as specified by an Extent object (or object from which an extent object can be extracted/created). If x is a Raster* object, the Extent is aligned to x. Areas included in y but outside the extent of x are ignored (see extend if you want a larger area)
#plot(dmpitaly, col=cl) 

#level plot dmp
#par(mfrow=c(2,1))
dmp2014 <- brick("c_gls_DMP300-RT5_QL_201406100000_GLOBE_PROBAV_V1.0.1.TIFF")
levelplot(dmp2014)
dmp2018 <- brick("c_gls_DMP300-RT5_QL_201806100000_GLOBE_PROBAV_V1.0.1.TIFF")
levelplot(dmp2018)
# just to have a look how much the dry matter production change in north emisphere in two weeks in summer time
dmp2019 <- brick("c_gls_DMP300-RT5_QL_201905310000_GLOBE_PROBAV_V1.0.1.TIFF")
levelplot(dmp2019)


### multivariate analysis of PCA
plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1) # most of the value under the curve
plot(snow.multitemp$snow2000r, snow.multitemp$snow2020r) 
abline(0,1,col="red")

#



#drop images
dmp2019 <- brick("c_gls_DMP300-RT5_QL_201905310000_GLOBE_PROBAV_V1.0.1.TIFF")
ext <- c(0, 20, 35, 50)
zoom(dmp2019, ext=ext)
zoom(dmp2019, ext=drawExtent()) 

# drop and stack images
rlistdmp <- list.files(pattern="DMP")
import <- lapply(rlistdmp, brick)
dmp.multitemp <- stack(import)
#cl <- colorRampPalette(c('yellow','light green','dark green'))(100)
#plot(dmp.multitemp, col=cl) 
levelplot(dmp.multitemp)

ext <- c(0, 20, 35, 50)
zoom(import, ext=ext)
zoom(dmp2019, ext=drawExtent())




















# 14. R_code_crop_image.r

setwd("C:/lab/")
library(raster)
library(ncdf4)                                                        # required to read our data format ".nc"

snow <- raster( "c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)        # colorRamp with the meanest colours
plot(snow, col=cl)
                                     # making a crop of our image i.g. Italy 
ext <- c(0, 20, 35, 50)              # this is the extent(polygon) in which we want to zoom on
                                     # extent function: Objects of class Extent are used to define the spatial extent (extremes) of objects of the BasicRaster and Raster* classes. we have to say 1' the image vector and 2' the extent 
zoom(snow, ext=ext)                  # zoom function: Zoom in on a map (plot) by providing a new extent, by default this is done by clicking twice on the map
                                     # but we can also cutting the image with crop function 
                                     # in this case the extent is direct related at the previous one
                                     # this is very useful to crop the area of interest
snowitaly <- crop(snow, ext)         # crop function: crop returns a geographic subset of an object as specified by an Extent object (or object from which an extent object can be extracted/created). If x is a Raster* object, the Extent is aligned to x. Areas included in y but outside the extent of x are ignored (see extend if you want a larger area)
plot(snowitaly, col=cl)             
                                     # you can drow one specific area in other way. drowing the extent                       
                                     # to give a geometry for example a rectangle # () this command to say we are working in the image as a vector 
zoom(snow, ext=drawExtent())         # drawExtent: Create an Extent object by drawing on a map, click on two points of a plot (map) to obtain an object of class Extent ('bounding box'). draw and zoom can be done with zoom and drawExtent function

# 12. R_code_snow.r

# R_code_snow.r 

# setwd("~/lab/") #linux
# setwd("/Users/utente/lab") #mac
setwd("C:/lab/") 
                                                                      # ncdf4 packages: returns a string that is the version number of the ncdf4 package: Interface to Unidata netCDF Format Data Files
install.packages("ncdf4")                                             # new library to read a different format data files
library(ncdf4)                                                        
library(raster)

snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.NC")     # giving the name for visualization process, raster function to import a rasterlayer
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)        # we make a colorramppalette for snow changes blue to white

# Exercise: plot snow cover with the cl palette
plot(snowmay,col=cl)  

                                                                      # import snow data # for several layers all together in one time 
                                                                      # he aggregates and clumped several pixels to make the images less heavy
                                                                      # you can create a new folder i.e. snow # needed to stack multitemporal images in few command with lapply and stack functions
                                                                      # slow manner to import dataset
                                                                      # new working directory setwd("C:/lab/snow")
# setwd("~/lab/snow") #linux
# setwd("/Users/utente/lab/snow") #mac

setwd("C:/lab/snow") # windows                                        # we want to show the multitemporal scale (20years) of snow cover in a specific geografic area

snow2000 <- raster("snow2000r.tif")                                   # import all the temporal images raster layer of interest
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))                                                     # we want plot all images together we used a multiframe par function
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

                                                                      # import the images set in faster way
                                                                      # lapply function: you can repet the same fuction for the whole set (i.g. raster func.): Apply a Function over a List or Vector, lapply returns a list of the same length as X, each element of which is the result of applying FUN to the corresponding element of X
                                                                      # you can list files of the all folder. i.e. snow_2000, snow_2005, ..... snow_2020
                                                                      # create a list of files, list.file function: List the Files in a Directory/Folder. you must say the pattern at the console where it can keep the files
rlist <- list.files(pattern="snow")                                   # lapply function repet the raster fuction to import the interest layer of the all set snow images
rlist                                               
import <- lapply(rlist, raster)                                       # we import the layer (raster func.) for all the single files (images) with the function lapply
                                                                      # this work procedure is call stack or raster stack
                                                                      # stack function: stack of all the dataset that we import: Stack or Unstack Vectors from a Data Frame or List: stacking vectors concatenates multiple vectors into a single vector along with a factor indicating where each observation originated. Unstacking reverses this operation
snow.multitemp <- stack(import)                                       # we give the proper name at the vector (snow multitemp, because we are analysing snow cover in a multitemporal scale)
 
plot(snow.multitemp, col=cl)                                          # in ecology as well as in science is important to make several analysis all together at the same time with much powerful information range and in faster way.
                                                                      # less codes and commands, less time involved, same result
#########################################################################################################################################
# make a prediction of snow cover, making a graph with time(x) and snow cover %(y)
# with a regression functions we can fit our set to predict snow cover for the next 5 years using the trend of the last 20 years
# Big Data structure of complex data, we want see how make the importation of data.
# lets look at the function "prediction"
# let make a prediction of snow cover with a simple line code
                                                                      # prediction is a R code already prepared and added in the working directory folder (snow)
source("prediction.r")                                                # with source function you can select a file like a R script and it will work by themself for all my analysis!!!!
                                 
############################################################ what there are inside the prediction function that we created?
# require(raster)
# require(rgdal)
# define the extent
# ext <- c(-180, 180, -90, 90)
# extension <- crop(snow.multitemp, ext)   
# make a time variable (to be used in regression)
# time <- 1:nlayers(snow.multitemp)
# run the regression
# fun <- function(x) {if (is.na(x[1])){ NA } else {lm(x ~ time)$coefficients[2] }} 
# predicted.snow.2025 <- calc(extension, fun) # time consuming: make a pause!
# predicted.snow.2025.norm <- predicted.snow.2025*255/53.90828

##################### day second
                                                                       # Read R Code from a File, a Connection or Expressions
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


load("R_code_snow.r.RData")                                            # load preview work/script

prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)

# the idea is to make a regretion model to understand the future snow cover in continuity with the trend until today
# export the output
# you made the calculation and you want to send the output to a collegue, you can use writeRaster func.

writeRaster(prediction, "final.tif")  # is good func for transform images in tif, It create new data 1. we import data (5images), 2. we stack all that ims and added prediction, 3. with writeRaster
# writeRaster is the oppost of Raster function!! that import data into R. at the opp. writeRaster you can export the image from R to your PC folder.
# final stack (to see if the extent of the image is different)
# https://gdal.org/ to have a look in which format will work my function


final.stack <- stack(snow.multitemp, prediction)   # we are going to stack the all multitemp+ prediction of snow cover
plot(final.stack, col=cl)
# export the R graph for a thesis in pdf or png format

pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()

png("my_final_exciting_graph.png")
plot(final.stack, col=cl)
dev.off()
#############################################################################################################################
#############################################################################################################################











