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

pdf("covid_density2.pdf")

                                                                        # higher number of intermediate colours (1000) to have more colours quality. The best option is a balance between the heavyness of the file and the quality of the colours
clrr <-colorRampPalette(c("light green", "yellow","orange","violet")) (1000)
plot(d, col=clrr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

pdf("covid_density3.pdf")                                               # this command is a wrapper for the pdf function

png("covid_density.png")                                                
                                                                        # or export in png data format. png function: graphics devices for BMP, JPEG, PNG and TIFF format bitmap files

dev.off()                                                               # dev.off function: Control Multiple Devices: these functions provide control over multiple graphics devices, with the command off, we say at the console to close all the plots in my environment

clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

dev.off()                                                                # to close all open plots

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
copNDVI10 <- aggregate(copNDVI, fact=10)                                               # levelplot function: Level plots and contour plots. Draws false color level plots and contour plots
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
                                                                       # par multiframe of 1 row and 4 columns, where i can plot all the 4 examples made in the R_code_reflectance                 
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









