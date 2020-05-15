#point_pattern_analysis: Density map

install.packages("spatstat")
library(spatstat)

setwd("C:/lab/")


covid <- read.table("covid_agg.csv", head=T)                           # inport and read covid data base
 
attach(covid)
head(covid)

                                                                       # set the coordinates of the vectors in covids in relation to the global map
covids <- ppp(lon, lat, c(-180,180), c(-90,90))                        # c it used to clastering all the variables/numeber together. i.g. dead <- 12,34,55,66,77,88,89 for the lat and lon we used min and max

                                                                       # attach covid at the coordinate : covids <- ppp(covid$lon, covid$lat, c(-180,180), c(-90,90))
d <- density(covids)                                                   # create a variable names d = density

plot (d)                                                               # plot this new variable

points(covids)                                                         # add points at the density plot

                                                                       #################################################### second part

setwd("C:/lab/")                                                       # set work directory

load("point_pattern_analysis.RData")                                   ## load previus work

                                                                       #to use vector format in coastline

library(spatstat)                                                      # call libraries requested
library(rgdal)

ls()                                                                   #list of objects I have 

                                                                       #plot density + points in covids
plot(d)
points(covids)

                                                                       # import in lab and then in R ne_10m_coastline.shp"
                                                                       # to use vector format in coastline
                                                                       # let’s input vector lines (x0y0, x1y1, x2y2..)
coastlines <- readOGR("ne_10m_coastline.shp")
                                                                       #install additional packages, alternative way to plot the coastline at the density plot
                                                                       #install.packages("rnaturalearth")
                                                                       #coastlines <- rnaturalearth::ne_download(scale = 10, type = 'coastline', category = 'physical')
plot(d)
points(covids)

                                                                       # plot density + coastline of the world for covid-19 num. of cases
plot(coastlines, add=T)

cl <-colorRampPalette(c("yellow","orange","red")) (100)                # change the colour and make the graph beautiful # build a new object names "cl"

                                                                       
plot(d, col=cl, main="Densities of covid-19")                          # Plot the new object with density and name of the label
points(covids)
plot(coastlines, add=T)

                                                                       # Export your Densities of covid-19 plot in Pdf
pdf("covid_density.pdf")

                                                                       ## number of colours: abrupt change of colours!!! making just an example
cll <- colorRampPalette(c("light green", "yellow","orange","violet")) (5)
plot(d, col=cll, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

                                                                       # Exercise: new colour ramp palette
clr <-colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)

pdf("covid_density2.pdf")

                                                                        # higher number of intermediate colours (1000)
clrr <-colorRampPalette(c("light green", "yellow","orange","violet")) (1000)
plot(d, col=clrr, main="Densities of covid-19")
points(covids)

plot(coastlines, add=T)
pdf("covid_density3.pdf")

png("covid_density.png")                                                # This command is a wrapper for the pdf function.
                                                                        # or export in png

dev.off()                                                               # to close all plots

clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")

points(covids)
plot(coastlines, add=T)

dev.off()
