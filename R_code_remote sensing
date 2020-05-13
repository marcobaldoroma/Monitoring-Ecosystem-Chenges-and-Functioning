# R code for  Remote Sensing 

setwd ("C:/lab/")

library (raster)

install.packages("RStoolbox")                                            #packages to make our analysis, the book that the Prof. wrote is about the algorith of this R function create by him


library (RStoolbox)                                                      #to be faster in the istallation of pkgs: install.packages(c("raster", "RStoolbox")

load("C:\\lab\\biomes_types.csv")                                        #load biomes_types and biomes data set
load("C:\\lab\\biomes.csv")

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

                                                                          #change the lay out of our graph in 4 rows and 1 collum and to make it different and stretched, could be really usefull for data analysis

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

dev.off ()                                                    # close others plots

                                                              # plotRGB function RGB colors

plotRGB (p224r63_2011, r=3, g=2, b=1, stretch="Lin")
                                                              # stretching as much possible the colors

                                                              # substitute red component with our NIR # in this way you can discriminete better the vegetation in the images, dark green is the forest, red and pink in general the pitch 
plotRGB (p224r63_2011, r=4, g=3, b=2, stretch="Lin")

                                                              # EX. NIR on top of the G component # NIR false colours
plotRGB (p224r63_2011, r=3, g=4, b=2, stretch="Lin")
                                                              # dark forest is forest with high ammount of water

                                                              # yellow area are bare in this way to set the colors
plotRGB (p224r63_2011, r=3, g=2, b=4, stretch="Lin")

################ day 2

#setwd("~/lab/") # linux
setwd("C:/lab/") # windows
# setwd("/Users/nome/Desktop/lab") # mac

load("remote_sensing.RData")

                                                               # list
ls()

p224r63_1988 <- brick("p224r63_1988_masked.grd")

plot(p224r63_1988)

p224r63_2011 <- brick("p224r63_2011_masked.grd")

ls()

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

                                                                           #giving colors at the div
cl <- colorRampPalette(c('"darkorchid3"','light blue','lightpink4'))(100) 
plot(dvi2011)

#or

cldvi <- colorRampPalette(c('light blue','light green','green'))(100) # 
plot(dvi2011, col=cldvi)

#or

cl2 <- colorRampPalette(c('"yellow"','light blue','lightpink4'))(100)
plot(dvi2011,col=cl2)

                                                                           # Exercise: dvi for 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre

cl3 <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100)
plot(dvi1988,col=cl3)

                                                                            #differances between 1988-2011 
diff <- dvi2011 - dvi1988
plot diff

                                                                            #aggregate 
                                                                            #changing the grain and resampling with factor=10
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)
 
                                                                            #EX. plot all together images with diff grain
par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin") 

                                                                             #information about image p224r63_2011
p224r63_201


