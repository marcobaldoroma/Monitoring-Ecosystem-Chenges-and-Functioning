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
