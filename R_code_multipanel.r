### Multipanel in R: seeing correlation amoung ecological variable

                                                 # install.packages("sp")
install.packages("GGally")                       # this is required to use the function ggpairs() we would graphy

                                                 #require those library
library(sp)

data(meuse)                                      # there is a dataset available named meuse

                                                 #to fix the data set
attach(meuse)

                                                 # Exercise: see the names of the variables and plot cadmium versus zinc
                                                 # There are two ways to se the names of the variables:
names(meuse)
head(meuse)

                                                 #plot my variables
plot(cadmium,zinc,pch=15,col="red",cex=2)

                                                 # Exercise: make all the possible pairs plots of the dataset
plot(x,cadmium)
plot(x,zinc)
plot(x,lead)
plot(x,copper)
                                                 # plot(cadmium,zinc,pch=15,col="red",cex=2) an amateur exercise to prettify plots
                                                 # plot is not a good idea sometimes
pairs(meuse)
                                                 # selection of variables we want to plot together
pairs(~ cadmium + copper + lead + zinc, data=meuse)

pairs(meuse[,3:6])

                                                 # Exercise: prettify this graph
pairs(meuse[,3:6], col="red", pch=18, cex=1.5)


                                                 
pairs(meuse[,3:6])

                                                 # from 3(cadmium) to 6(zinc)
                                                 # 16 and 15 are the same function!

pairs(meuse[,3:6],pch=19)

library(GGally)                                  # GGally package will show infographic even better
ggpairs(meuse[,3:6])
