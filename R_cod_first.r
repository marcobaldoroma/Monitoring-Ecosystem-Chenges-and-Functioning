install.packages("sp")

library(sp)
data(meuse)

# Let's see how meuse dataset is structure:
meuse

# let's look at the firsts row of the set
head(meuse)

#let's lot two variables
#let's see if the zin concentration is realate to that  copper
attach(meuse)
plot(zinc,copper)
plot(zinc,copper,col="green")
#riempimento
plot(zinc,copper,col="green",pch=19)
#big plot mod, cex caracter exageration= 2, if you put 0-1 to have less size
plot(zinc,copper,col="green",pch=19,cex=2)

