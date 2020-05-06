# R code for multivariate analysis

setwd ("C:/lab/")

#packages for vegetation analysis already installed. call them
library (vegan)

#give the name and the linkage, give the command to read exel or grid dataset files, give the header name "first line" of the collums "intestazione", in the end the separator of the collums = comma
#biomes             <-          read.table     biomes.csv csv #tipe of the data file   header = T #first line of my table                                  sep = ","          

biomes <- read.table ("biomes.csv", header=TRUE, sep = ",")

#have a look at the dataset
head (biomes) # view (biomes), biomes

#for multivariate analysis in many programs is bit difficult trick built it, in R scientist already create a function to facilities it, names DECORANA
#DEtreanded COrrespondence ANAlysis = DECORANA
multivar <- decorana (biomes)

plot (multivar)

#put an exageration
plot (multivar, cex=1.2)

plot (multivar)

# we can see in the R macchine the percentage of the reality in 2 dimentions. our algorithm give us the 82% of our data just in 2 dimantion, this is not bad 
multivar    
# (52 + 30) % of real dataset of 20 dimentions in 2, x,y . we can see it trought calling our new variable #multivar, and having a look at Eigenvalues = DCA1 0.5117 #the x axis and  DCA2 0.3036 # the y axis, our second axis   

#let liked every same biomes per plot
#biomes_types, codify for R the new data set

biomes_types <- read.table ("biomes_types.csv", header=TRUE, sep = ",")
head ("biomes_types")

#now we are link every point = plot, one each other of the same biomes type
#first attach the dataset 

attach (biomes_types)

#order inside an ellipse our biomes/plot, the variables, collum types, color,  form of the shape,  size of the line


ordiellipse (multivar, type, col = 1:4, kind = "ehull", lwd = 3)

# we use it to connect the species at the ellipse, this can be a way to add a dimention at our graph
ordispider(multivar, type, col=1:4, label = T)

#dev_off ()
#save.image("C:\\lab\\R_code_multivar")
