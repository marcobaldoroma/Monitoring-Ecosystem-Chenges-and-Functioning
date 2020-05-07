# R_code_PCA_remote_sensing

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
