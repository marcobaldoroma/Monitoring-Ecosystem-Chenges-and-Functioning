# R_code_radiance

library(raster)

# create a new raster dataset 

toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)
values(toy) <- c(1.13,1.44,1.55,3.4)

plot(toy)
text(toy, digits=2)      # put the information on the pixel

toy2bits <- stretch(toy,minv=0,maxv=3)    # we traslate the information = radiance in bits information

storage.mode(toy2bits[]) = "integer"      # the way to starage our information 

plot(toy2bits)
text(toy2bits, digits=2)

toy4bits <- stretch(toy,minv=0,maxv=15)     # we can change the information from 2bits to 4bits
 

storage.mode(toy4bits[]) = "integer"
 
plot(toy4bits)
text(toy4bits, digits=2)
                                # passing from 4 bits to 8
                                                       
                                
toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer"

 

plot(toy8bits)
text(toy8bits, digits=2)                              
                                
                                
par(mfrow=c(1,4))                           
                                
                                
