
# Raster sampler -- write xy randomly selected

# This function writes xy coordinates taken at random from a particular raster file using it as a mask or extent.
# Author: Chris Hensz, Date: May 1, 2013

# First read a raster layer (raster object) that will be the reference to extract coordinates from

raster_sampler=function(rasterlayer,numberpoints,replications){
output=list()
for(i in 1:replications){
	output[[i]]=
		cbind(
			run=rep(i,numberpoints),
			sampleRandom(rasterlayer, numberpoints, na.rm=TRUE, xy=TRUE)[,1:2]
			)
	}
output2=NULL
for(i in 1:replications){output2=rbind(output2,output[[i]])}
write.csv(output2,file=paste(names(rasterlayer),"_rand",".csv",sep=""))
}

# raster_sampler(rasterobject,100,100)
