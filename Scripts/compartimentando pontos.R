# compartimento de pontos

# nao funciona
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

raster('Bioclims/bio1.asc') -> bio1

raster_sampler(bio1,1,10000) -> bg
# nao funciona

#bg

install.packages('dismo')

bg <- dismo::randomPoints(envs.bg[[9]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(occs)

get.checkerboard1(amb.10.asc, envs, aggregation.factor = 5, gridSampleN = 10000) -> checkerboard1.amb

