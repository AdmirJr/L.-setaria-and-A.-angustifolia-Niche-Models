#crop rasters with Atlantic Forest shape and fix mismatches

library(raster)
library(rgdal)
library(maptools)

setwd("C:/Users/Victor/Documents/PPG-BOT/Tainara/planilhas/limpas/Atlantic_Forest")
x <- read.csv("C:/Users/Victor/Documents/PPG-BOT/Tainara/planilhas/limpas/sp_menos4.csv", sep=";")
#x<-read.csv("F:/Niche_Models/species.csv", stringsAsFactors=FALSE)
#x<-x[-68,]
str(x)
names(x)
splist = unique(x[,1])
splist

raster0=raster("C:/Users/Victor/Documents/PPG-BOT/Tainara/AllZero_NT.asc")
#raster_m=raster("F:/Niche_Models/to_merge.asc")

ext1<-"_atl"	
Shp1 = readOGR("C:/Users/Victor/Documents/Shapes/Mata_atl_toda_oficial_merge.shp")
#i=1
  for (i in 1:length(splist))
  {
    path<-paste0("C:/Users/Victor/Documents/PPG-BOT/Tainara/planilhas/limpas/",
                       splist[i], "_ECO.asc")

    r1 = raster(path)
    cr1 = crop(r1,Shp1)
    cr2 = mask(cr1,Shp1)
    #cr3 = merge(cr2,raster0)
    #cr4 = cr3+raster0
    #cr5 = merge(cr4,raster_m)
    ### For file with extension.
    #FileName = paste(substr(filelist[i],1,nchar(filelist[i])-4),substr(filelist[i],nchar(filelist[i])-3,nchar(filelist[i])),ext1,sep="")
  #  species <- as.character(splist[i])
   # FileName = paste(species,ext1,sep="")
    ### For files without extension
    writeRaster(cr2,FileName, "ascii") # filename = as.character()
    plot(cr4, main=splist[i])
    print(splist[i])
  }
rm(list = ls()); gc() #Liberar memoria RAM
