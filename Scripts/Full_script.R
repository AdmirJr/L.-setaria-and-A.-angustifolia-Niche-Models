# Downloading Leptasthenura setaria occurence data -----------------------------
# Data bases: GBIF, iDigBio, SpeciesLink, VertNet and iNat
# Packages

library(readxl)

# Downloading the data

"Leptasthenura setaria" -> sp # to make it easier

## Creating a function to download the occ data

download.occ <- function(species) {
  
  require(tidyverse)
  require(rgbif)
  require(ridigbio)
  require(spocc)
  
  
  dt <- data.frame()
  
  # GBIF
  gbif.data <- occ_data(scientificName = species,
                        hasCoordinate = TRUE,
                        limit = 100000)
  
  gbif.data <- select(gbif.data$data, 
                      "scientificName", 
                      "decimalLatitude", 
                      "decimalLongitude") 
  
  colnames(gbif.data) <- c("name", "latitude", "longitude")
  gbif.data[,'prov'] = 'gbif'
  
  
  rbind(dt, gbif.data) -> dt
  rm(gbif.data)
  
  # iDigBio 
  idigbio.data <- idig_search_records(rq = list(scientificname = species, 
                                                geopoint = list(type = "exists"))) %>%
    select(scientificname, geopoint.lat, geopoint.lon)
  
  colnames(idigbio.data) <- c("name", "latitude", "longitude")
  idigbio.data[,'prov'] = 'idigbio'
  
  rbind(dt, idigbio.data) -> dt
  rm(idigbio.data)
  
  # VertNet e iNat 
  databases <- c('inat', 'vertnet')
  
  spocc.data <- occ(query = species,
                    from = databases,
                    limit = 100000,
                    has_coords = TRUE) %>% 
    occ2df() %>% 
    select(name, latitude, longitude, prov)
  
  rbind(dt, spocc.data) -> dt
  rm(spocc.data)
  
  return(dt)
  
}

## Downloading the data

dados.brutos <- download.occ(sp)

## loading SpeciesLink data (downloaded manualy)

specieslink.data <- read_excel(path = 'SpeciesLink/L_set_occ.xlsx')

specieslink.data <- specieslink.data[, c("scientificname", "latitude", "longitude")]
specieslink.data[, "prov"] = "specieslink"
colnames(specieslink.data) <- c("name", "latitude", "longitude", "prov")

rbind(dados.brutos, specieslink.data) -> dados.brutos
rm(specieslink.data)

# Adding a ID column

dados.brutos <- tibble::rowid_to_column(dados.brutos, "ID")

# Salvando the dataset

write.csv(dados.brutos, file = 'Datasets/Uncleaned data.csv')

# Cleaning Leptasthenura setaria occurence data --------------------------------
# Packages

library(maps)
library(CoordinateCleaner)

# Removing NAs

dados.brutos <- dados.brutos %>% 
  filter(!is.na(longitude)) %>% 
  filter(!is.na(latitude))

# Plotting the points on a map

plot(latitude ~ longitude, dados.brutos)
map(add = T)

# Cleaning the coordinates

## Converting 'latitude' and 'longitude' columns from 'chr' to 'dbl'

dados.brutos[,c('latitude', 'longitude')] <- sapply(dados.brutos[,c('latitude','longitude')], as.numeric)

## Then actually cleaning the coordinates

dados.coord <- clean_coordinates(dados.brutos,
                                 lon = 'longitude',
                                 lat = 'latitude',
                                 species = 'name',
                                 value = 'clean')

## Plotting the points on a map again

plot(latitude ~ longitude, dados.coord)
map(add = T)

## Manualy removing a weird an unlikely point

dados.coord <- subset(dados.coord, (latitude < -20))

## Plotting the points on a map again

plot(latitude ~ longitude, dados.coord)
map(add = T) ### it's ready

# "Lepthastenuria setaria" is the currently accepted name for the species 
# so let's rename the column "name"

dados.coord[, 'name'] = sp

# Removing duplicates

dados.coord <- dados.coord %>% 
  distinct(latitude, longitude, .keep_all = TRUE)

# Plotting the map once more

plot(latitude ~ longitude, dados.coord)
map(add = T)

# Saving the data

write.csv(dados.coord, file = 'Datasets/Cleaned coordinates.csv')

# Spatial filtering ------------------------------------------------------------
## Packages

library(spThin)
library(maps)

# Filtering

reps <- c(10, 50, 100)
distances <- c(5, 10, 15, 20, 25, 30)

for (k in distances) {
  
  for (i in reps) {
    
    message('Starting for ', k,' kms and ', i,' reps')
    
    thin(loc.data = dados.coord,
         lat.col = 'latitude',
         long.col = 'longitude',
         spec.col = 'name',
         thin.par = k,
         reps = i,
         locs.thinned.list.return = FALSE,
         write.files = TRUE,
         max.files = 1,
         out.dir = 'Spatial_thinned_data/',
         out.base = k,
         write.log.file = FALSE)
    
    save.image('.RData')
    
  }
  
}


# Loading the several .csv at once

temp <- list.files("Spatial_thinned_data", 
                   pattern = ".csv",
                   full.names = TRUE)

spatial.thinned <- lapply(temp, read.csv)

# Saving the plots for later consult

pdf('maps/map.pdf')

for (p in c(1:length(spatial.thinned))) {
  
  plot(latitude ~ longitude, spatial.thinned[[p]])
  map(add = T)
  title(main = p)
  
}

dev.off()

write.csv(spatial.thinned[[6]], file = 'Datasets/spat_filt_10km.csv')

# Environmental filtering ------------------------------------------------------
## Packages
library(flexsdm)
library(raster)
library(terra)

# Downloading environmental data

climate <- getData("worldclim", var="bio",res=2.5)

# Selecting only the South America

plot(climate[[1]])
drawExtent() -> e ### Select the SA manualy 
climate.focus <- crop(climate, e)
stack(climate.focus) -> climate.focus
plot(climate.focus)
writeRaster(climate.focus, 'SA_climate/SA_climate')

climate.focus <- rast(climate.focus)

# Actually filtering the data

amb.10 <- occfilt_env(data = dados.coord,
                         x = 'longitude',
                         y = 'latitude',
                         id = 'ID',
                         env_layer = climate.focus,
                         nbins = 10)

amb.20 <- occfilt_env(data = dados.coord,
                         x = 'longitude',
                         y = 'latitude',
                         id = 'ID',
                         env_layer = climate.focus,
                         nbins = 20)

amb.30 <- occfilt_env(data = dados.coord,
                      x = 'longitude',
                      y = 'latitude',
                      id = 'ID',
                      env_layer = climate.focus,
                      nbins = 30)

amb.40 <- occfilt_env(data = dados.coord,
                      x = 'longitude',
                      y = 'latitude',
                      id = 'ID',
                      env_layer = climate.focus,
                      nbins = 40)

# Calibrating M ----------------------------------------------------------------

# Packages

library(spocc)
library(rgdal)
library(rgeos)
library(ellipsenm)

# Downloading ecorregions

download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
              destfile = file.path(getwd(), "wwf_ecoregions.zip"))

unzip(file.path(getwd(), "wwf_ecoregions.zip"), 
      exdir = file.path(getwd(), "WWF_ecoregions"))

file.remove(file.path(getwd(), "wwf_ecoregions.zip"))

ecor <- readOGR("WWF_ecoregions/official", layer = "wwf_terr_ecos")

# Drawing areas

## Buffering points
M_buffer <- buffer_area(dados.coord, longitude = "longitude", latitude = "latitude", 
                          buffer_distance = 100)

## Using convex hulls
M_convex <- convex_area(dados.coord, longitude = "longitude", latitude = "latitude", 
                          buffer_distance = 75)

## Selecting polygons
M_ecorreg <- polygon_selection(dados.coord, longitude = "longitude", latitude = "latitude",
                                 polygons = ecor, buffer_distance = 25)


M_intersect <- gIntersection(M_buffer, M_convex)
M_intersect <- gIntersection(M_intersect, M_ecorreg)

writeOGR(M_intersect, ".", "M_intersection", driver = "ESRI Shapefile") #error

### trying to solve the error

p.df <- data.frame( ID=1:length(M_intersect))
rownames(p.df)

pid <- sapply(slot(M_intersect, "polygons"), function(x) slot(x, "ID"))

p.df <- data.frame( ID=1:length(M_intersect), row.names = pid)

M_intersect <- SpatialPolygonsDataFrame(M_intersect, p.df)
class(M_intersect)

# Croping raster --------------------------------------------------------------

library(maptools)

CropRaster<-function(filelist=NA,ShapeFile=NA)
{
  if(is.na(filelist)){
    filelist = choose.files(caption="Select ASCII files to crop: ")
  }
  if(is.na(ShapeFile)){
    ShapeFile = file.choose("Select shape file: ")
  }
  ext1 = readline("Enter Sufix to output file name: ")
  Shp1 = readShapePoly(ShapeFile)
  for (i in 1:length(filelist))
  {
    r1 = raster(filelist[i])
    cr1 = crop(r1,Shp1)
    cr2 = mask(cr1,Shp1)
    ### For file with extension.
    FileName = paste(substr(filelist[i],1,nchar(filelist[i])-4),ext1,substr(filelist[i],nchar(filelist[i])-3,nchar(filelist[i])),sep="")
    ### For files without extension
    writeRaster(cr2,FileName, "ascii")
    plot(cr2)
    print(i)
  }
}

CropRaster()


# n sei 

library(raster)
library(rgdal)
library(maptools)
library(dplyr)

# Creating an object for each bioclim

bio1 <- climate[[1]]
bio2 <- climate[[2]]
bio3 <- climate[[3]]
bio4 <- climate[[4]]
bio5 <- climate[[5]]
bio6 <- climate[[6]]
bio7 <- climate[[7]]
bio8 <- climate[[8]]
bio9 <- climate[[9]]
bio10 <- climate[[10]]
bio11 <- climate[[11]]
bio12 <- climate[[12]]
bio13 <- climate[[13]]
bio14 <- climate[[14]]
bio15 <- climate[[15]]
bio16 <- climate[[16]]
bio17 <- climate[[17]]
bio18 <- climate[[18]]
bio19 <- climate[[19]]

list(
  bio1,
  bio2,
  bio3,
  bio4,
  bio5,
  bio6,
  bio7,
  bio8,
  bio9,
  bio10,
  bio11,
  bio12,
  bio13,
  bio14,
  bio15,
  bio16,
  bio17,
  bio18,
  bio19) -> bio.clim


# crop() -> mask() -> writeRaster('asc')


shp1 <- readOGR("Shapefiles/M.shp")


crop(bio1, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio1.asc")

crop(bio2, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio2.asc")

crop(bio3, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio3.asc")

crop(bio4, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio4.asc")

crop(bio5, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio5.asc")

crop(bio6, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio6.asc")

crop(bio7, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio7.asc")

crop(bio8, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio8.asc")

crop(bio9, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio9.asc")

crop(bio10, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio10.asc")

crop(bio11, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio11.asc")

crop(bio12, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio12.asc")

crop(bio13, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio13.asc")

crop(bio14, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio14.asc")

crop(bio15, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio15.asc")

crop(bio16, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio16.asc")

crop(bio17, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio17.asc")

crop(bio18, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio18.asc")

crop(bio19, shp1) -> cr1
mask(cr1, shp1) -> cr2
writeRaster(cr2, "Bioclims/bio19.asc")
