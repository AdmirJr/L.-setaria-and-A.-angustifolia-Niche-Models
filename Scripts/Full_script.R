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

dados.amb <- occfilt_env(data = dados.coord,
                         x = 'longitude',
                         y = 'latitude',
                         id = 'ID',
                         env_layer = climate.focus,
                         nbins = 5)

