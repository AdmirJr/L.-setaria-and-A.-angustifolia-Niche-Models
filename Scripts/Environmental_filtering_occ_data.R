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


dados.amb <- occfilt_env(data = dados.coord,
                         x = 'longitude',
                         y = 'latitude',
                         id = 'ID',
                         env_layer = climate.focus,
                         nbins = 5)

