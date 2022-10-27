# Filtragem ambiental de novo mas com .asc

library(flexsdm)
library(raster)
library(terra)

file_list <- list.files("Bioclims/",
                        pattern = ".asc", 
                        full.names = T) 
envs <- stack(file_list)
climate.focus <- rast(climate.focus)


amb.5.asc <- occfilt_env(data = dados.coord,
                          x = 'longitude',
                          y = 'latitude',
                          id = 'ID',
                          env_layer = climate.focus,
                          nbins = 5)

amb.10.asc <- occfilt_env(data = dados.coord,
                      x = 'longitude',
                      y = 'latitude',
                      id = 'ID',
                      env_layer = climate.focus,
                      nbins = 10)

amb.20.asc <- occfilt_env(data = dados.coord,
                      x = 'longitude',
                      y = 'latitude',
                      id = 'ID',
                      env_layer = climate.focus,
                      nbins = 20)

amb.30.asc <- occfilt_env(data = dados.coord,
                      x = 'longitude',
                      y = 'latitude',
                      id = 'ID',
                      env_layer = climate.focus,
                      nbins = 30)

amb.40.asc <- occfilt_env(data = dados.coord,
                      x = 'longitude',
                      y = 'latitude',
                      id = 'ID',
                      env_layer = climate.focus,
                      nbins = 40)

write.csv(amb.5.asc, file = 'Datasets/amb_filt_5bin.csv')
