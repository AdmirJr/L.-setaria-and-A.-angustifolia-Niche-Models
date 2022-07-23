# Filtragem ambiental dos dados

library(flexsdm)
library(raster)
library(dplyr)
library(terra)

## Baixando os dados ambientais

climate <- getData("worldclim", var="bio",res=2.5)

## Selecionando apenas a América do Sul

plot(climate[[1]])
drawExtent() -> e # Seleciona a América do Sul manualmente 
climate.focus <- crop(climate, e)
stack(climate.focus) -> climate.focus
plot(climate.focus)
writeRaster(climate.focus, 'SA_climate/SA_climate')

climate.focus <- rast(climate.focus)

# Filtrando os dados ambientalmente


dados.amb <- occfilt_env(data = dados.coord,
                         x = 'longitude',
                         y = 'latitude',
                         id = 'ID',
                         env_layer = climate.focus,
                         nbins = 5)

