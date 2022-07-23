# Limpandos os dados

# Pacotes ----------------------------------------------------------------------

library(dplyr)
library(maps)
library(CoordinateCleaner)

# Removendo toda linha que possui NA -------------------------------------------

dados.brutos <- dados.brutos %>% 
  filter(!is.na(longitude)) %>% 
  filter(!is.na(latitude))

# Vizualizando em um mapa ------------------------------------------------------

plot(latitude ~ longitude, dados.brutos)
map(add = T)

# Limpando os dados baseado nas coordenadas ------------------------------------

## Convertendo as colunas 'latitude' e 'longitude' de 'chr' para 'dbl'

dados.brutos[,c('latitude', 'longitude')] <- sapply(dados.brutos[,c('latitude','longitude')], as.numeric)

## Utilizando a função clean_coordinates() para limpar as coordenadas

dados.coord <- clean_coordinates(dados.brutos,
                                 lon = 'longitude',
                                 lat = 'latitude',
                                 species = 'name',
                                 value = 'clean')

## Vizualizando os dados

plot(latitude ~ longitude, dados.coord)
map(add = T)

## Retirando manualmente um ponto improvável

dados.coord <- subset(dados.coord, (latitude < -20))

## Vizualizando os dados

plot(latitude ~ longitude, dados.coord)
map(add = T)

# Renomeando a coluna 'name' ---------------------------------------------------

## A checagem manual confirma que "Leptasthenura setaria" é o nome atual da espécie

dados.coord[, 'name'] = sp

# Removendo duplicatas ---------------------------------------------------------

dados.coord <- dados.coord %>% 
                    distinct(latitude, longitude, .keep_all = TRUE)

## Vizualizando os dados

plot(latitude ~ longitude, dados.coord)
map(add = T)

# Salvando os dados em um .csv

write.csv(dados.coord, file = 'Datasets/Cleaned coordinates.csv')
