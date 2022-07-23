# Vamos baixar dados de ocorrência do pássaro Leptasthenura setaria, o grimpeiro.
## Para isso vamos procurar nas seguintes bases de dados de ocorrência:
  ### GBIF
  ### iDigBio
  ### SpeciesLink
  ### VertNet
  ### iNat

# Pacotes ----------------------------------------------------------------------

library(tidyverse)
library(rgbif)
library(ridigbio)
library(spocc)
library(readxl)

# Baixando os dados e organizando os dados -------------------------------------

"Leptasthenura setaria" -> sp # Para facilitar

dados.brutos <- data.frame() # Para receber os dados

## GBIF ------------------------------------------------------------------------

gbif.data <- occ_data(scientificName = sp,
         hasCoordinate = TRUE,
         limit = 100000)
  

gbif.data <- select(gbif.data$data, 
                    "scientificName", 
                    "decimalLatitude", 
                    "decimalLongitude") 
  
colnames(gbif.data) <- c("name", "latitude", "longitude")
gbif.data[,'prov'] = 'gbif'


rbind(dados.brutos, gbif.data) -> dados.brutos
rm(gbif.data)

## iDigBio ---------------------------------------------------------------------

idigbio.data <- idig_search_records(rq = list(scientificname = sp, 
                                        geopoint = list(type = "exists"))) %>%
                select(scientificname, geopoint.lat, geopoint.lon)

colnames(idigbio.data) <- c("name", "latitude", "longitude")
idigbio.data[,'prov'] = 'idigbio'

rbind(dados.brutos, idigbio.data) -> dados.brutos
rm(idigbio.data)

## VertNet e iNat --------------------------------------------------------------

databases <- c('inat', 'vertnet')

spocc.data <- occ(query = sp,
                  from = databases,
                  limit = 100000,
                  has_coords = TRUE) %>% 
              occ2df() %>% 
              select(name, latitude, longitude, prov)

rbind(dados.brutos, spocc.data) -> dados.brutos
rm(spocc.data)

## SpeciesLink -----------------------------------------------------------------

### Para o SpeciesLink, é muito mais fácil baixar e carregar manualmente os
### dados do que utilizar API's

specieslink.data <- read_excel(path = 'SpeciesLink/L_set_occ.xlsx')

specieslink.data <- specieslink.data[, c("scientificname", "latitude", "longitude")]
specieslink.data[, "prov"] = "specieslink"
colnames(specieslink.data) <- c("name", "latitude", "longitude", "prov")

rbind(dados.brutos, specieslink.data) -> dados.brutos
rm(specieslink.data)

## Adicionando uma coluna ID ---------------------------------------------------

dados.brutos <- tibble::rowid_to_column(dados.brutos, "ID")

## Salvando o dataset ----------------------------------------------------------

write.csv(dados.brutos, file = 'Datasets/Uncleaned data.csv')
