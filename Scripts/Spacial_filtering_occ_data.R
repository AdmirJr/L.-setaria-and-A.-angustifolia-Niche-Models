# Filtragem espacial dos dados -------------------------------------------------

## Packages

library(spThin)
library(maps)

## Filtrando

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

# Carregando os vários .csv criados de uma só vez ------------------------------

temp <- list.files("Spatial_thinned_data", 
                   pattern = ".csv",
                   full.names = TRUE)

spatial.thinned <- lapply(temp, read.csv)

# Salvando os plots para consulta ----------------------------------------------



pdf('maps/map.pdf')

for (p in c(1:length(spatial.thinned))) {
  
  plot(latitude ~ longitude, spatial.thinned[[p]])
  map(add = T)
  title(main = p)
  
}

dev.off()

