## Model calibration

# Areas by buffering records (100 km buffer)
M_buffer <- buffer_area(dados.coord, longitude = "longitude", latitude = "latitude", 
                        buffer_distance = 50)

#####
# Areas using convex hulls (including 75 km buffer)
M_convex <- convex_area(dados.coord, longitude = "longitude", latitude = "latitude", 
                        buffer_distance = 50)

#####
# Areas by selecting polygons (including 25 km buffer)
M_ecorreg <- polygon_selection(dados.coord, longitude = "longitude", latitude = "latitude",
                               polygons = ecor, buffer_distance = 50)

M_intersect <- gIntersection(M_buffer, M_convex)
M_intersect <- gIntersection(M_intersect, M_ecorreg)

par(mfrow = c(2, 2), cex = 0.6, mar = rep(0.3, 4))
plot(M_buffer); points(dados.coord[, 2:3]); legend("topleft", legend = "Buffer", bty = "n")
plot(M_convex); points(dados.coord[, 2:3]); legend("topleft", legend = "Convex hull", bty = "n")
plot(M_ecorreg); points(dados.coord[, 2:3]); legend("topleft", legend = "Ecorregions", bty = "n")
plot(M_intersect); points(dados.coord[, 2:3]); legend("topleft", legend = "Intersection", bty = "n")

plot(M_intersect)
map(add = T)

  
## Model calibration
  
# Areas by buffering records (100 km buffer)
M_buffer.1 <- buffer_area(dados.coord, longitude = "longitude", latitude = "latitude", 
                          buffer_distance = 50)

#####
# Areas using convex hulls (including 75 km buffer)
M_convex.1 <- convex_area(dados.coord, longitude = "longitude", latitude = "latitude", 
                        buffer_distance = 75)

#####
# Areas by selecting polygons (including 25 km buffer)
M_ecorreg.1 <- polygon_selection(dados.coord, longitude = "longitude", latitude = "latitude",
                               polygons = ecor, buffer_distance = 25)

M_intersect.1 <- gIntersection(M_buffer.1, M_convex.1)
M_intersect.1 <- gIntersection(M_intersect.1, M_ecorreg.1)



plot(M_intersect.1)
map(add = T)
plot(M_intersect)
map(add = T)
