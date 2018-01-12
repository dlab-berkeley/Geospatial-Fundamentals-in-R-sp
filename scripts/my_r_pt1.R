setwd("~/Documents/Dlab/workshops/2018/r-geospatial-workshop")
sfhomes <- read.csv('data/sf_properties_25ksample.csv')  
head(sfhomes)
class(sfhomes)
dim(sfhomes)
str(sfhomes)
summary(sfhomes$totvalue)
plot(sfhomes$lon, sfhomes$lat)

library(ggplot2)

ggplot() + geom_point(data=sfhomes, aes(lon,lat), col="red", size=1)

sfhomes15 <- subset(sfhomes, as.numeric(SalesYear) == 2015)

nrow(sfhomes15) # How many records?

ggplot() + geom_point(data=sfhomes15, aes(lon,lat,col=totvalue), size=3)

library(ggmap)

sf_map <- get_map("San Francisco, CA")  
ggmap(sf_map)

###########################
library(sp)
getClass("Spatial")

# First make a copy of the data frame
sfhomes15_sp <- sfhomes15

coordinates(sfhomes15_sp) <- c('lon','lat') # ORDER MATTERS!!

class(sfhomes15_sp) # check it
coordinates(sfhomes15_sp)

##
summary(sfhomes15_sp)
head(sfhomes15_sp@data)
class(sfhomes15_sp@data)

sfhomes15_sp@bbox
bbox(sfhomes15_sp)

head(sfhomes15_sp@coords)
head(sfhomes15_sp$lat)
head(sfhomes15_sp$lon)

sfhomes15_sp@proj4string
proj4string(sfhomes15_sp)

dir("data", pattern="sf_boundary")
sfboundary <- readOGR(dsn="data",layer="sf_boundary")

class(sfboundary)
str(sfboundary) 
head(sfboundary@data)  
###
plot(sfboundary)
points(sfhomes15_sp, col="red")
proj4string(sfhomes15_sp) <- CRS("+init=epsg:4326") 
proj4string(sfhomes15_sp)

sfboundary_lonlat <- spTransform(sfboundary, CRS("+init=epsg:4326"))

plot(sfboundary_lonlat)
points(sfhomes15_sp, col="red")
points(sfhomes15_sp[sfhomes15_sp$totvalue<1000000,], col="green")

library(tmap)
qtm(sfhomes15_sp, title="SF Homes")

tmap_mode("view")
m<-qtm(sfhomes15_sp)

save_tmap(m, "map.html")
tmap_mode("plot")

qtm(sfboundary)
tm_shape(sfboundary_lonlat) + tm_polygons(col="beige", border.col="black")
tm_shape(sfboundary_lonlat) + tm_lines(col="black")
?`tmap-element`

map2<- tm_shape(sfboundary_lonlat) + tm_polygons(col="black", border.col="grey") + 
  tm_shape(sfhomes15_sp) + tm_dots(col="totvalue", size=.5, title="SF Homes")

save_tmap(map2, "map2.png")
qtm(sfhomes15)
