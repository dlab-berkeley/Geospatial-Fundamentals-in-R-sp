## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=F-------------------------------------------------------------
## install.packages(
##   c("sp","rgdal","tmap","classInt","RColorBrewer",
##     "ggplot2","leaflet", "ggmap"), dependencies=TRUE
## )

## ---- eval=F-------------------------------------------------------------
## install.packages(
##   c("sp","rgdal","tmap","classInt","RColorBrewer",
##     "ggplot2","leaflet", "ggmap"), dependencies=TRUE
## )

## ---- eval=F-------------------------------------------------------------
## sfhomes <- read.csv('data/sf_properties_25ksample.csv',
##                     stringsAsFactors = FALSE)
## 
## # Take a look at first 5 rows and a few of the columns
## sfhomes[1:5,c("YearBuilt","totvalue","AreaSquareFeet","Neighborhood",
##               "NumBedrooms")]
## 

## ---- eval=T-------------------------------------------------------------
sfhomes <- read.csv('data/sf_properties_25ksample.csv', 
                    stringsAsFactors = FALSE)

# Take a look at first 5 rows and a few of the columns
sfhomes[1:5,c("YearBuilt","totvalue","AreaSquareFeet","Neighborhood",
              "NumBedrooms")]


## ---- eval=F-------------------------------------------------------------
## class(sfhomes)            # what is the data object type?
## dim(sfhomes)              # how many rows and columns
## str(sfhomes)              # display the structure of the object
## head(sfhomes)             # take a look at the first 10 records
## summary(sfhomes)          # explore the range of values
## summary(sfhomes$totvalue) # explore the range of values for one column
## hist(sfhomes$totvalue)    # histogram for the totvalue column

## ---- eval=FALSE---------------------------------------------------------
## plot(sfhomes$lon, sfhomes$lat) # using base plot function

## ------------------------------------------------------------------------
plot(sfhomes$lon, sfhomes$lat) # using base plot function

## ------------------------------------------------------------------------
library(ggplot2)

## ---- eval=F-------------------------------------------------------------
## library(ggplot2)
## 
## ggplot() + geom_point(data=sfhomes, aes(lon,lat))
## 

## ------------------------------------------------------------------------
ggplot() + geom_point(data=sfhomes, aes(lon,lat), size=1)

## ------------------------------------------------------------------------
ggplot() + geom_point(data=sfhomes, aes(lon,lat), size=1) + coord_map()

## ---- eval=F-------------------------------------------------------------
## ggplot() + geom_point(data=sfhomes, aes(lon,lat, col=totvalue)) +
##   coord_map()

## ---- echo=F-------------------------------------------------------------
ggplot() + 
  geom_point(data=sfhomes, aes(lon,lat, col=totvalue), size=1) + 
  coord_map()

## ---- eval=F-------------------------------------------------------------
## sfhomes_low2high <- sfhomes[order(sfhomes$totvalue, decreasing = FALSE),]
## 
## ggplot() +
##   geom_point(data=sfhomes_low2high, aes(lon,lat, col=totvalue)) +
##   coord_map()

## ---- echo=F-------------------------------------------------------------
sfhomes_low2high <- sfhomes[order(sfhomes$totvalue, decreasing = FALSE),]

ggplot() + 
  geom_point(data=sfhomes_low2high, aes(lon,lat, col=totvalue)) + 
  coord_map()

## ------------------------------------------------------------------------
sfhomes_high2low <- sfhomes[order(sfhomes$totvalue, decreasing = T),]
ggplot() + geom_point(data=sfhomes_high2low, aes(lon,lat, col=totvalue)) + 
  coord_map()

## ---- eval=F-------------------------------------------------------------
## sfhomes2010_15 <- subset(sfhomes_low2high, as.numeric(SalesYear) > 2009)
## 
## ggplot() +
##   geom_point(aes(lon, lat, col=totvalue), data = sfhomes2010_15 )  +
##   facet_wrap(~ SalesYear)

## ---- echo=F-------------------------------------------------------------
sfhomes2010_15 <- subset(sfhomes_low2high, as.numeric(SalesYear) > 2009)

ggplot() +
  geom_point(aes(lon, lat, col=totvalue), data = sfhomes2010_15 )  +
  facet_wrap(~ SalesYear)

## ------------------------------------------------------------------------
library(ggmap)

## ---- eval=F-------------------------------------------------------------
## #register_google(key="AIzXXXXXXXXXXXXXXXXXxPSE") # your key here
## register_google(key="AIzaSyBE23PnI3bzQBJkhXuhXTenSfwwOVH6xUo")

## ---- eval=F-------------------------------------------------------------
## #devtools::install_github("dkahle/ggmap")
## #library(ggmap)

## ---- eval=F-------------------------------------------------------------
## sf_map <- get_map("San Francisco, CA")

## ------------------------------------------------------------------------
sf_map <- get_map("San Francisco, CA")  

## ---- eval=F-------------------------------------------------------------
## ggmap(sf_map)

## ------------------------------------------------------------------------
ggmap(sf_map)

## ---- eval=F-------------------------------------------------------------
## # ggplot() +
## 
## ggmap(sf_map) +
##   geom_point(data=sfhomes, aes(x=lon, y=lat, col=totvalue))
## 

## ------------------------------------------------------------------------
ggmap(sf_map) +
  geom_point(data=sfhomes, aes(x=lon, y=lat, col=totvalue))

## ---- eval=F-------------------------------------------------------------
## 
## # FIRST - subset the data
## sfhomes15 <- subset(sfhomes, as.numeric(SalesYear) == 2015)
## 
## # Get the center point of the data
## sf_ctr <- c(lon = mean(sfhomes$lon), lat = mean(sfhomes$lat))
## sf_ctr  # take a look
## 
## # create the map
## sf_basemap <- get_map(sf_ctr, zoom=12, scale=1)

## ------------------------------------------------------------------------
# FIRST - subset the data
sfhomes15 <- subset(sfhomes, as.numeric(SalesYear) == 2015)

# Get the center point of the data
sf_ctr <- c(lon = mean(sfhomes15$lon), lat = mean(sfhomes15$lat))
sf_ctr  # take a look

# create the map
sf_basemap <- get_map(sf_ctr, zoom=12, scale=1)

## ------------------------------------------------------------------------
ggmap(sf_basemap) +
  geom_point(data=sfhomes15, aes(x=lon, y=lat, col=totvalue))

## ------------------------------------------------------------------------
bart <- read.csv("./data/bart.csv")
# take a look
head (bart)

## ---- eval=F-------------------------------------------------------------
## ggmap(sf_basemap) +
##   geom_point(data=sfhomes15, aes(x=lon, y=lat))  +
##   geom_point(data=bart, aes(x=X,y=Y), col="red")

## ---- echo=F-------------------------------------------------------------
ggmap(sf_basemap) +
  geom_point(data=sfhomes15, aes(x=lon, y=lat))  +
  geom_point(data=bart, aes(x=X,y=Y), col="red")

## ---- eval=T-------------------------------------------------------------
police_stations <- read.csv("data/sf_police_addresses.csv", 
                            stringsAsFactors = F)
head(police_stations)

## ---- eval=F-------------------------------------------------------------
## #?geocode
## geocode(police_stations)

## ---- eval=T-------------------------------------------------------------
police_stations$full_address <- paste(police_stations$Address, 
                                      police_stations$City, police_stations$State)
head(police_stations)

## ------------------------------------------------------------------------
station_coords <- geocode(police_stations$full_address)

## ------------------------------------------------------------------------
police_stations <- cbind(police_stations, station_coords)
police_stations[1:5,c("PoliceDistrict","full_address","lon","lat")]

## ---- eval=F-------------------------------------------------------------
## ggmap(sf_basemap) +
##   geom_point(data=sfhomes15, aes(x=lon, y=lat))  +
##   geom_point(data=bart, aes(x=X,y=Y), col="red", size=3) +
##   geom_point(data=police_stations, aes(x=lon,y=lat), shape=22,
##              col="black", fill="grey", size=4)

## ---- echo=F-------------------------------------------------------------
ggmap(sf_basemap) +
  geom_point(data=sfhomes15, aes(x=lon, y=lat))  +
  geom_point(data=bart, aes(x=X,y=Y), col="red", size=3) +
  geom_point(data=police_stations, aes(x=lon,y=lat), shape=22, 
             col="black", fill="grey", size=4)

## ------------------------------------------------------------------------
landmarks <- read.csv("./data/landmarks.csv")
head(landmarks)

## ---- eval=F-------------------------------------------------------------
## ggmap(sf_basemap) +
##   geom_point(data=sfhomes15, aes(x=lon, y=lat))  +
##   geom_point(data=bart, aes(x=X,y=Y), col="red", size=3) +
##   geom_point(data=landmarks, aes(x=X,y=Y), shape=22,
##              col="black", fill="grey", size=4)

## ---- echo=F-------------------------------------------------------------
ggmap(sf_basemap) +
  geom_point(data=sfhomes15, aes(x=lon, y=lat))  +
  geom_point(data=bart, aes(x=X,y=Y), col="red", size=3) +
  geom_point(data=landmarks, aes(x=X,y=Y), shape=22, 
             col="black", fill="grey", size=4)

## ---- echo=F-------------------------------------------------------------
head(landmarks)

## ---- eval=F-------------------------------------------------------------
## library(sp)
## getClass("Spatial")

## ---- echo=F-------------------------------------------------------------
library(sp)
getClass("Spatial") 

## ---- eval=F-------------------------------------------------------------
## 
## sfhomes15_sp <- sfhomes15  # Make a copy - why?
## 
## class(sfhomes15_sp)  # check the class of the object
## 
## coordinates(sfhomes15_sp) <- c('lon','lat')  #  Make it spatial - ORDER MATTERS!!
## 
## class(sfhomes15_sp)  # check the class of the object

## ------------------------------------------------------------------------

sfhomes15_sp <- sfhomes15  # Make a copy - why?

class(sfhomes15_sp)  # check the class of the object

coordinates(sfhomes15_sp) <- c('lon','lat')  #  Make it spatial - ORDER MATTERS!!

class(sfhomes15_sp)  # check the class of the object

## ------------------------------------------------------------------------
coordinates(sfhomes15_sp)

## ------------------------------------------------------------------------
str(sfhomes15) # the data frame

## ------------------------------------------------------------------------
str(sfhomes15_sp) # the SPDF

## ---- eval=F-------------------------------------------------------------
## 
## summary(sfhomes15_sp)
## head(sfhomes15_sp@data)
## class(sfhomes15_sp@data)
## 
## sfhomes15_sp@bbox
## bbox(sfhomes15_sp)
## 
## head(sfhomes15_sp@coords)
## head(sfhomes15_sp$lat)
## head(sfhomes15_sp$lon)
## 
## sfhomes15_sp@proj4string
## proj4string(sfhomes15_sp)
## 

## ------------------------------------------------------------------------

proj4string(sfhomes15_sp) # get a CRS object


## ------------------------------------------------------------------------
plot(sfhomes15_sp)  # using sp::plot

## ------------------------------------------------------------------------

library(rgdal)

# See what file types are supported by rgdal drivers
# ogrDrivers()$name


## ------------------------------------------------------------------------

dir("data", pattern="sf_boundary")


## ------------------------------------------------------------------------

sfboundary <- readOGR(dsn="data",layer="sf_boundary")

# or
# sfboundary <- readOGR("data","sf_boundary")
# but not
#sfboundary <- readOGR(dsn="data/",layer="sf_boundary")


## ---- eval=F-------------------------------------------------------------
## class(sfboundary)
## str(sfboundary)
## head(sfboundary@data)

## ------------------------------------------------------------------------

plot(sfboundary)


## ---- eval=F, echo=F-----------------------------------------------------
## ## Demonstration - A more complex `SpatialPolygonsDataFrame`
## library(tigris)
## calcounties <- counties(state="California", cb=T)
## class(calcounties)
## sf_cen <- subset(calcounties, COUNTYFP == "075")
## plot(sf_cen)
## head(sf_cen@data)
## # Now look in Envi Window
## 

## ---- eval=F-------------------------------------------------------------
## 
## proj4string(sfboundary)
## 

## ------------------------------------------------------------------------

plot(sfboundary)
points(sfhomes15_sp, col="red")


## ---- eval=F-------------------------------------------------------------
## plot(sfboundary)
## points(sfhomes15_sp, col="red")

## ---- eval=F-------------------------------------------------------------
## 
## proj4string(sfboundary)
## proj4string(sfhomes15_sp)
## proj4string(sfboundary) == proj4string(sfhomes15_sp)
## 

## ------------------------------------------------------------------------

proj4string(sfboundary)
proj4string(sfhomes15_sp)
proj4string(sfboundary) == proj4string(sfhomes15_sp)


## ------------------------------------------------------------------------
sfboundary@bbox
sfhomes15_sp@bbox

## ------------------------------------------------------------------------
bbox(sfhomes15_sp)

## ------------------------------------------------------------------------

proj4string(sfhomes15_sp) <- CRS("+proj=longlat 
                               +ellps=WGS84 +datum=WGS84 +no_defs")  

## ---- eval=F-------------------------------------------------------------
## 
## # use an EPSG code for WGS84
## proj4string(sfhomes15_sp) <- CRS("+init=epsg:4326")
## 
## # or enter the parameter string
## # proj4string(sfhomes15_sp) <- CRS("+proj=longlat
## #                               +ellps=WGS84 +datum=WGS84 +no_defs")

## ------------------------------------------------------------------------
proj4string(sfhomes15_sp)

## ---- eval=F-------------------------------------------------------------
## # Define the CRS with the parameter string
## proj4string(sfhomes15_sp) <- CRS("+proj=longlat
##                               +ellps=WGS84 +datum=WGS84 +no_defs")

## ---- eval=F-------------------------------------------------------------
## # Define the CRS with an EPSG code for WGS84
## proj4string(sfhomes15_sp) <- CRS("+init=epsg:4326")

## ------------------------------------------------------------------------

proj4string(sfboundary)
proj4string(sfhomes15_sp)
proj4string(sfboundary) == proj4string(sfhomes15_sp)


## ------------------------------------------------------------------------

sfboundary_lonlat <- spTransform(sfboundary, CRS("+init=epsg:4326"))


## ------------------------------------------------------------------------
proj4string(sfhomes15_sp) == proj4string(sfboundary_lonlat)


## ---- eval=F-------------------------------------------------------------
## plot(sfboundary_lonlat)
## points(sfhomes15_sp, col="red")
## 
## 

## ---- echo=F-------------------------------------------------------------
plot(sfboundary_lonlat)
points(sfhomes15_sp, col="red")


## ---- eval=F-------------------------------------------------------------
## # write transformed data to a new shapefile
## writeOGR(sfboundary_lonlat,
##           dsn = "data",
##           layer = "sfbounary_lonlat",
##           driver="ESRI Shapefile")
## 
## # is it there?
## dir("data")

## ---- eval=F-------------------------------------------------------------
## # Convert the DF to a SPDF
## coordinates(landmarks) <-c("X","Y")
## # Define the CRS with an EPSG code for Web mercator
## proj4string(landmarks) <- CRS("+init=epsg:3857")
## # Transform the CRS to WGS84
## landmarks_lonlat <- spTransform(landmarks, CRS("+init=epsg:4326"))
## # map it
## plot(sfboundary_lonlat)
## points(sfhomes15_sp, col="red")
## points(landmarks_lonlat, col="green")

## ---- echo=F-------------------------------------------------------------
# Convert the DF to a SPDF
coordinates(landmarks) <-c("X","Y")
# Define the CRS with an EPSG code for Web mercator
proj4string(landmarks) <- CRS("+init=epsg:3857") 
# Transform the CRS to WGS84
landmarks_lonlat <- spTransform(landmarks, CRS("+init=epsg:4326"))
# map it
plot(sfboundary_lonlat)
points(sfhomes15_sp, col="red")
points(landmarks_lonlat, col="green")

## ---- eval=F-------------------------------------------------------------
## highways <- readOGR(dsn="data", layer="sf_highways")
## class(highways)
## proj4string(highways)
## highways_lonlat <- spTransform(highways, CRS("+init=epsg:4326"))
## 
## plot(sfboundary_lonlat)
## lines(highways_lonlat, col="black")
## points(sfhomes15_sp, col="red")
## points(landmarks_lonlat, col="green")
## 

## ---- echo=F-------------------------------------------------------------
highways <- readOGR(dsn="data", layer="sf_highways")
class(highways)
proj4string(highways)
highways_lonlat <- spTransform(highways, CRS("+init=epsg:4326"))

plot(sfboundary_lonlat)
lines(highways_lonlat, col="black")
points(sfhomes15_sp, col="red")
points(landmarks_lonlat, col="green")


## ------------------------------------------------------------------------
# map of the sfhomes data by totalvaue
spplot(sfhomes15_sp,"totvalue")

## ------------------------------------------------------------------------
library(tmap)


## ------------------------------------------------------------------------
qtm(sfhomes15_sp)

## ------------------------------------------------------------------------
tmap_mode("view")
qtm(sfhomes15_sp)

## ------------------------------------------------------------------------
tmap_mode("plot")

## ---- eval=F-------------------------------------------------------------
## 
## ttm()       # Toggle Tmap mode between "interactive" and "plot"
## 
## last_map()  # display the last map
## 

## ------------------------------------------------------------------------
tmap_mode("plot")

## ------------------------------------------------------------------------
tm_shape(sfboundary_lonlat) + 
  tm_polygons(col="beige", border.col="black")

## ------------------------------------------------------------------------
tm_shape(sfhomes15_sp) + 
  tm_dots(col="red", size=.25)

## ------------------------------------------------------------------------
tm_shape(highways_lonlat) + 
  tm_lines(col="black")

## ---- eval=F-------------------------------------------------------------
## tm_shape(sfhomes15_sp) +
##   tm_dots(col="totvalue", size=.25)    # column names must be quoted

## ------------------------------------------------------------------------
tm_shape(sfhomes15_sp) + 
  tm_dots(col="totvalue", size=.25)  # columnn names **must** be quoted!

## ---- eval=F-------------------------------------------------------------
## tm_shape(sfboundary_lonlat) +
##   tm_polygons(col="beige", border.col="black") +
## tm_shape(highways_lonlat) +
##   tm_lines(col="black") +
## tm_shape(sfhomes15_sp) +
##   tm_dots(col="totvalue", size=.25)
## 

## ---- echo=F-------------------------------------------------------------
tm_shape(sfboundary_lonlat) + 
  tm_polygons(col="beige", border.col="black") +
tm_shape(highways_lonlat) + 
  tm_lines(col="black") +
tm_shape(sfhomes15_sp) + 
  tm_dots(col="totvalue", size=.25) 
 

## ---- eval=F-------------------------------------------------------------
## tm_shape(sfboundary_lonlat) +
##   tm_polygons(col="beige", border.col="black") +
## tm_shape(highways_lonlat) +
##   tm_lines(col="black") +
## tm_shape(sfhomes15_sp) +
##   tm_dots(col="totvalue", size=.25, title = "San Francisco Property Values (2015)") +
## tm_layout(inner.margins=c(.05, .2, .15, .05)) # bottom, left, top, right

## ---- echo=F-------------------------------------------------------------
tm_shape(sfboundary_lonlat) + 
  tm_polygons(col="beige", border.col="black") +
tm_shape(highways_lonlat) + 
  tm_lines(col="black") +
tm_shape(sfhomes15_sp) + 
  tm_dots(col="totvalue", size=.25, title = "San Francisco Property Values (2015)") + 
tm_layout(inner.margins=c(.05, .2, .15, .05)) # bottom, left, top, right

## ---- eval=F-------------------------------------------------------------
## tm_shape(sfboundary_lonlat) +
##   tm_polygons(col="beige", border.col="black") +
## tm_shape(highways_lonlat) +
##   tm_lines(col="black") +
## tm_shape(sfhomes15_sp) +
##   tm_dots(col="totvalue", size=.25,
##           title = "San Francisco Property Values (2015)") +
## tm_shape(landmarks) +
##   tm_markers(col="black", size=.4, text="name",text.just = "left") +
## tm_layout(inner.margins=c(.05, .2, .15, .05)) # bottom, left, top, right

## ---- echo=F-------------------------------------------------------------
tm_shape(sfboundary_lonlat) + 
  tm_polygons(col="beige", border.col="black") +
tm_shape(highways_lonlat) + 
  tm_lines(col="black") +
tm_shape(sfhomes15_sp) + 
  tm_dots(col="totvalue", size=.25, title = "San Francisco Property Values (2015)") + 
tm_shape(landmarks) +
  tm_markers(col="black", size=.4, text="name",text.just = "left") +
tm_layout(inner.margins=c(.05, .2, .15, .05)) # bottom, left, top, right

## ---- eval=F-------------------------------------------------------------
## ttm()       # Toggle Tmap mode between "interactive" and "plot"
## last_map()  # display the last map

## ---- echo=F-------------------------------------------------------------
ttm()       # Toggle Tmap mode between "interactive" and "plot"
last_map()  # display the last map 

## ---- eval=F-------------------------------------------------------------
## tm_shape(sfboundary_lonlat) +
##   tm_polygons(col="beige", border.col="black") +
## tm_shape(highways_lonlat) +
##   tm_lines(col="black") +
## tm_shape(sfhomes15_sp) +
##   tm_dots(col="totvalue", size=.25,
##           title = "San Francisco Property Values (2015)",
##           popup.vars=c("SalesYear","totvalue","NumBedrooms",
##                        "NumBathrooms","AreaSquareFeet")) +
## tm_shape(landmarks) +
##   tm_markers(col="black", size=.4, text="name",text.just = "left") +
## tm_layout(inner.margins=c(.05, .2, .15, .05)) # bottom, left, top, right

## ---- echo=F-------------------------------------------------------------
tm_shape(sfboundary_lonlat) + 
  tm_polygons(col="beige", border.col="black") +
tm_shape(highways_lonlat) + 
  tm_lines(col="black") +
tm_shape(sfhomes15_sp) + 
  tm_dots(col="totvalue", size=.25, 
          title = "San Francisco Property Values (2015)",
          popup.vars=c("SalesYear","totvalue","NumBedrooms",
                       "NumBathrooms","AreaSquareFeet")) + 
tm_shape(landmarks) +
  tm_markers(col="black", size=.4, text="name",text.just = "left") +
tm_layout(inner.margins=c(.05, .2, .15, .05)) # bottom, left, top, right

## ------------------------------------------------------------------------
map1 <- last_map()
map1 # then display it

## ---- eval=F-------------------------------------------------------------
## save_tmap(map1, "sf_properties.png", height=6) # Static image file
## save_tmap(map1, "sf_properties.html") # interactive web map

## ---- echo=F-------------------------------------------------------------
save_tmap(map1, "sf_properties.png", height=6) # Static image file
save_tmap(map1, "sf_properties.html") # interactive web map

## ---- eval=F-------------------------------------------------------------
## tmap_mode("interactive")
## map1 #  in console

## ---- eval=F-------------------------------------------------------------
## 
## # write transformed data objects to a new shapefile
## coordinates(bart) <- c("X","Y")
## proj4string(bart) <- CRS(proj4string(sfboundary_lonlat))
## 
## save_layers <- c("sfhomes15_sp", "sfboundary_lonlat", "highways_lonlat", "bart", "landmarks")
## for (i in save_layers){
##   print(paste("Saving [",i,"] to file..."))
##   writeOGR(get(i),
##           dsn = "saved_data",
##           layer = i,
##           driver="ESRI Shapefile")
## }
## # is it there?
## dir("saved_data")
## 

## ---- eval=F, echo=F-----------------------------------------------------
## ## Output code to script
## library(knitr)
## purl("r-geospatial-workshop-pt1.Rmd",
##      output = "scripts/r-geospatial-workshop-pt1.R", documentation = 1)
## 

