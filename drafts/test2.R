#' ---
#' title: "r-geospatial-pt1-fall2017"
#' author: "Patty Frontiera"
#' date: "October 16, 2017"
#' output: ioslides_presentation
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = T)

#' 
#' ## Before we begin
#' 
#' 1. Download the workshop files
#' 
#' 2. Install any required libraries
#' 
#' 
#' ## Creating Maps in R
#' 
#' Who are you?
#' 
#' Why are you here?
#' 
#' 
#' ## Geographic Data
#' 
#' 
#' ## Geographic Data
#' 
#' Observations about locations on or near the surface of the Earth.
#' 
#' <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/ee/Anatone_WA.jpg/640px-Anatone_WA.jpg"></img>
#' 
#' That tells you some useful information about this place but it doesn't tell you where it is on the surface of the Earth.
#' 
#' ## Anatone?
#' 
#' <img width="800px" src="images/anatone_google.png"></img>
#' 
#' ## Geographic Coordinates
#' 
#' `46.130479, -117.134167`
#'   
#' <img width="600px" src="images/anatone_google.png"></img>
#' 
#' Which coordinate value is latitude?  which one is longitude?
#' 
#'   
#' ## Geospatial data
#' 
#' Geospatial data are geographic data that represent location geometrically with coordinates, such as latitude and longitude. 
#' 
#' These coordinates are referenced to specific locations on the Earth using a a coordinate reference system. 
#' 
#' <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/62/Latitude_and_Longitude_of_the_Earth.svg/640px-Latitude_and_Longitude_of_the_Earth.svg.png"></img>
#' 
#' 
#' 
#' 
#' ## Anatone, WA   
#' 
#' `46.130479, -117.134167`
#' 
#' Latitude is zero at the equator and decreases to -90 at the South Pole and increases to + 90 at the North pole.
#' 
#' Longitude is zero at the Prime Meridian (near London, England) and increases to -180 as you move west and to + 180 as you move east.
#' 
#' So, since latitude is never > 90 or less than < -90 the order must be latitude, longitude.
#' 
#' 
#' ## Spatial Data
#' 
#' Spatial data is a more generic term that is not just for geographic data. 
#' 
#' Spatial data are powerful because software can dynamically determine spatial metrics like area and length, characteristics like distance and direction, and relationships like inside and intersects from these data.
#' 
#' ## TODO
#' 
#' TODO Add more slides about vector & raster data
#' with examples
#' 
#' ## Geospatial Data in R
#' 
#' ## Geospatial Data in R
#' 
#' There are many approaches to and packages for working with geospatial data in R.
#' 
#' One approach is to keep it simple and store geospatial data in a data frame
#' 
#' ##
#' 
## ------------------------------------------------------------------------
cafes <- read.csv('data/cafes.csv')
head(cafes)


#' ## Base plots of points
#' 
## ------------------------------------------------------------------------
plot(cafes$long,cafes$lat)

#' 
#' ## Plotting points with ggplot2
## ------------------------------------------------------------------------
library(ggplot2)
ggplot() + geom_point(data=cafes, aes(long,lat))

#' ## Plots with ggmap
## ---- eval=F-------------------------------------------------------------
## library(ggmap)
## berkeley_map <- get_map(location="Berkeley, CA", zoom=14)
## ggmap(berkeley_map) + geom_point(data=cafes, aes(long,lat))

#' 
#' ## WAIT!
#' 
#' For a number of reasons, it's not always great idea to treat spatial data like non-spatial data.
#' 
#' ## Why?
#' 
#' ## Problem 1: Complexity
#' 
#' Complexity of spatial data representations
#' 
#' Plus the attribute data that describe each feature
#' 
#' ## For example
#' 
#' ADD some complex geometries
#' 
#' 
#' ## For example
#' ```
#' # Spatial Line
#' l1 = cbind(c(1,2,3),c(3,2,2))
#' myline1 = Line(l1)
#' myline2 = Lines(list(myline1), ID="a")
#' myline3 = SpatialLines(list(myline2))
#' myline4 <- SpatialLinesDataFrame(myline3,data.frame(c("Road 1")), match.ID = F)
#' plot(myline4, col = c("red"))
#' str(myline4)
#' 
#' ```
#' 
#' 
#' 
#' ##  Problem 2. The Earth
#' 
#' The software needs to know how to reference spatial data to the surface of the Earth.
#' 
#' - data about Earth coordinate reference systems
#' - spatial objects that include CRS definitions
#' - methods for using this information, eg transformations
#' 
#' 
#' ## `sp` package
#' 
#' ### The `SP` Package
#' 
#' The `SP` package is most commonly used to provide support for spatial data objects in R. 
#' Other R packages that do things with spatial data typically build on these SP objects.
#' 
#' 
## ------------------------------------------------------------------------
library(sp)
 

#' 
#' ## SP Objects
#' <style>
#'    th,td{
#'      padding:5px 5px 5px 5px;
#'    }
#' </style>
#' <table border=1>
#' <tbody>
#' <tr><th>Vector Data </th><th>SP Spatial Class </th><th>SP Spatial Class with Attributes ></th></tr>
#' <tr><td>Points</td><td>SpatialPoints</td><td>SpatialPointsDataFrame</td></tr>
#' <tr><td>Lines</td><td>SpatialLines</td><td>SpatialLinesDataFrame</td></tr>
#' <tr><td>Polygons</td><td>SpatialPolygons</td><td>SpatialPolygonsDataFrame</td></tr>
#' </tbody>
#' </table>
#' <br> 
#' 
#' 
#' ## Let's take a look
#' 
#' ## Read in point Data
#' 
## ------------------------------------------------------------------------
rentals <- read.csv('data/sf_airbnb_2bds.csv')
class(rentals)
dim(rentals)


#' 
#' ## Examine Data Structure
## ------------------------------------------------------------------------
str(rentals)

#' 
#' ## Examine Data Content
## ------------------------------------------------------------------------
head(rentals)

#' 
#' ## Visualize Data
## ------------------------------------------------------------------------
hist(rentals$price)

#' 
#' ## Process Data
## ------------------------------------------------------------------------
cheap <- subset(rentals, price < 301)
hist(cheap$price)

#' 
#' ## Process Data Some More
## ------------------------------------------------------------------------
cheap_good <- subset(cheap, review_scores_rating > 98)
hist(cheap$price)

#' 
#' 
#' ## Make it spatial
#' 
## ------------------------------------------------------------------------
library(sp)

#' 
#' ## Create a SpatialPointsDataFrame
#' 
#' Use the `sp::coordinates()` method
#' 
#' Requires a vector indicating the x,y columns
#' 
#' SP will create a SpatialPointsDataFrame from csv
#' 
#' ## SpatialPointsDataFrame (SPDF)
## ------------------------------------------------------------------------

#First make a copy
cheap_good_orig <- cheap_good

coordinates(cheap_good) <- c('longitude','latitude')
class(cheap_good)

#' 
#' ## Compare SPDF to DF
## ------------------------------------------------------------------------
str(cheap_good_orig)

#' 
#' ## SPDF
## ------------------------------------------------------------------------
str(cheap_good)

#' 
#' ## S*DF Slots
#' 
#' You can see from **str(cheap_good)** that a S*DF object is a collection of slots or components. The key ones are:
#' 
#' - `@data` data frame of attributes that describe each location
#' - `@coords` the coordinates for each location (or centroid if polygon)
#' - `@bbox` the min and max lon(x) and lat(y) coordinates tthat together define the minimum bounding box around the locations
#' - `@proj4string` the coordinate reference system defintion as a string
#' 
#' ## S*DF Slots
#'  
#' Review the output of each of these:
#' 
## ---- eval=F-------------------------------------------------------------
## head(cheap_good@coords)
## head(cheap_good@data)
## cheap_good@bbox
## cheap_good@proj4string

#' 
#' ## S*DF Slots
#'  
#' `@coords`
#' 
## ---- eval=F-------------------------------------------------------------
## head(cheap_good@coords)
## 
## #cheap_coords$latitude
## #cheap_coords$longitude

#' 
#' ## S*DF Slots
#'  
#' `@data`
#' 
## ---- eval=F-------------------------------------------------------------
## head(cheap_good@data)

#' 
#' ## S*DF Slots
#'  
#' `@bbox`
#' 
## ---- eval=F-------------------------------------------------------------
## 
## cheap_good@bbox
## 

#' 
#' ## S*DF Slots
#'  
#' `@proj4string`
#' 
## ---- eval=F-------------------------------------------------------------
## 
##  cheap_good@proj4string
## 

#' 
#' 
#' ## What's missing
#' 
#' Are all the columns that were present in the DF now in the SPDF?
#' 
#' Is there a slot without data?
#' 
#' 
#' ## What is the CRS of the data?
## ------------------------------------------------------------------------

cheap_good@proj4string # get a CRS object
# or 
proj4string(cheap_good) # method to get or set the CRS


#' 
#' ## Define and Assign a CRS
#' 
#' Assign (or set) a CRS to a spatial data object 
#' - locates it on the surface of the Earth
#' 
#' You need to know 
#' 
#' - the CRS for the data
#' - how to define the CRS object
#' - how to assign it to the spatial data
#' 
#' ## ## Define and Assign a CRS
#' 
#' Known as `defining a projection` in ArcGIS
#' 
#' **Defining a CRS != Transforming CRS**
#' 
#' We will get to transformations soon!
#' 
#' ## CRS Objects
#' 
## ------------------------------------------------------------------------

# Create a CRS object
WGS84_CRS <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 

# Set the CRS of the SPDF
proj4string(cheap_good) <- WGS84_CRS

# check it
cheap_good@proj4string

#' 
#' ## Another way
#' 
## ------------------------------------------------------------------------

proj4string(cheap_good) <- CRS("+proj=longlat 
                               +ellps=WGS84 +datum=WGS84 +no_defs")  


# or use the EPSG code
proj4string(cheap_good) <- CRS("+init=epsg:4326") 


#' 
#' ## Define and Assign - Incorrectly?
#' 
#' What happens if we assign the wrong CRS?
## ------------------------------------------------------------------------
proj4string(cheap_good) <- CRS("+init=epsg:26910") # The EPSG code for UTM10 NAD83

#' ## Define and Assign - Incorrectly?
#' 
#' What happens if we assign the wrong CRS?
#' 
#' Software doesn't care but you need to!
#' 
#' ## Finding CRS Codes
#' 
#' See [http://spatialreference.org/](http://spatialreference.org/)
#' 
#' Use this site to find EPSG codes and proj4 CRS strings
#' 
#' ## Common CRS Codes
#' 
#' * `4326` Geographic, WGS84
#' 
#' * `4269` Geographic, NAD83
#' 
#' * `3310` Projected, CA ALbers
#' 
#' * `26910` UTM Zone 10, NAD83 (Northern Cal)
#' 
#' ## Geographic vs Projected CRSs
#' 
#' Geographic CRSs
#' 
#' - longitude, latitude
#' - specify the shape of the earth ellipsoid
#' - identify the origin (equator & prime meridean)
#' - specify units, eg decimal degrees
#' - "attach" the CRS to the Earth
#'   - global: center of the earth
#'   - local: local point of tangency
#' 
#' ## Projected CRS
#' 
#' - Geographic CRS
#' - Map projection to transform to 2D
#' - Parameters to adjust origins
#' - units (typically meters)
#' - X, Y coordinates
#' 
#' ## CRS
#' 
#' 
#' ## Challenge
#' 
#' Use [http://spatialreference.org/](http://spatialreference.org/) to make an educated guess as to the CRS of these coordinates:
#' 
#' X = 549228.058249, Y = 4176578.444299
#' 
#' Strategy:
#' - review the bounding box coordinates for the CRSs referenced by the above codes.
#' 
#' ## Challenge 2
#' 
#' What are the Units for that CRS?
#' 
#' 
#' ## Mapping Spatial Objects
#' 
#' ## spplot
#' 
#' `sp` includes a plotting method `spplot`
#' 
#' You can use it to create great maps but it is very low level
#' 
#' which means, complex, non-intuitive syntax, long code
#' 
#' ## `spplot` the Data
## ------------------------------------------------------------------------
spplot(cheap_good)
 

#' 
#' ## `spplot` the Data
## ------------------------------------------------------------------------

spplot(cheap_good,"price")


#' 
#' ## Challenge
#' 
#' Use `spplot` to create data maps from some of the other columns in the @data slot
#' 
#' Getting help:
#' 
#' `?spplot`
#' 
#' 
#' ## Examples
#' ```
#' spplot(cheap_good,"bathrooms")
#' 
#' spplot(cheap_good,"accommodates")
#' 
#' spplot(cheap_good,"property_type")
#' 
#' spplot(cheap_good,"neighbourhood")
#' 
#' ```
#' 
#' ## What's special about spatial objects?
#' 
#' ## What do you think this code does?
#' 
#' Think about it
#' 
#' Try it and see
#' 
#' Check the help `?spDist`
#' 
## ---- eval=F-------------------------------------------------------------
## 
## coit_tower <- c("-122.405837,37.802032")
## 
## cheap_good$coit_dist <-
##   spDistsN1(cheap_good,c(-122.405837,37.802032), longlat = T)
## 
## head(cheap_good@data)
## 

#' 
#' ## RECAP
#' 
#' 
#' ## Polygons and Lines
#' 
#' Point data often in CSV files.
#' 
#' Line, Polygon and raster data are more often in a spatial data file format.
#' 
#' ## Spatial Data File formats
#' 
#' Vector points, lines & polygons:
#' 
#' * CSV
#' * [ESRI Shapefile](https://en.wikipedia.org/wiki/Shapefile)
#' 
#' Raster grids
#' 
#' * GeoJSON
#' * TIFF, JPEG
#' 
#' ## ESRI Shapefile
#' 
#' <img src="images/shapefile.png"></img>
#' 
#' ## Reading in Shapefile
#' 
#' There's an R package for that! 
#' 
#' ## `rgdal`
#' 
#' `rgdal` is an R port of the powerful and widely used [GDAL](http://gdal.org) library.
#' 
#' It is the most commonly used R library for importing and exporting spatial data. 
#' 
#' * `OGR`: for vector data: readOGR() and writeOGR()
#' * `GDAL` for raster data: readGDAL() and writeGDAL()
#' 
#' ## Reading in an ESRI Shapefile
#' 
## ------------------------------------------------------------------------
library(rgdal)
ogrDrivers()$name
sfboundary <- readOGR(dsn="data",layer="sf_boundary")
# or
# sfboundary <- readOGR("data","sf_boundary")
# but not
#sfboundary <- readOGR(dsn="data/",layer="sf_boundary")



#' 
#' ## Getting help
#' 
#' - gdal.org
#' - `?readOGR
#' -  For more info on working with `rgdal` to load different types of spatial data in R see this excellent [tutorial](http://zevross.com/blog/2016/01/13/tips-for-reading-spatial-files-into-r-with-rgdal/) by Zev Ross.
#' 
#' ## Make a quick plot to check the data
## ------------------------------------------------------------------------

plot(sfboundary)
points(cheap_good, col="red")
points(cheap_good[cheap_good$price<100,], col="green")


#' 
#' ## Take a look
#' 
## ------------------------------------------------------------------------
str(sfboundary)  # not proj4string

#' 
#' 
#' ## Take a look at the attribute data
#' 
## ------------------------------------------------------------------------
head(sfboundary@data)   

#' 
#' ## Take a look at the coordinate data
## ------------------------------------------------------------------------
sfboundary@bbox

#' 
#' ## Make a quick plot to check the data
## ------------------------------------------------------------------------

plot(sfboundary)
text(sfboundary@polygons[[1]]@labpt[1],sfboundary@polygons[[1]]@labpt[2], "San Francisco")


#' 
#' ## Plot with Rentals
## ------------------------------------------------------------------------
plot(sfboundary)
points(cheap_good, col="red")

#' 
#' ## What's Wrong?
## ------------------------------------------------------------------------
proj4string(sfboundary)
proj4string(cheap_good)
proj4string(sfboundary) == proj4string(cheap_good)

#' 
#' ## What do the coord data look like?
#' 
## ------------------------------------------------------------------------
sfboundary@bbox
cheap_good@bbox

#' 
#' 
#' ## Reprojecting the Data
#' 
#' Data transformations super common
#' Most common is projection transfromation, or `reprojection`
#' 
## ------------------------------------------------------------------------
sf_lonlat <- spTransform(sfboundary, WGS84_CRS)

#What other syntax can use to transform?

#' 
#' ## spTransform
#' ```
#' # Set it to the same as another data layer
#' sf_lonlat <- spTransform(sfboundary, CRS(proj4string(cheap_good)))
#' 
#' # Use CRS parameters
#' sf_lonlat <- spTransform(sfboundary, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#' 
#' # USE CRS code
#' sf_lonlat <- spTransform(sfboundary, CRS("+init=epsg:4326"))
#' 
#' # use a CRS object
#' wgs84 <- CRS(proj4string(cheap_good))
#' sf_lonlat <- spTransform(sfboundary, WS84_CRS)
#' 
#' ```
#' 
#' ## Did it work?
## ------------------------------------------------------------------------
plot(sf_lonlat)
points(cheap_good, col="red")
points(cheap_good[cheap_good$price<100,], col="green")

#' 
#' ## Adding Lines
## ------------------------------------------------------------------------
sf_streets <- readOGR(dsn='data/sf_highways.geojson', layer="OGRGeoJSON")

#' ## Take a look
## ------------------------------------------------------------------------
plot(sf_streets)


#' 
#' ## 
## ------------------------------------------------------------------------
plot(sf_lonlat)
points(cheap_good, col="red")
points(cheap_good[cheap_good$price<100,], col="green")
lines(sf_streets)

#' 
#' ## 
## ------------------------------------------------------------------------

#' 
#' ## What can we do with our spatial data?
## ------------------------------------------------------------------------
coit_tower <- c("-122.405837,37.802032") 
#?spDist
cheap_good$coit_dist <- spDistsN1(cheap_good,c(-122.405837,37.802032), longlat = T) # Try False - what do you get?
head(cheap_good@data)
summary(cheap_good$coit_dist)
hist(cheap_good$coit_dist)
spplot(cheap_good,"coit_dist")

#' 
#' ## Tmap
#' 
#' `Tmap` stands for Thematic maps
#' 
#' Create pretty good maps with less code than the alternatives
#' 
#' Syntax should be familar to ggplot2 users, but simpler
#' 
#' 
#' ## Tmap
## ------------------------------------------------------------------------
library(tmap)
?tmap

#' 
#' ## Simple tmap 
## ------------------------------------------------------------------------
tm_shape(sf_lonlat) + tm_polygons()
tm_shape(sfboundary) + tm_polygons(col="beige")


#' 
#' 
#' ## Adding layers 
## ------------------------------------------------------------------------
tm_shape(sfboundary) + tm_polygons(col="beige") + tm_shape(cheap_good) + tm_symbols(col="red")


#' 
#' ## Data Driven Maps
#' 
#' Thematic Maps / Data Maps
#' 
#' What do we mean?
#' 
#' ## Symbol Maps
## ------------------------------------------------------------------------
tm_shape(sfboundary) + 
  tm_polygons(col="beige") + 
tm_shape(cheap_good) + 
  tm_symbols(size="coit_dist", title.size="Distance to Coit Tower (KM)")

#' 
#' ## Color Maps
## ------------------------------------------------------------------------
tm_shape(sfboundary) + 
  tm_polygons(col="beige") + 
tm_shape(cheap_good) + 
  tm_symbols(col="coit_dist", title.col="Distance to Coit Tower (KM)")

#' 
#' 
#' 
#'   
#' ## Data Classification
## ------------------------------------------------------------------------

#TODO - add info

#' 
#' 
#' ## Graduate Symbol Maps
## ------------------------------------------------------------------------
tm_shape(sfboundary) + 
  tm_polygons(col="white") + 
tm_shape(cheap_good) + tm_symbols(col="price", style="quantile", size=.5, 
      palette="Reds", auto.palette.mapping=F, 
      border.alpha=0, alpha=0.75, title.col="Airbnb 2bd Price")

#' 
#' ## Graduate Symbol Maps
## ------------------------------------------------------------------------
tm_shape(sfboundary) + 
  tm_polygons(col="white") + 
  tm_shape(cheap_good) + 
  tm_symbols(col="price", style="jenks", palette="Reds", size=.5, 
        border.alpha=0, alpha=0.75, title.col="Airbnb 2bd Price")

#' 
#' 
#' ## Style - classification method
#' 
#' Try different values for Style
## ------------------------------------------------------------------------
tm_shape(sfboundary) + 
  tm_polygons(col="white") + 
tm_shape(cheap_good) + tm_symbols(col="price", style="jenks",size=.5, 
      palette="Reds", auto.palette.mapping=F, 
      border.alpha=0, alpha=0.75, title.col="Airbnb 2bd Price")

#'  
#' ```
#' 
#' ## Color Palettes
#' Try different values for palette
## ------------------------------------------------------------------------

 

#' 
#' 
#' ## Choropleth maps
#' 
#' Map areas by values
#' 
#' ## SF Census Tract Data
#' 
## ------------------------------------------------------------------------
sftracts <- readOGR("./data", "sf_pop_by_tracts")
plot(sftracts)

#' 
#' ## Fixins
## ------------------------------------------------------------------------
head(sftracts@data)
sftracts <- subset(sftracts, pop14 >1)
plot(sftracts)
#?writeOGR
#writeOGR(sftracts,"data", "sftracts_wpop", driver="ESRI Shapefile")

#' 
#' ## Repeat
## ------------------------------------------------------------------------
sftracts <- readOGR("./data", "sftracts_wpop")
plot(sftracts)
spplot(sftracts,"pop14")

#' 
#' ## Plot with Tmap
## ------------------------------------------------------------------------
tm_shape(sftracts) + tm_polygons()
 

#' 
#' ## 
## ------------------------------------------------------------------------

head(sftracts@data)

#' 
#' 
#' ## 
## ------------------------------------------------------------------------
tm_shape(sftracts) + tm_polygons("pop14", title="Population")

#' 
#' ## 
## ------------------------------------------------------------------------
sftracts$popdens <- (sftracts$pop14 / sftracts$ALAND)
tm_shape(sftracts) + tm_polygons("popdens", title="Population")

#' ## 
## ------------------------------------------------------------------------
sftracts$popdens <- sftracts$pop14 / (sftracts$ALAND / (1000*1000))
tm_shape(sftracts) + tm_polygons("popdens", title="Population per KM2")


#' ## Try a different class ification method (style)
#' 
#' What is the default style?
## ------------------------------------------------------------------------
tm_shape(sftracts) + tm_polygons("popdens", title="Population per KM2", style="jenks")

#' 
#' ## Multivariate
## ------------------------------------------------------------------------
bigmap <- tm_shape(sfboundary) + 
   tm_polygons(col="beige") + 
 tm_shape(cheap_good) + 
    tm_symbols(size="coit_dist", title.size="Distance to Coit Tower (KM)", col="price", title.col="Price", shape="property_type", title.shape="Property Type") +
   tm_layout( legend.bg.color="white",inner.margins=c(.05,.05, .15, .25), title="Airbnb 2 Bedroom Rentals, San Francisco Fall 2017", legend.position=c("right","center"))
 
 bigmap


#' 
#' ## Interactive
## ------------------------------------------------------------------------
?tmap_mode
tmap_mode('view')
bigmap


#' 
#' ## 
#' 
#' ## Multivariate
## ------------------------------------------------------------------------
tm_shape(sftracts) + tm_polygons("popdens", title="Population per KM2", style="jenks") +
tm_shape(cheap_good) + 
   tm_symbols(size="price", title.size="Price") +
  tm_layout( legend.bg.color="white",inner.margins=c(.05,.05, .15, .25), title="Airbnb 2 Bedroom Rentals, San Francisco Fall 2017", legend.position=c("right","center"))

 


#' 
#' ## leaflet
## ------------------------------------------------------------------------
tmap_mode('plot')
library(leaflet)
?leaflet

#' 
#' ## Simple leaflet map
#' 
## ------------------------------------------------------------------------
leaflet(cheap_good) %>% addTiles() %>%
    addCircleMarkers(data = cheap_good, radius = 5, stroke=F,
    color = "purple", fillOpacity = 0.75
  )

#' 
#' ## Color Palettes
#' ```
#' library(RColorBrewer)
#' 
#' display.brewer.all()
#' 
#' ```
#' ## Add palette
#' 
#' 
## ------------------------------------------------------------------------
pal <- colorQuantile("Reds",NULL,5)
leaflet(cheap_good) %>% addTiles() %>%
    addCircleMarkers(
      data = cheap_good,
      radius = 6,
      color = ~pal(price),
      stroke = F, 
      fillOpacity = 0.75
  )

#' 
#' # Add popup
## ------------------------------------------------------------------------
popup_content <- cheap_good$name
popup_content <- paste0(popup_content, "<br>Price per night: $", cheap_good$price)
popup_content <- paste0(popup_content, "<br><a href=",cheap_good$listing_url,">More info...</a>")


## ------------------------------------------------------------------------
 
leaflet(cheap_good) %>% addTiles() %>%
    addCircleMarkers(
      data = cheap_good,
      radius = 6,
      color = ~pal(price),
      stroke = F, 
      fillOpacity = 0.75,
      popup = popup_content)
  

#' 
#' 
#' ## References
#' - https://data.cdrc.ac.uk/tutorial/an-introduction-to-spatial-data-analysis-and-visualisation-in-r
#' - http://neondataskills.org/tutorial-series/vector-data-series/
#' - http://www.nickeubank.com/gis-in-r/
#' - http://www.rspatial.org/spatial/rst/3-vectordata.html
#' - https://dl.dropboxusercontent.com/u/9577903/broomspatial.pdf
#' - https://github.com/Robinlovelace/Creating-maps-in-R/raw/master/intro-spatial-rl.pdf
#' - https://rstudio.github.io/leaflet
#' - http://zevross.com/blog/2015/10/14/manipulating-and-mapping-us-census-data-in-r-using-the-acs-tigris-and-leaflet-packages-3/
#' - zev
#' - http://rstudio-pubs-static.s3.amazonaws.com/6577_3b66f8d8f4984fb2807e91224defa854.html
#' - https://cengel.github.io/rspatial/
#' - http://robinlovelace.net/geocompr/ (The future is SF)
