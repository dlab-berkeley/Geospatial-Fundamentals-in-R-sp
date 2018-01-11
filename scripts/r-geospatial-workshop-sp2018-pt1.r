## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = T)

## ---- eval=F-------------------------------------------------------------
## install.packages(
##   c("sp","rgdal","tmap","classInt","RColorBrewer",
##     "ggplot2","leaflet", "ggmap"), dependencies=TRUE
## )

## ---- eval=F-------------------------------------------------------------
## 
## install.packages(
##   c("sp","rgdal","tmap","classInt","RColorBrewer",
##     "ggplot2","leaflet", "ggmap"), dependencies=TRUE
## )
## 

## ---- echo=F-------------------------------------------------------------
sfhomes <- read.csv('data/sf_properties_25ksample.csv')

## ---- eval=F-------------------------------------------------------------
## sfhomes <- read.csv('data/sf_properties.csv')
## head(sfhomes,6)
## 

## ---- eval=F-------------------------------------------------------------
## class(sfhomes) # what type of data object?
## dim(sfhomes) # how many rows and columns
## str(sfhomes) # display the structure of the object
## head(sfhomes) # take a look at the first 10 records
## summary(sfhomes) # explore the range of values
## hist(sfhomes$totvalue)  # plot he range of values for the totvalue column

## ------------------------------------------------------------------------
plot(sfhomes$lon, sfhomes$lat) # using base plot function

## ------------------------------------------------------------------------
library(ggplot2)

ggplot() + geom_point(data=sfhomes, aes(lon,lat), col="red", size=1)

## ------------------------------------------------------------------------
ggplot() + geom_point(data=sfhomes, aes(lon,lat, col=totvalue))

## ------------------------------------------------------------------------
sfhomes15 <- subset(sfhomes, as.numeric(SalesYear) == 2015)

nrow(sfhomes15) # How many records?


## ------------------------------------------------------------------------
hist(sfhomes15$totvalue) # What is the distribution of totvalue?

## ------------------------------------------------------------------------
ggplot() + geom_point(data=sfhomes15, aes(x=lon, y=lat, col=totvalue))

## ------------------------------------------------------------------------
library(ggmap)

# fetch map data (default=Google) to plot
# ?get_map
sf_map <- get_map("San Francisco, CA")  

## ------------------------------------------------------------------------
ggmap(sf_map)

## ------------------------------------------------------------------------
ggmap(sf_map) +
  geom_point(data=sfhomes15, aes(x=lon, y=lat, col=totvalue))


## ------------------------------------------------------------------------

# Get the center point of the data
sf_ctr <- c(lon = mean(sfhomes15$lon), lat = mean(sfhomes15$lat))
sf_ctr  # take a look

# create the map
sf_basemap <- get_map(sf_ctr, zoom=12, scale=1)

## ------------------------------------------------------------------------
ggmap(sf_basemap) +
  geom_point(data=sfhomes15, aes(x=lon, y=lat, col=totvalue))

## ------------------------------------------------------------------------
sf_basemap_lite <- get_map(sf_ctr, zoom=12, scale=1, 
                            maptype = "toner-lite", source="stamen")

## ------------------------------------------------------------------------
ggmap(sf_basemap_lite) +
  geom_point(data=sfhomes15, aes(x=lon, y=lat, col=totvalue))


## ---- eval=F-------------------------------------------------------------
## # Let's look at last 5 years
## sfhomes2010_2015 <- subset(sfhomes, as.numeric(SalesYear) > 2009)
## 
## ggmap(sf_basemap_lite) +
##   geom_point(aes(lon, lat, col=totvalue), data = sfhomes2010_15 )  +
##   facet_wrap(~ SalesYear)
## 

## ---- echo=F-------------------------------------------------------------
# DO DONT SHOW
# Let's look at last 5 years
sfhomes2010_15 <- subset(sfhomes, as.numeric(SalesYear) > 2005)

## ------------------------------------------------------------------------
ggmap(sf_basemap_lite) +
  geom_point(aes(lon, lat, col=totvalue), data = sfhomes2010_15 )  +
  facet_wrap(~ SalesYear)


## ------------------------------------------------------------------------
sfhomes1995_1999 <- subset(sfhomes, (as.numeric(SalesYear) >= 1995) & (as.numeric(SalesYear) <= 1999))

ggmap(sf_basemap_lite) +
  geom_point(aes(lon, lat, col=totvalue), data = sfhomes1995_1999 )  +
  facet_wrap(~ SalesYear)

## ------------------------------------------------------------------------
library(sp)
getClass("Spatial") 

## ------------------------------------------------------------------------

# First make a copy of the data frame
sfhomes15_sp <- sfhomes15

coordinates(sfhomes15_sp) <- c('lon','lat') # ORDER MATTERS!!

class(sfhomes15_sp) # check it

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
## proj4strings(sfhomes15_sp)
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


## ---- eval=F-------------------------------------------------------------
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

# use an EPSG code
proj4string(sfhomes15_sp) <- CRS("+init=epsg:4326") 

# or enter the proj4 string
# proj4string(sfhomes15_sp) <- CRS("+proj=longlat 
#                               +ellps=WGS84 +datum=WGS84 +no_defs")  

## ------------------------------------------------------------------------
proj4string(sfhomes15_sp)

## ------------------------------------------------------------------------

proj4string(sfboundary)
proj4string(sfhomes15_sp)
proj4string(sfboundary) == proj4string(sfhomes15_sp)


## ------------------------------------------------------------------------

sfboundary_lonlat <- spTransform(sfboundary, CRS("+init=epsg:4326"))

# or
# sfboundary_lonlat <- spTransform(sfboundary, 
#                              CRS(proj4string(sfhomes15_sp)))


## ---- eval=F-------------------------------------------------------------
## proj4string(sfhomes15_sp) == proj4string(sfboundary_lonlat)
## 

## ---- eval=F-------------------------------------------------------------
## plot(sfboundary_lonlat)
## points(sfhomes15_sp, col="red")
## points(sfhomes15_sp[sfhomes15_sp$totvalue<1000000,], col="green")
## 

## ---- echo=F-------------------------------------------------------------
plot(sfboundary_lonlat)
points(sfhomes15_sp, col="red")
points(sfhomes15_sp[sfhomes15_sp$totvalue<1000000,], col="green")


## ---- eval=F-------------------------------------------------------------
## # write transformed data to a new shapefile
## writeOGR(sfboundary_lonlat,
##           dsn = "data",
##           layer = "sfbounary_lonlat",
##           driver="ESRI Shapefile")
## 
## # is it there?
## dir("data")

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

## ------------------------------------------------------------------------
qtm(sfboundary_lonlat, borders="black", fill="beige")

## ------------------------------------------------------------------------
tm_shape(sfboundary_lonlat) + tm_polygons(col="beige", border.col="black")

## ------------------------------------------------------------------------
tm_shape(sfhomes15_sp) + tm_dots(col="totvalue", size=.5)

## ------------------------------------------------------------------------
tmap_mode("plot")

## ------------------------------------------------------------------------

tm_shape(sfboundary_lonlat) + tm_polygons(col="black", border.col="grey") + 
   tm_shape(sfhomes15_sp) + tm_dots(col="totvalue", size=.5)

## ------------------------------------------------------------------------
us_states <- readOGR("./data", "us_states_pop")

## ------------------------------------------------------------------------
qtm(us_states)

## ------------------------------------------------------------------------
tm_shape(us_states) + tm_polygons(col="grey", border.col = "white")

## ---- echo=F-------------------------------------------------------------
tm_shape(us_states) + tm_polygons(col="grey", border.col = "white")

## ------------------------------------------------------------------------
tm_shape(us_states, projection="+init=epsg:5070") + tm_polygons(col="grey", border.col = "white")

## ------------------------------------------------------------------------
tm_shape(us_states, projection="+init=epsg:5070") + tm_polygons(col="grey", border.col = "white") +
tm_shape(us_states) + tm_borders(col="purple") 

## ------------------------------------------------------------------------
us_states_5070 <- spTransform(us_states, CRS("+init=epsg:5070"))


## ------------------------------------------------------------------------
tm_shape(us_states_5070) + tm_polygons(col="beige") +
  tm_shape(us_states) + tm_borders(col="purple")

## ---- eval=F-------------------------------------------------------------
## library(knitr)
## purl("r-geospatial-workshop-sp2018-pt1.Rmd", output = "scripts/r-geospatial-workshop-sp2018-pt1.r", documentation = 1)

###
# > sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] bindrcpp_0.2  ggmap_2.7     ggplot2_2.2.1
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.14      bindr_0.1         compiler_3.4.1    plyr_1.8.4        bitops_1.0-6      tools_3.4.1       digest_0.6.13     memoise_1.1.0    
# [9] evaluate_0.10.1   tibble_1.3.4      gtable_0.2.0      lattice_0.20-35   pkgconfig_2.0.1   png_0.1-7         rlang_0.1.4.9000  mapproj_1.2-5    
# [17] yaml_2.1.16       proto_1.0.0       withr_2.1.1       stringr_1.2.0     dplyr_0.7.4       knitr_1.17        devtools_1.13.4   RgoogleMaps_1.4.1
# [25] maps_3.2.0        rprojroot_1.3-1   grid_3.4.1        glue_1.2.0        R6_2.2.2          jpeg_0.1-8        rmarkdown_1.8     sp_1.2-5         
# [33] reshape2_1.4.3    magrittr_1.5      backports_1.1.2   scales_0.5.0      htmltools_0.3.6   assertthat_0.2.0  colorspace_1.3-2  geosphere_1.5-7  
# [41] stringi_1.1.6     lazyeval_0.2.1    munsell_0.4.3     rjson_0.2.15     
