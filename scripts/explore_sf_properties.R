library(ggplot2)
library(ggmap)
library(sp)
library(rgdal)
#rspatial tutorial code
workshop_dir <- "~/Documents/Dlab/workshops/2018/r-geospatial-workshop"
setwd(workshop_dir)
dir()

d1 <- read.csv("data/sf_properties.csv", strip.white = T, stringsAsFactors = F)


# ----------------
# plots
# ----------------

 

# plot all properties
ggplot(d1, aes(x=lon, y=lat)) + geom_point()

# plot properties sold in 2015
d2015 <- subset(d1, as.numeric(SalesYear) == 2015)
nrow(d2015)
hist(d1$totvalue)
ggplot(d2015, aes(x=lon, y=lat, col=totvalue)) + geom_point()

# facet plot properties sold 2010 - 2016
d2010_16 <- subset(d1, as.numeric(SalesYear) > 2005)
ggplot(d2010_16, aes(x=lon, y=lat, col=totvalue, alpha=1)) +
  geom_point() +
  facet_wrap(~ SalesYear)
# facet plot properties sold 1995 to 2000
d1995_99 <- subset(d1, (as.numeric(SalesYear) >= 1995) & (as.numeric(SalesYear) < 2000))

ggplot(d1995_99, aes(x=lon, y=lat, col=totvalue, alpha=1)) +
  geom_point() +
  facet_wrap(~ SalesYear)


sf_map <- get_map("San Francisco, CA")
ggmap(sf_map)

sf_ctr <- c(lon = mean(d1$lon), lat = mean(d1$lat))
sf_ctr
sf_basemap <- get_map(sf_ctr, zoom=12, scale=1)
ggmap(sf_basemap)


ggmap(sf_basemap) +
  geom_point(aes(lon, lat), data = d2015 )

ggmap(sf_basemap) +
  geom_point(aes(lon, lat, col=totvalue), data = d2015 )

ggmap(sf_basemap) +
  geom_point(aes(lon, lat, col=Neighborhood), data = d2015 )

# Change basemap
# ?get_map to see options
sf_basemap_lite <- get_map(sf_ctr, zoom=12, scale=1, maptype = "toner-lite", source="stamen")

ggmap(sf_basemap_lite) +
  geom_point(aes(lon, lat, col=totvalue), data = d2015 )

# Facet maps with ggmap
ggmap(sf_basemap_lite) +
  geom_point(aes(lon, lat, col=totvalue), data = d2010_16 )  +
  facet_wrap(~ SalesYear)


# Use base_layer argument to ggmap() to specify data and x, y mappings
ggmap(sf_basemap_lite,
      base_layer = ggplot(d2010_16, aes(lon, lat))) +
  geom_point(aes(color = totvalue)) +
  facet_wrap(~ SalesYear)

ggmap(sf_basemap_lite,
      base_layer = ggplot(d2015, aes(lon, lat))) +
  geom_point(aes(color = totvalue)) +
  facet_wrap(~ Neighborhood)

ggmap(sf_basemap_lite,
      base_layer = ggplot(d2015, aes(lon, lat))) +
  geom_point(aes(color = totvalue)) +
  facet_wrap(~ SupeDistrict)

# GEOCODING with GGMAP
# Add a point for Union Square
city_hall <- geocode("City Hall, San Francisco, CA")
city_hall$lat
city_hall$lon

ggmap(sf_basemap_lite,
      base_layer = ggplot(d2015, aes(lon, lat))) +
  geom_point(aes(color = totvalue)) +
  geom_point(aes(x=city_hall$lon, y=city_hall$lat), color="red", size=4)

# heat map of density of POINTS (not values)
ggmap(sf_basemap_lite, extent = "device") + 
  #geom_density2d(data = d2015, aes(x = lon, y = lat), size = 0.3) + 
  stat_density2d(data = d2015, aes(x = lon, y = lat, fill = ..level.., alpha = ..level..), size = 0.01, 
       bins = 16, geom = "polygon") + 
  scale_fill_gradient(low = "green", high = "red") + 
  scale_alpha(range = c(0, 0.3), guide = FALSE)

# See
# https://github.com/dkahle/ggmap

## Limits

Oops - want to add a shapefile

Ooops - want to add non-lat/lon data

Oops - what houes are within walking distance (.25 miles to BartAZ)
Ooops want to do spatial summaries

#Add streets
sf_streets <- readOGR(dsn='data/sf_highways.geojson', layer="OGRGeoJSON")
proj4string(sf_streets)
sf_streets_utm10 <-spTransform(sf_streets, CRS("+init=epsg:26910"))


sf_nhoods <- readOGR(dsn="./data", layer="sf_nhoods", 
                     stringsAsFactors = F)
proj4string(sf_nhoods)

plot(sf_streets)

bart<- readOGR("data","BART_Sta_13")
plot(bart, add=T)
head(bart@data)
bart_utm <- spTransform(bart,CRS("+init=epsg:26910"))
install.packages("broom")
library(broom)

bart_utm_df <- tidy(bart_utm)
#-----------------------------
# More stuff


# NOT SURE IF THE FOLLOWING IS OBSOLETE
# You've seen you can add layers to a ggmap() plot by adding geom_***() layers and
# specifying the data and mapping explicitly, but this approach has two big downsides:
# further layers also need to specify the data and mappings, and facetting won't work at all.
# Luckily ggmap() provides a way around these downsides: the base_layer argument.
# You can pass base_layer a normal ggplot() call that specifies the default data and mappings for all layers.
# By moving aes(x, y) and data from the initial geom_point() function to the ggplot() call within the ggmap() call,
# you can add facets, or extra layers, the usual ggplot2 way.


# Exploring sp slots
one <- countries_spdf@polygons[[169]]
str(one, max.level=2)

# typical
# read in readOGR, summary(), plot(), head()

# explore crs
# proj4string() on nyc_tracts and neighborhoods
proj4string(nyc_tracts)
proj4string(neighborhoods)
# coordinates() on nyc_tracts and neighborhoods
head(coordinates(nyc_tracts))
head(coordinates(neighborhoods))


# before merge:
#setdiff(nyc_income$tract, nyc_tracts$TRACTCE)

a<- c(1,2,4)
b<- c(2,6,4)
setdiff(a,b)
intersect(a,b)
a == b

# --------------------------------------------------------------------------------------
# Outline
# --------------------------------------------------------------------------------------
# ADD: why R for geo?
# already use it - cutting edge - reproducibility - free - shiny - rleaflet - github (pipeline/workflows)
# Plotting points with ggplot2
# pros & cons
# Spatial objects - pros & cons
# SP objects - points
# Sp objects - polygons
# mapping the points
# mapping polys
# crs
# color palettes
# classification

# geoprocessing
# subsetting
# aggregating
# merge by attribute
# dissolve
# point in poly - count, mean
# compute area
# map density, counts, rates

# Aspatial aggregate and Join
x<- aggregate(totvalue ~ Neighborhood, d2015, mean)
x
nrow(x) == length(unique(d1$Neighborhood))

# Spatial aggregate


# --------------------------------------------------------------------------------------

