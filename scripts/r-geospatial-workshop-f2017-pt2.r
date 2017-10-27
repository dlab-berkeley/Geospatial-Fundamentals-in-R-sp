## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=F-------------------------------------------------------------
## # Clean environment
## rm(list=ls())
## 
## # Set the working directory
## setwd("~/Documents/Dlab/dlab_workshops/r-geospatial-f2017")

## ---- eval=F-------------------------------------------------------------
## library(sp)
## library(rgdal)
## library(rgeos)
## library(tmap)
## library(RColorBrewer)
## library(dplyr)
## ## leaflet, ggplot2, ggmap, maptools, ??

## ---- echo=F-------------------------------------------------------------
library(sp)
library(rgdal)
library(rgeos)
library(tmap)
library(RColorBrewer)
library(ggplot2)
library(ggmap)

## ------------------------------------------------------------------------
parks <- read.csv("./data/sf_neighborhood_parks.csv", 
                  stringsAsFactors = FALSE)

## ---- eval=F-------------------------------------------------------------
## head(parks)
## class(parks)
## str(parks)
## summary(parks)

## ------------------------------------------------------------------------
?coordinates

coordinates(parks) <- c("lon","lat") # or ~lon+lat

## ---- eval=F-------------------------------------------------------------
## head(parks)
## class(parks)
## str(parks)
## summary(parks)

## # For more info

## ---- eval=F-------------------------------------------------------------
## parks@bbox
## head(parks@coords)
## head(parks@data)
## head(parks$ParkName)

## ---- eval=T-------------------------------------------------------------

plot(parks)


## ---- eval=F-------------------------------------------------------------
## 
## tmap_mode("view")
## 
## # Start with a quick tmap or qtm
## qtm(parks)
## 

## ---- echo=F-------------------------------------------------------------

tmap_mode("view")
qtm(parks)


## ------------------------------------------------------------------------
#library(rgdal)

# read the sf_nhoods.shp file from the current working directory
sf_nhoods <- readOGR(dsn="./data", layer="sf_nhoods", 
                     stringsAsFactors = F)


## ---- eval=F-------------------------------------------------------------
## class(sf_nhoods)
## summary(sf_nhoods)
## str(sf_nhoods@data)
## head(sf_nhoods@data)

## ------------------------------------------------------------------------

tm_shape(sf_nhoods) + tm_polygons(col="beige") + tm_shape(parks) + tm_symbols(col="green")


## ---- eval=F-------------------------------------------------------------
## parks_w_nhoods <- over(parks, sf_nhoods)

## ---- eval=F-------------------------------------------------------------
## 
## parks_w_nhoods <- over(parks, sf_nhoods)
## 

## ------------------------------------------------------------------------
# What is the CRS of the parks?
parks@proj4string # or proj4string(parks)

# What is the CRS of the sf_nhoods?
sf_nhoods@proj4string

## ------------------------------------------------------------------------
# Set the CRS for parks to be the same as that for sf_nhoods
proj4string(parks) <- CRS(proj4string(sf_nhoods))

# make sure the CRSs are the same
proj4string(parks) == proj4string(sf_nhoods) 

## ------------------------------------------------------------------------
# Now try the overlay operation again
# For each feature in parks give me info about the sf_nhoods 
# at the corresponding location
parks_w_nhoods <-over(parks,sf_nhoods)

## ---- eval=F-------------------------------------------------------------
## head(parks_w_nhoods) # take a look at the output
## class(parks_w_nhoods)
## nrow(parks_w_nhoods)
## nrow(parks)

## ------------------------------------------------------------------------
# Take a look at the data before we change it
# head(parks@data)

# Now combine the information about neighborhoods with the spatial data 
parks@data <- cbind(parks@data,parks_w_nhoods)  
  ## NOTE - binding to the parks@data not parks!!!

# Review and note the change
head(parks@data)

## ------------------------------------------------------------------------

tm_shape(sf_nhoods) + tm_polygons(col="beige") + tm_shape(parks) + 
  tm_symbols(col="green", size=0.04)

## ---- eval=F-------------------------------------------------------------
## # attribute query because we already made the spatial association
## parks[parks$nhood=='Noe Valley',]$ParkName  ## Why doesn't this work?
## 
## ## what about this syntax?
## parks[which(parks$nhood=='Noe Valley'),]$ParkName
## 

## ------------------------------------------------------------------------
# What park doesn't have a value for nhood?
parks[is.na(parks$nhood=='Noe Valley'),]$ParkName  

## ------------------------------------------------------------------------
parks@data[parks$ParkID==48,]$nhood <- "Bayview Hunters Point"

# What parks are in Noe Valley?
parks[parks$nhood=='Noe Valley',]$ParkName  ## Now it works work?

# How many Parks are in Noe Valley?
length(parks[parks$nhood=='Noe Valley',]$ParkName)

## ------------------------------------------------------------------------

nhood_park_count0.df <- aggregate(x = parks, by = sf_nhoods, 
                                  FUN = length)

nhood_park_count.df <- aggregate(x = parks["ParkID"], by = sf_nhoods, 
                                 FUN = length)

## ------------------------------------------------------------------------
# Some checks first
#head(nhood_park_count.df@data)
#nrow(sf_nhoods)
#nrow(nhood_park_count.df)

sf_nhoods$park_count <- nhood_park_count.df$ParkID

head(sf_nhoods@data)

## ---- eval=F-------------------------------------------------------------
## tm_shape(sf_nhoods) + tm_polygons(col="park_count", palette="Reds",
##                                   popup.vars=c("nhood","park_count"))

## ---- echo=F-------------------------------------------------------------
tm_shape(sf_nhoods) + tm_polygons(col="park_count", palette="Reds", popup.vars=c("nhood","park_count")) 

## ------------------------------------------------------------------------

# compute the park acreage within each neighborhood
nhood_park_acres.df <- aggregate(x = parks["Acreage"], 
                                 by = sf_nhoods, FUN = sum)

# Add the acreage sums back to the SF Neighborhoods layer
sf_nhoods$park_acreage <- nhood_park_acres.df$Acreage

## ------------------------------------------------------------------------
tm_shape(sf_nhoods) + tm_polygons(col="park_acreage", palette="Reds", 
                        popup.vars=c("nhood","park_count", "park_acreage"))

## ------------------------------------------------------------------------
parks_utm10 <- spTransform(parks, CRS("+init=epsg:26910"))
sfhoods_utm10 <- spTransform(sf_nhoods, CRS("+init=epsg:26910"))

## ---- eval=F-------------------------------------------------------------
## 
## gArea(sfhoods_utm10)
## 
## gArea(sfhoods_utm10, byid=TRUE)
## 

## ------------------------------------------------------------------------
# Add the output from gArea to sfhoods_utm10 as `nhood_area_m2`
sfhoods_utm10$nhood_area_m2 <- gArea(sfhoods_utm10, byid=T)

# Convert to acreage (1 sqmeter = 0.000247105 acres) as `nhood_acreage`
sfhoods_utm10$nhood_acreage <- sfhoods_utm10$nhood_area_m2 * 0.000247105

# Create a new column that is proportion of park_acreage / nhood_acreage 
sfhoods_utm10$prop_park_acres <- sfhoods_utm10$park_acreage / sfhoods_utm10$nhood_acreage

## ------------------------------------------------------------------------
tm_shape(sfhoods_utm10) + tm_polygons(col="park_acreage", palette="Reds",
                                popup.vars=c("nhood","park_count", "park_acreage", 
                                             "nhood_acreage", "prop_park_acres"))

## ------------------------------------------------------------------------
sftracts <- readOGR(dsn="./data", layer="sf_pop_by_tracts")
head(sftracts@data)
proj4string(sftracts)

## ------------------------------------------------------------------------
sftracts_utm10 <- spTransform(sftracts, CRS("+init=epsg:26910"))

head(sftracts_utm10)

## ------------------------------------------------------------------------
tm_shape(sftracts_utm10) + tm_polygons(border.col="blue", 
                                       col="pop14", palette="Reds") + 
  tm_shape(sfhoods_utm10) + tm_polygons(alpha=0, border.col="yellow", lwd=2)


## ------------------------------------------------------------------------

# Convert the SpatialPolygonsDataFrame to a SpatialPointsDataFrame
sftract_ctrs <- SpatialPointsDataFrame(sftracts_utm10, 
                                       data=sftracts_utm10@data)

# define the CRS
proj4string(sftract_ctrs) <- CRS(proj4string(sftracts_utm10))

## ------------------------------------------------------------------------
tm_shape(sftracts_utm10) + tm_polygons(border.col="blue", col="pop14",
                                       palette="Reds") +
  tm_shape(sfhoods_utm10) + tm_polygons(alpha=0, border.col="yellow", lwd=2) +
  tm_shape(sftract_ctrs) + tm_symbols(col="black", size=.5)

## ------------------------------------------------------------------------
# over
nhood_and_tract_over_output <- over(sftract_ctrs,sfhoods_utm10)

# We only want the nhood column
tract_nhood_only <-nhood_and_tract_over_output[c("nhood")]

# check it
nrow(tract_nhood_only) == nrow(sftract_ctrs)

# Add the neighborhood for each tract vector to the sftracts SPDF
sftracts_utm10@data <- cbind(sftracts_utm10@data,tract_nhood_only)

## ------------------------------------------------------------------------

pop_by_hood <- aggregate(pop14 ~ nhood, sftracts_utm10, sum)

# Take a look
#View(pop_by_hood)


## ------------------------------------------------------------------------
#global sum
sum(pop_by_hood$pop14)

#is it the same as the sum of pop14 in the sftracts_utm10 layer?
sum(sftracts_utm10$pop14)

## ------------------------------------------------------------------------

# Add the acreage sums back to the SF Neighborhoods layer
sfhoods_utm10$pop14 <- pop_by_hood$pop14

# compute the pop density
sfhoods_utm10$people_per_park <- sfhoods_utm10$pop14/ 
                          sfhoods_utm10$park_acreage

## ------------------------------------------------------------------------
tm_shape(sfhoods_utm10) + tm_polygons(col="people_per_park", 
          popup.vars=c("nhood","pop14","people_per_park",
                       "park_acreage","park_count")) 


## ---- eval=F-------------------------------------------------------------
## writeOGR(sfhoods_utm10, ".","sf_park_nhoods", driver="ESRI Shapefile")

## ------------------------------------------------------------------------
library(ggmap)
library(ggplot2)
coit_tower <- geocode('Coit Tower, San Francisco, CA')
coit_tower_pt <- c(coit_tower$lon, coit_tower$lat) 
coit_tower_pt


## ------------------------------------------------------------------------
dist_matrix <- spDistsN1(parks,coit_tower_pt, longlat=T)
parks$coit_distkm <- dist_matrix

#View(parks@data)

## ------------------------------------------------------------------------

coordinates(coit_tower) <- c("lon","lat")
proj4string(coit_tower)<- CRS("+init=epsg:4326")

coit_tower_utm <- spTransform(coit_tower, CRS("+init=epsg:26910"))
coit_km_buffer <- gBuffer(coit_tower_utm, width=1000)

#map check
#qtm(coit_tower_utm)
#qtm(coit_km_buffer)


## ------------------------------------------------------------------------

near_coit <- over(coit_km_buffer,parks_utm10, returnList = T)

# Take a look at output
head(near_coit)


## ---- eval=F-------------------------------------------------------------
## library(knitr)
## purl("r-geospatial-workshop-f2017-pt2.Rmd", output = "scripts/r-geospatial-workshop-f2017-pt2.r", documentation = 1)

