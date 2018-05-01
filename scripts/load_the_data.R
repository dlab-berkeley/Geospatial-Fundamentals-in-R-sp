
#
# Load libraries
#
library(sp)     # spatial objects and methods
library(rgdal)  # read/write from file; manage CRSs
library(rgeos)  # geometric operations
library(tmap)   # mapping spatial objects

#######################################################################
# Load Part I data
#######################################################################
# Read in from CSV file
sfhomes <- read.csv('data/sf_properties_25ksample.csv', 
                    stringsAsFactors = FALSE)

# subset the data
sfhomes15 <- subset(sfhomes, as.numeric(SalesYear) == 2015)

# Make a copy
sfhomes15_sp <- sfhomes15  

# Make it spatial
coordinates(sfhomes15_sp) <- c('lon','lat')


# Read in the Bart data from CSV file
bart <- read.csv("./data/bart.csv")

# Read in the Landmarks data from CSV File
landmarks <- read.csv("./data/landmarks.csv")

## RGDAL
sfboundary <- readOGR(dsn="data",layer="sf_boundary")


# Transform the CRS
sfboundary_lonlat <- spTransform(sfboundary, CRS("+init=epsg:4326"))

sfboundary_lonlat <- readOGR(dsn="data",layer="sf_boundary_lonlat")

# Convert the DF to a SPDF
coordinates(landmarks) <-c("X","Y")
# Define the CRS with an EPSG code for Web mercator
proj4string(landmarks) <- CRS("+init=epsg:3857") 
# Transform the CRS to WGS84
landmarks_lonlat <- spTransform(landmarks, CRS("+init=epsg:4326"))

highways <- readOGR(dsn="data", layer="sf_highways")
class(highways)
proj4string(highways)
highways_lonlat <- spTransform(highways, CRS("+init=epsg:4326"))

#######################################################################
# Load Part II data
#######################################################################


