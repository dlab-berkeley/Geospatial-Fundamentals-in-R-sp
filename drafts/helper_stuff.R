setwd("~/Documents/Dlab/dlab_workshops/r-geospatial-f2017")
rentals <- read.csv('data/sf_airbnb_2bds.csv')
head(rentals)
hist(rentals$price)
cheap <- subset(rentals, price < 500)
hist(cheap$price)
plot(cheap)
library(sp)
library(tmap)
cheap_sp <- cheap #make a copy
coordinates(cheap_sp) <- c("longitude","latitude")
plot(cheap_sp)

str(cheap_sp)

x <- as.SpatialPolygonsDataFrame(cheap_sp@bbox)

str(x)
x
coit_tower <- c("-122.405837,37.802032") #37.802032, -122.405837
ct <- SpatialPoints(coordinates(coit_tower))

set.seed(1331)
pts = cbind(1:5, 1:5)
dimnames(pts)[[1]] = letters[1:5]
df = data.frame(a = 1:5)
row.names(df) = letters[5:1]
df
pts
data(package="tmap")
data(metro)
#data(land)
data(World)
tm_shape(cheap_sp) + tm_dots()

head(metro)
head(pop2010)
proj4string(cheap_sp) <- CRS("+proj=longlat + ellps=WGS84")

tm_shape(metro) + tm_bubbles("pop2010","red")
tm_shape(cheap_sp) + tm_bubbles("price", "red")
tm_shape(cheap_sp) + tm_bubbles("price", "red",size.lim = c(0, 400))
tm_shape(cheap_sp) + tm_dots("price",style="quantile", palette="-RdYlBu")

qtm(cheap_sp, fill="price", size=20)

tm_shape(metro) +
  tm_symbols(size = "pop2010", col="pop2010", shape="pop2010") +
  tm_layout(legend.outside = TRUE, legend.outside.position = "bottom", legend.stack = "horizontal")

cheap_sp<- subset(cheap_sp, review_scores_rating > 98)
cheap_sp <-subset(cheap_sp, price < 301)
tm_shape(cheap_sp) + tm_symbols(size="price", col="white", border.co="blue")
tm_shape(cheap_sp) + tm_symbols(size="price", col="white", border.co="blue", style="jenks")
tm_shape(cheap_sp) + tm_symbols(col="price", size=.5)
tm_shape(cheap_sp) + tm_symbols(col="price", size=.5, palette="-RdYlGn", style="kmeans",auto.palette.mapping=F, border.col="white")
tm_shape(cheap_sp) + tm_symbols(col="price", size=.5, palette="Reds", style="kmeans",auto.palette.mapping=F, border.col="white", title.col="Airbnb 2bd Price")
tm_shape(cheap_sp) + tm_symbols(col="price", size=.5, palette="Reds", style="kmeans",auto.palette.mapping=F, border.alpha=0, title.col="Airbnb 2bd Price")
tm_shape(cheap_sp) + tm_symbols(col="price", size=.5, palette="Reds", style="kmeans",auto.palette.mapping=F, border.alpha=0, alpha=0.75, title.col="Airbnb 2bd Price")

map1 <- tm_shape(cheap_sp) + tm_symbols(col="price", size=.5, palette="Reds", style="kmeans",auto.palette.mapping=F, border.alpha=0, alpha=0.75, title.col="Airbnb 2bd Price")

lmap1 <- tmap_leaflet(map1)

lmap1



x_coords <- c(0,1,2,3,3)
y_coords <- c(0,1,2,3,0)
plot(x_coords,y_coords, col="red") # show points
lines(x_coords, y_coords, col="green") # connect for lines
x_coords <- c(0,1,2,3,3,0) # close for polygons
y_coords <- c(0,1,2,3,0,0)
polygon(x_coords,y_coords, col="lightblue")



# from the sp vignette:
l1 = cbind(c(1,2,3),c(3,2,2))
rownames(l1) = letters[1:3]
l2 = cbind(c(1,2,3),c(1,1.5,1))
rownames(l2) = letters[1:3]
Sl1 = Line(l1)
Sl2 = Line(l2)
S1 = Lines(list(Sl1), ID="a")
S2 = Lines(list(Sl2), ID="b")
Sl = SpatialLines(list(S1,S2))
summary(Sl)
plot(Sl, col = c("red", "blue"))

 

# Spatial POlygon
Sr1 = Polygon(cbind(c(2,4,4,1,2),c(2,3,5,4,2)))
Srs1 = Polygons(list(Sr1), "s1")
SpP = SpatialPolygons(list(Srs1))
plot(SpP, col="blue")

# Spatial Line
l1 = cbind(c(1,2,3),c(3,2,2))
myline1 = Line(l1)
myline2 = Lines(list(myline1), ID="a")
myline3 = SpatialLines(list(myline2))
myline4 <- SpatialLinesDataFrame(myline3,data.frame(c("myline")), match.ID = F)
plot(myline4, col = c("red"))
str(myline4)

sftracts <- readOGR("./data", "sf_pop_by_tracts")
str(sftracts)
head(sftracts@data)
proj4string(sftracts)
sftracts <- subset(sftracts, pop14 >1)
plot(sftracts)
?writeOGR
writeOGR(sftracts,"data", "sftracts_wpop", driver="ESRI Shapefile")
sftracts <- readOGR("./data", "sftracts_wpop")
plot(sftracts)
spplot(sftracts,"pop14")
