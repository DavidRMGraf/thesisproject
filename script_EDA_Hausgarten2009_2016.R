## EDA picoplankton data from AWI Hausgarten stations (Arctic) provided by Katja Metfies
## Doro Hodapp, March 2019

rm(list=ls(all=TRUE)) # clear variables from previous

# load packages
library(vegan)
library(dplyr)
library(ggmap)
library(plotrix)
library(sp)
#options(stringsAsFactors = FALSE)


# set working directory 
setwd("~/Studium/19SS/BA/diffmaps_BA")

# load data 
dat <- read.csv("Sequences_Hausgarten2009-2016.csv", header=F)
dim(dat)
dat[1:5,1:10]

# convert data to wide format
data <- as.data.frame(t(dat[,])) # transpose data to have OTUs as column and sampleID as row
colnames(data) <- as.character(unlist(data[1,])) # use OTUs as column names
data <- data[-1,] # remove first row (original headers in csv file)
data <- data %>% rename(Station = Proben_ID_intern) # rename Proben ID for merging with coordinate file

# load station coordinates
coord <- read.csv("Koordinaten_Arctic_NorthSea_all_cor.csv")
head(coord)
coord <- coord[,c("Cruise","Station","lat.calc","long.calc")]
coord$Station <- as.character(coord$Station)

# merge information on community composition and station coordinates
otu_coord <- merge(data,coord, by = ("Station"))
dim(otu_coord)


# plot station location 
div_coord <- otu_coord[,c("Station","sample_ID","lat.calc","long.calc")]

div_coord_sp <- div_coord[,c("lat.calc","long.calc")] # select longitude and latitude columns
div_coord_sp <- na.omit(div_coord_sp) # remove missing values
coordinates(div_coord_sp) <- ~long.calc+lat.calc # spatial coordinate information 



# create map

library(rgdal)                                                                                                      
library(raster)
library(ggplot2)


# get world map
wm <- rworldmap::getMap(resolution ="high")
wm <- crop(wm, extent(0, 30, 70, 90)) # adjust extent of map to sampling region

# simple plot
plot(wm, axes=T) 
points(coordinates(div_coord_sp), pch=16, col="red") # add sampling stations to map
#plot(div_coord_sp, pch=16, col="red", add=T) # alternativeadd sampling stations to map
